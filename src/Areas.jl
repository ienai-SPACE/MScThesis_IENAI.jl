include("MTalgorithm.jl")
include("Origins.jl")
include("Convex.jl")
include("NonConvex.jl")
include("materialInputs.jl")


"""
    analyze_areas(geometry::AbstractGeometry, viewpoint::Viewpoint)

# Inputs
- `geometry::AbstractGeometry`
- `viewpoint::Viewpoint`
# Outputs
- `areas_convex(geometry, viewpoint)`
- `areas_nonconvex(geometry, viewpoint)`
"""
function analyze_areas(geometry::AbstractGeometry, viewpoint::Viewpoint)
    if is_convex(geometry)
        return areas_convex(geometry, viewpoint)
    else
        return areas_nonconvex(geometry, viewpoint)
    end
end

"""
    areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint)

# Inputs
-`geometry::AbstractGeometry`
-`viewpoint::Viewpoint`
# Outputs
- `Aproj`
- `Aref`
- `intercept_info::Vector{InteractionGeometry{Float64}}`                  : areas and incidence angles of intercepted triangles, ordered in increasing index
- `_face_normals::Vector{StaticArraysCore.SVector{3, Float64}}`           : normal vector of all intercepted triangles
- `culling_ratio`
"""
function areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint)
    #view.direction --> dir (particle beam direction)
    samplerG = GridFilter(1e5)
    samplerF = FibonacciSampler(1e5)
    samplerMC = MonteCarloSampler(1e5)
    sampler = samplerF

    #Solar cells identification + ray-tracing
    solarCellsGeo = getSolarCellsGeo(geometry, viewpoint, sampler)

    #PCSA: back-face culling + ray-tracing
    rti_vec, filtered_geometry, Aray = raytrace(geometry, viewpoint, sampler) #culling +  ray tracing
    valid_rti = rti_vec |> Filter(rti -> rti.mode == FrontFaceIntersection) |> tcollect
    Aproj = Aray * length(valid_rti)

    faces_hit_idx_nonunique = sort(valid_rti .|> rti -> rti.face_index)
    hit_idx = unique(faces_hit_idx_nonunique) |> Filter(idx -> idx > 0) |> collect
    ray_per_index = [count(x -> x == ii, faces_hit_idx_nonunique) for ii in hit_idx]
    # face_areas = [face_area(filtered_geometry, idx) for idx in hit_idx]
    # _face_vertices = [face_vertices(filtered_geometry, idx) for idx in hit_idx]
    _face_normals = [face_normal(filtered_geometry, idx) for idx in hit_idx]
    _face_angles = angle_V_n(viewpoint.direction, _face_normals)
    _face_areas = [Aray * ray_per_index[ii] / cos(_face_angles[ii]) for ii in 1:length(hit_idx)]
    Aref = sum([_face_areas[ii] for ii in 1:length(hit_idx)])
    _barycenters = getBarycenters(filtered_geometry, hit_idx)

    if is_homogeneous(geometry)
        intercept_info = map(ii -> InteractionGeometryHomo(_face_areas[ii], _face_angles[ii]), 1:lastindex(hit_idx))

    else
        intercept_info = map(ii -> InteractionGeometryHetero(
                _face_areas[ii],
                _face_angles[ii],
                filtered_geometry.faces[hit_idx[ii]].properties.m_srf
            ), eachindex(hit_idx))
    end

    culling_ratio = n_faces(filtered_geometry) / n_faces(geometry)

    return Aproj, Aref, intercept_info, _face_normals, culling_ratio, _barycenters, solarCellsGeo
end


"""
    areas_convex(geometry::AbstractGeometry, viewpoint::Viewpoint)

# Inputs
-`geometry::AbstractGeometry`
-`viewpoint::Viewpoint`
# Outputs
- `Aproj`
- `Aref`
- `intercept_info::Vector{InteractionGeometry{Float64}}`                  : areas and angles of intercepted triangles, ordered in increasing index
- `_face_normals::Vector{StaticArraysCore.SVector{3, Float64}}`           : normal vector of all intercepted triangles
- `culling_ratio`
"""
function areas_convex(geometry::AbstractGeometry, viewpoint::Viewpoint)

    #solar cells geometry calculations
    solarCell_geometry = getSolarCellsGeo(geometry::AbstractGeometry)

    #back-face culling
    filtered_geometry = filter_backfaces(geometry, viewpoint)
    #Number of triangles
    Ntri = n_faces(filtered_geometry)

    _face_areas = [face_area(filtered_geometry, idx) for idx in 1:Ntri]
    _face_normals = [face_normal(filtered_geometry, idx) for idx in 1:Ntri]
    _face_angles = angle_V_n(viewpoint.direction, _face_normals)
    _barycenters = getBarycenters(filtered_geometry, 1:Ntri)

    #Calculate areas
    Aproj = sum(projection(filtered_geometry, viewpoint, Ntri))
    Aref = sum([_face_areas[ii] for ii in 1:Ntri])

    if is_homogeneous(geometry)
        intercept_info = map(ii -> InteractionGeometryHomo(_face_areas[ii], _face_angles[ii]), 1:Ntri)

    else
        intercept_info = map(ii -> InteractionGeometryHetero(
                _face_areas[ii],
                _face_angles[ii],
                filtered_geometry.faces[ii].properties.m_srf
            ), 1:Ntri)
    end

    culling_ratio = n_faces(filtered_geometry) / n_faces(geometry)

    return Aproj, Aref, intercept_info, _face_normals, culling_ratio, _barycenters, solarCell_geometry
end

"""
    raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)

Perform back-face culling and ray-tracing

# Inputs
- `geometry::AbstractGeometry`      : storage of information about the 3D CAD model
- `viewpoint::Viewpoint`
- `sampler`                         : identification of Fibonacci, Monte Carlo or Grid Filter sampling methods
# Outputs
- `rt_vec`                          : distance from origin of local ref.frame to the perpendicular source plane
- `filtered_geometry`               : only faces that have been intercepted by the rays
- `Aray`                            : area of each ray [m^2]
"""
function raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    Ntri = n_faces(geometry)
    Ntri_preculling = Ntri

    #Preculling: Meshed satellite
    println("Ntri_preculling=", Ntri_preculling)
    println("max in MeshVertices (preculling)=", viewpoint.rmax)
    println("viewpoint.direction=", viewpoint.direction)

    #back-face culling
    filtered_geometry = filter_backfaces(geometry, viewpoint)
    if typeof(filtered_geometry) == Bool
        return 0, 0, 0
    else

        #Number of triangles
        Ntri = n_faces(filtered_geometry)

        new_viewpoint = shrink_viewpoint(filtered_geometry, viewpoint)
        println("culling ratio =", Ntri / Ntri_preculling)

        rt_vec, Aray = _raytrace(filtered_geometry, new_viewpoint, sampler)

        return rt_vec, filtered_geometry, Aray
    end
end

"""
    _raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)

This is where the actual ray-tracer performs its tasks

# Inputs
- `filtered_geometry::AbstractGeometry`      : storage of information about the 3D CAD model
- `new_viewpoint::Viewpoint`
- `sampler`                         : identification of Fibonacci, Monte Carlo or Grid Filter sampling methods
# Outputs
- `rt_vec`                          : distance from origin of local ref.frame to the perpendicular source plane
- `Aray`                            : area of each ray [m^2]
"""
function _raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    dir = viewpoint.direction
    rmax = viewpoint.rmax
    O, Norig, Aray = generate_ray_origins(sampler, dir, rmax, viewpoint.distance)        #coordinates of ray origins, number of origins

    Ntri = n_faces(geometry)
    println("Ntri in _raytrace=", Ntri)

    m = Map(jj -> begin
        orig = O[jj]
        ray = Ray(SV3(orig), dir)
        rti = ray_mesh_intersection(geometry, ray)
    end)

    rti_vec = 1:Norig |> m |> tcollect

    return rti_vec, Aray
end

export analyze_areas