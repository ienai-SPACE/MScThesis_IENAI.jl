include("MTalgorithm.jl")
include("Origins.jl")
include("Convex.jl")
include("NonConvex.jl")
include("materialInputs.jl")

"""
    OutGeometry{T}

- `area::Vector{T}`
- `angle::Vector{T}`
"""

struct OutGeometry{T}
    area::Vector{T}
    angle::Vector{T}
    normals::Vector{SVector{3,Float64}}
end


_eltype(::OutGeometry{T}) where {T} = T


homogeneity_dict = Dict("homogeneous" => 1, "heterogeneous" => 2)

""" 
    areas(rmax, distance, dir, triangles, convexFlag)

Obtaining all the areas and perpendicular angles of the triangular mesh elements that have been intercepted by the rays 
        - It has two inner functions `areasConvex(vertices, dir)` and `areasConcave(dir, rmax, distance, triangles, Ntri)`, for convex and non-convex shapes, respectively
#INPUT:
- `rmax`            : radius of the circular plane from where rays originate
- `distance`        : distance at which the circular plane is located (it should be out from the satellite body)
- `dir`             : direction of the oncoming particle
- `triangles`       : vertices coordinates of the triangular/quad mesh element
- `convexFlag`      :1 for convex, 0 for non-convex
#OUTPUT:
- `OutLMNTs:: OutGeometry{T}`                                   : struct of 3 fields (`area`, `angle`, and `normals`) each containing a vector with the respective magnitude of all intercepted surfaces
- `InteractionGeometry_v::Vector{InteractionGeometry{T}}`       : vector of struct storing the areas and angles of all intercepted surfaces       
- `Aproj`                                                       : projection of the intercepted triangular areas onto the selected direction 
- `Aref`                                                        : sum of all intercepted triangular areas
"""
function areas(rmax, distance, dir, triangles, convexFlag)
    if convexFlag == 1
        areas_convex(dir, triangles)
    else
        areas_nonconvex(rmax, distance, dir, triangles)
    end
end

"""
    areas_convex(Vdir, triangles)

For convex shapes: calculation of reference and projected areas, and normals and areas of each forward-facing element

#INPUTS
- `Vdir`             : direction of the oncoming particle
- `triangles`       : vertices coordinates of the triangular/quad mesh element
#OUTPUTS
- `OutLMNTs:: OutGeometry{T}`                                   : struct of 3 fields (`area`, `angle`, and `normals`) each containing a vector with the respective magnitude of all intercepted surfaces
- `InteractionGeometry_v::Vector{InteractionGeometry{T}}`       : vector of struct storing the areas and angles of all intercepted surfaces       
- `Aproj`                                                       : projection of the intercepted triangular areas onto the selected direction 
- `Aref`                                                        : sum of all intercepted triangular areas
"""

function areas_convex(Vdir, triangles)
    #Number of triangles
    Ntri = size(triangles, 1)

    counter = 0
    for jj ∈ 1:Ntri #for loop to iterate over all triangles/quads
        vertices = triangles[jj, 1:9]
        OutAreaConvex = areasConvex(vertices, Vdir)
        if OutAreaConvex[1] != 0.0
            counter += 1
            if counter == 1
                OutFacets = @SVector [jj, OutAreaConvex[1], OutAreaConvex[2], OutAreaConvex[3], OutAreaConvex[4], OutAreaConvex[5]]
            else
                OutFacets = hcat(OutFacets, [jj, OutAreaConvex[1], OutAreaConvex[2], OutAreaConvex[3], OutAreaConvex[4], OutAreaConvex[5]])
            end
        end
    end


    #store normal vector components into vector
    normal_v = Vector{SVector{3,Float64}}(undef, size(OutFacets, 2))
    for ii ∈ 1:lastindex(OutFacets[2, :])
        normal_v[ii] = SVector(OutFacets[4, ii], OutFacets[5, ii], OutFacets[6, ii])
    end


    #sum of all intercepted triangular areas
    Aref = sum(OutFacets[2, :])

    #Project all intercepted triangular areas and find projected total area

    #------pre-allocation-------------------
    Aproj = zeros(size(OutFacets, 2), 1)
    #---------------------------------------

    for ii ∈ 1:Int(size(OutFacets, 2))
        Aproj[ii] = abs(OutFacets[2, ii] * cos(OutFacets[3, ii]))
    end
    Aproj = sum(Aproj)                 #sum of all intercepted triangular projected areas

    OutLMNTs = OutGeometry(OutFacets[2, :], OutFacets[3, :], normal_v[:])
    T = _eltype(OutLMNTs)

    #pre-allocation of the vector to be populated by structs
    InteractionGeometry_v = Vector{InteractionGeometry{T}}(undef, length(OutFacets[2, :]))

    for ii ∈ 1:length(OutFacets[2, :])
        InteractionGeometry_v[ii] = InteractionGeometry(OutFacets[2, ii], OutFacets[3, ii])
    end

    return Aproj, Aref, OutLMNTs, InteractionGeometry_v
end

"""
    areas_nonconvex(rmax, distance, dir, triangles)

For convex shapes: calculation of reference and projected areas, and normals and areas of each forward-facing element

#INPUTS
- `rmax`            : radius of the circular plane from where rays originate
- `distance`        : distance at which the circular plane is located (it should be out from the satellite body)
- `dir`             : direction of the oncoming particle
- `triangles`       : vertices coordinates of the triangular/quad mesh element
#OUTPUTS
- `OutLMNTs:: OutGeometry{T}`                                   : struct of 3 fields (`area`, `angle`, and `normals`) each containing a vector with the respective magnitude of all intercepted surfaces
- `InteractionGeometry_v::Vector{InteractionGeometry{T}}`       : vector of struct storing the areas and angles of all intercepted surfaces       
- `Aproj`                                                       : projection of the intercepted triangular areas onto the selected direction 
- `Aref`                                                        : sum of all intercepted triangular areas
"""

function areas_nonconvex(rmax, distance, dir, triangles)

    #Number of triangles
    Ntri = size(triangles, 1)

    OutFacets = areasConcave(dir, rmax, distance, triangles, Ntri)


    #store normal vector components into vector
    normal_v = Vector{SVector{3,Float64}}(undef, size(OutFacets, 2))
    for ii ∈ 1:lastindex(OutFacets[2, :])
        normal_v[ii] = SVector(OutFacets[4, ii], OutFacets[5, ii], OutFacets[6, ii])
    end


    #sum of all intercepted triangular areas
    Aref = sum(OutFacets[2, :])

    #Project all intercepted triangular areas and find projected total area

    #------pre-allocation-------------------
    Aproj = zeros(size(OutFacets, 2), 1)
    #---------------------------------------

    for ii ∈ 1:Int(size(OutFacets, 2))
        Aproj[ii] = abs(OutFacets[2, ii] * cos(OutFacets[3, ii]))
    end
    Aproj = sum(Aproj)                 #sum of all intercepted triangular projected areas

    OutLMNTs = OutGeometry(OutFacets[2, :], OutFacets[3, :], normal_v[:])
    T = _eltype(OutLMNTs)

    #pre-allocation of the vector to be populated by structs
    InteractionGeometry_v = Vector{InteractionGeometry{T}}(undef, length(OutFacets[2, :]))

    for ii ∈ 1:length(OutFacets[2, :])
        InteractionGeometry_v[ii] = InteractionGeometry(OutFacets[2, ii], OutFacets[3, ii])
    end


    #int_geos = [InteractionGeometry(OutFacets[2, ii], OutFacets[3, ii]) for ii ∈ 1:length(OutFacets[2, :])]
    #int_geos = map(ii -> InteractionGeometry(OutFacets[2, ii], OutFacets[3, ii]), 1:length(OutFacets[2, :]))

    #print(OutLMNTs)
    return Aproj, Aref, OutLMNTs, InteractionGeometry_v
end

export OutGeometry


"""
    analyze_areas(geometry::AbstractGeometry, viewpoint::Viewpoint)

#INPUTS
-`geometry::AbstractGeometry`
-`viewpoint::Viewpoint`
#OUTPUTS
- `areas_convex(geometry, viewpoint)`
- `areas_nonconvex(geometry, viewpoint)`
"""

function analyze_areas(geometry::AbstractGeometry, viewpoint::Viewpoint, materialDistribution::String)
    if is_convex(geometry)
        return areas_convex(geometry, viewpoint, materialDistribution::String)
    else
        return areas_nonconvex(geometry, viewpoint, materialDistribution::String)
    end
end

"""
    areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint)

#INPUTS
-`geometry::AbstractGeometry`
-`viewpoint::Viewpoint`
#OUTPUTS
- `Aproj`
- `Aref`
- `intercept_info::Vector{InteractionGeometry{Float64}}`                  : areas and incidence angles of intercepted triangles, ordered in increasing index
- `_face_normals::Vector{StaticArraysCore.SVector{3, Float64}}`           : normal vector of all intercepted triangles
- `culling_ratio`
"""

function areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint, matDist::String)
    #view.direction --> dir (particle beam direction)
    samplerG = GridFilter(1e5)
    samplerF = FibonacciSampler(1e5)
    samplerMC = MonteCarloSampler(1e5)
    sampler = samplerF

    rti_vec, w, filtered_geometry, Aray = raytrace(geometry, viewpoint, sampler) #culling +  ray tracing
    valid_rti = rti_vec |> Filter(rti -> rti.mode == FrontFaceIntersection) |> tcollect
    Aproj = w * length(valid_rti)

    faces_hit_idx_nonunique = sort(valid_rti .|> rti -> rti.face_index)
    hit_idx = unique(faces_hit_idx_nonunique) |> Filter(idx -> idx > 0) |> collect

    ray_per_index = [count(x -> x == ii, faces_hit_idx_nonunique) for ii in hit_idx]
    face_areas = [face_area(filtered_geometry, idx) for idx in hit_idx]
    # _face_vertices = [face_vertices(filtered_geometry, idx) for idx in hit_idx]
    _face_normals = [face_normal(filtered_geometry, idx) for idx in hit_idx]
    _face_indices = [face_idx(filtered_geometry, idx) for idx in hit_idx]
    _face_angles = angle_V_n(viewpoint.direction, _face_normals)
    _face_areas = [Aray * ray_per_index[ii] / cos(_face_angles[ii]) for ii in 1:length(hit_idx)]
    Aref = sum([_face_areas[ii] for ii in 1:length(hit_idx)])

    if homogeneity_dict[matDist] == 1 #homogeneous --> single material
        intercept_info = map(ii -> InteractionGeometryHomo(_face_areas[ii], _face_angles[ii]), 1:lastindex(hit_idx))

        #TODO: CHECK THAT _face_indices IS DOING THE JOB RIGHT
    elseif homogeneity_dict[matDist] == 2 #heterogeneous --> more than one material 
        mesh_facets_material, materials = loadMaterialProperties()
        intercept_info = map(ii -> InteractionGeometryHetero(_face_areas[ii], _face_angles[ii], materials[mesh_facets_material[_face_indices[ii]]["material_index"]+1]["atomic mass"]), 1:lastindex(hit_idx))
    end

    culling_ratio = n_faces(filtered_geometry) / n_faces(geometry)

    return Aproj, Aref, intercept_info, _face_normals, culling_ratio, filtered_geometry, _face_indices, hit_idx
end


"""
    areas_convex(geometry::AbstractGeometry, viewpoint::Viewpoint)

#INPUTS
-`geometry::AbstractGeometry`
-`viewpoint::Viewpoint`
#OUTPUTS
- `Aproj`
- `Aref`
- `intercept_info::Vector{InteractionGeometry{Float64}}`                  : areas and angles of intercepted triangles, ordered in increasing index
- `_face_normals::Vector{StaticArraysCore.SVector{3, Float64}}`           : normal vector of all intercepted triangles
"""
function areas_convex(geometry::AbstractGeometry, viewpoint::Viewpoint, matDist::String)
    #back-face culling
    filtered_geometry = filter_backfaces(geometry, viewpoint)
    #Number of triangles
    Ntri = n_faces(filtered_geometry)
    #Calculate areas
    Aproj = sum(projection(filtered_geometry, viewpoint, Ntri))
    Atot = sum([filtered_geometry.faces[idx].area for idx in 1:Ntri])

    _face_areas = [face_area(filtered_geometry, idx) for idx in 1:Ntri]
    _face_normals = [face_normal(filtered_geometry, idx) for idx in 1:Ntri]
    _face_indices = [face_idx(filtered_geometry, idx) for idx in 1:Ntri]
    _face_angles = angle_V_n(viewpoint.direction, _face_normals)

    # intercept_info = map(ii -> InteractionGeometry(_face_areas[ii], _face_angles[ii]), 1:Ntri)
    if homogeneity_dict[matDist] == 1 #homogeneous --> single material
        intercept_info = map(ii -> InteractionGeometryHomo(_face_areas[ii], _face_angles[ii]), 1:Ntri)

        #TODO: CHECK THAT _face_indices IS DOING THE JOB RIGHT
    elseif homogeneity_dict[matDist] == 2 #heterogeneous --> more than one material 
        mesh_facets_material, materials = loadMaterialProperties()
        intercept_info = map(ii -> InteractionGeometryHetero(_face_areas[ii], _face_angles[ii], materials[mesh_facets_material[_face_indices[ii]]["material_index"]+1]["atomic mass"]), 1:Ntri)
    end

    return Aproj, Atot, intercept_info, _face_normals
end

"""
    raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)

Perform back-face culling and ray-tracing

#INPUTS
- `geometry::AbstractGeometry`      : storage of information about the 3D CAD model
- `viewpoint::Viewpoint`
- `sampler`                         : identification of Fibonacci, Monte Carlo or Grid Filter sampling methods
#OUTPUTS
- `rt_vec`                          : distance from origin of local ref.frame to the perpendicular source plane
- `w`                               : weight = (pi * rmax^2) / (number of rays)
- `filtered_geometry`               : only faces that have been intercepted by the rays
- `Aray`                            : area of each ray [m^2]
"""

function raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    Ntri = n_faces(geometry)
    Ntri_preculling = Ntri

    #Preculling: Meshed satellite
    # _face_vertices_preCulling = [face_vertices(geometry, idx) for idx in 1:Ntri]
    # _face_normals_preCulling = [face_normal(geometry, idx) for idx in 1:Ntri]
    println("Ntri_preculling=", Ntri_preculling)
    println("max in MeshVertices (preculling)=", viewpoint.rmax)
    println("viewpoint.direction=", viewpoint.direction)

    #back-face culling
    filtered_geometry = filter_backfaces(geometry, viewpoint)

    #Number of triangles
    Ntri = n_faces(filtered_geometry)
    # _face_vertices_preRT = [face_vertices(filtered_geometry, idx) for idx in 1:Ntri]
    # _face_normals_preRT = [face_normal(filtered_geometry, idx) for idx in 1:Ntri]

    new_viewpoint = shrink_viewpoint(filtered_geometry, viewpoint)
    println("culling ratio =", Ntri / Ntri_preculling)

    rt_vec, w, Aray = _raytrace(filtered_geometry, new_viewpoint, sampler)

    return rt_vec, w, filtered_geometry, Aray
end

function _raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    dir = viewpoint.direction
    rmax = viewpoint.rmax
    O, Norig, Aray = generate_ray_origins(sampler, dir, rmax, viewpoint.distance)        #coordinates of ray origins, number of origins

    Ntri = n_faces(geometry)
    println("Ntri in _raytrace=", Ntri)
    # _face_vertices_preRT2 = [face_vertices(geometry, idx) for idx in 1:Ntri]
    # _face_normals_preRT2 = [face_normal(geometry, idx) for idx in 1:Ntri]

    m = Map(jj -> begin
        orig = O[jj]
        ray = Ray(SV3(orig), dir)
        rti = ray_mesh_intersection(geometry, ray)
    end)

    rti_vec = 1:Norig |> m |> tcollect

    return rti_vec, (pi * rmax^2) / (Norig), Aray
end

export analyze_areas