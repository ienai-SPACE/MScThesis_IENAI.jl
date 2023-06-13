# using LinearAlgebra
# using StaticArrays

include("MTalgorithm.jl")
include("Origins.jl")
include("Convex.jl")
include("NonConvex.jl")

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
        areas_convex(rmax, distance, dir, triangles)
    else
        areas_nonconvex(rmax, distance, dir, triangles)
    end
end

function areas_convex(rmax, distance, Vdir, triangles)  #this is actually dir
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


    #int_geos = [InteractionGeometry(OutFacets[2, ii], OutFacets[3, ii]) for ii ∈ 1:length(OutFacets[2, :])]
    #int_geos = map(ii -> InteractionGeometry(OutFacets[2, ii], OutFacets[3, ii]), 1:length(OutFacets[2, :]))

    #print(OutLMNTs)
    return Aproj, Aref, OutLMNTs, InteractionGeometry_v
end

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



function analyze_areas(geometry::AbstractGeometry, viewpoint::Viewpoint)
    if is_convex(geometry)
        return nothing #areas_convex(viewpoint, geometry)
    else
        return areas_nonconvex(geometry, viewpoint)
    end
end

function areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint)
    #view.direction --> dir (particle beam direction)
    samplerG = GridFilter(5e5)
    samplerF = FibonacciSampler(1e5)
    samplerMC = MonteCarloSampler(1e5)
    sampler = samplerF
    # rti_vec, w, _face_vertices_preRT, _face_normals_preRT, _face_vertices_preCulling, _face_normals_preCulling, indices, filtered_geometry = raytrace(geometry, viewpoint, sampler) #culling +  ray tracing
    rti_vec, w, filtered_geometry = raytrace(geometry, viewpoint, sampler) #culling +  ray tracing
    valid_rti = rti_vec |> Filter(rti -> rti.mode == FrontFaceIntersection) |> tcollect
    println("length(valid_rti)=", length(valid_rti))
    Aproj = w * length(valid_rti) # / sampler.ray_density
    faces_hit_idx_nonunique = valid_rti .|> rti -> rti.face_index
    faces_hit_idx = unique(faces_hit_idx_nonunique) |> Filter(idx -> idx > 0) |> collect
    sorted_hit_idx = sort(faces_hit_idx)
    # println("length(filter(!iszero, indices))=", length(filter(!iszero, indices)))
    # indices = Int.(unique(filter(!iszero, indices)))


    Aref = sum(face_area(filtered_geometry, idx) for idx in sorted_hit_idx)
    face_areas = [face_area(filtered_geometry, idx) for idx in sorted_hit_idx]
    _face_vertices = [face_vertices(filtered_geometry, idx) for idx in sorted_hit_idx]
    _face_normals = [face_normal(filtered_geometry, idx) for idx in sorted_hit_idx]

    # return Aproj, Aref, sorted_hit_idx, face_areas, _face_vertices, _face_normals, _face_vertices_preRT, _face_normals_preRT, _face_vertices_preCulling, _face_normals_preCulling, rti_vec, valid_rti, indices, geometry
    return Aproj, Aref, sorted_hit_idx, face_areas, _face_vertices, _face_normals, rti_vec, valid_rti, geometry
end


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
    # println("length filtered geo", length(filtered_geometry.faces))

    #Number of triangles
    Ntri = n_faces(filtered_geometry)
    # _face_vertices_preRT = [face_vertices(filtered_geometry, idx) for idx in 1:Ntri]
    # _face_normals_preRT = [face_normal(filtered_geometry, idx) for idx in 1:Ntri]

    new_viewpoint = shrink_viewpoint(filtered_geometry, viewpoint)
    println("culling ratio =", Ntri / Ntri_preculling)
    # println("max in MeshVertices=", new_viewpoint.rmax)
    # rt_vec, w, indices = _raytrace(filtered_geometry, new_viewpoint, sampler)
    rt_vec, w = _raytrace(filtered_geometry, new_viewpoint, sampler)

    # return rt_vec, w, _face_vertices_preRT, _face_normals_preRT, _face_vertices_preCulling, _face_normals_preCulling, indices, filtered_geometry
    return rt_vec, w, filtered_geometry
end

function _raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    dir = viewpoint.direction
    rmax = viewpoint.rmax
    O, Norig = generate_ray_origins(sampler, dir, rmax, viewpoint.distance)        #coordinates of ray origins, number of origins

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

    #------------------------------------------------
    # indices = zeros(Norig * Ntri)
    # counter = 0
    # for jj ∈ 1:Norig       #iterate over the set of ray origins
    #     orig = O[jj]
    #     ray = Ray(SV3(orig), dir)
    #     rti = ray_mesh_intersection(geometry, ray)
    #     if mode(rti) ∈ (BackFaceIntersection, FrontFaceIntersection)   #the triangle is intercepted by the ray
    #         ii = rti.face_index
    #         counter += 1
    #         indices[counter] = ii
    #     end
    # end

    # return rti_vec, (pi * rmax^2) / (Norig), indices

    return rti_vec, (pi * rmax^2) / (Norig)
end

export analyze_areas