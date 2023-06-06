

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

function areas_convex(rmax, distance, Vdir, triangles)
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

function areas_nonconvex(rmax, distance, Vdir, triangles)

    #Number of triangles
    Ntri = size(triangles, 1)
    OutFacets = areasConcave(Vdir, rmax, distance, triangles, Ntri)


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

#TESTING
#---------------------------------------------------------------------------------------------------------------------------------

#=
convexFlag = 0
rmax = 1                            #radius of the circular plane from where rays originate
distance = 10                       #distance at which the circular plane

#triangle vetices coordinates
triangles = @SMatrix[4.394897596952136 -1.3063587207678875 4.012655067923802 3.3442061435039823 0.7676371202562053 4.251153740481179 5.445214679778381 1.8739750984535304 5.439657623886108
    1.3410577524314113 0.6227469916184274 1.5604711027511295 1.4074299081230686 -0.5929514915580713 2.0821525791882842 2.850415168046063 0.5144988358467968 1.1637316942088742]
#Vrel = defined as global variable in 'EnvironmentalInputs'
#direction vector
dir = [-0.6988494037821777, -0.137916655763437, -0.7018464981007773]

triangles = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
dir = @SVector [1, 1, 0]

convexFlag = 1
rmax = 1                            #radius of the circular plane from where rays originate
distance = 10                       #distance at which the circular plane is located (it should be out from the satellite body)





Aproj, Aref, OutFacets = areas(rmax, distance, dir, triangles, convexFlag)
=#

function analyze_areas(geometry::AbstractGeometry, viewpoint::Viewpoint)
    if is_convex(geometry)
        return nothing #areas_convex(viewpoint, geometry)
    else
        return areas_nonconvex(geometry, viewpoint)
    end
end

function areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint)
    samplerG = GridFilter(50000)
    samplerF = FibonacciSampler(50000)
    samplerMC = MonteCarloSampler(50000)
    sampler = samplerF
    rti_vec, w = raytrace(geometry, viewpoint, sampler)
    valid_rti = rti_vec |> Filter(rti -> rti.mode != NoIntersection) |> tcollect
    Aproj = w * length(valid_rti)
    faces_hit_idx_nonunique = rti_vec .|> rti -> rti.face_index
    faces_hit_idx = unique(faces_hit_idx_nonunique) |> Filter(idx -> idx > 0) |> collect
    Aref = sum(face_area(geometry, idx) for idx in faces_hit_idx)
    Aproj, Aref
end


function raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    Ntri = n_faces(geometry)
    Ntri_preculling = Ntri
    rmax = viewpoint.rmax
    println("Ntri_preculling=", Ntri_preculling)
    println("max in MeshVertices (preculling)=", viewpoint.rmax)
    #back-face culling
    filtered_geometry = filter_backfaces(geometry, viewpoint)
    #Number of triangles
    Ntri = n_faces(filtered_geometry)
    new_viewpoint = shrink_viewpoint(filtered_geometry, viewpoint)
    println("culling ratio =", Ntri / Ntri_preculling)
    println("max in MeshVertices=", new_viewpoint.rmax)
    return _raytrace(filtered_geometry, new_viewpoint, sampler)
end

function _raytrace(geometry::AbstractGeometry, viewpoint::Viewpoint, sampler)
    dir = viewpoint.direction
    rmax = viewpoint.rmax
    O, Norig = generate_ray_origins(sampler, dir, rmax, viewpoint.distance)        #coordinates of ray origins, number of origins
    m = Map(jj -> begin
        orig = O[jj]
        ray = Ray(SV3(orig), dir)
        rti = ray_mesh_intersection(geometry, ray)
    end)

    rti_vec = 1:Norig |> m |> tcollect
    return rti_vec, (π * rmax^2) / Norig
end

export analyze_areas