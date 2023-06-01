

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

    #Number of triangles
    Ntri = size(triangles, 1)

    if convexFlag == 1
        counter = 0
        for jj ∈ 1:Ntri #for loop to iterate over all triangles/quads
            vertices = triangles[jj, 1:9]
            OutAreaConvex = areasConvex(vertices, dir)

            if OutAreaConvex[1] != 0.0
                counter += 1
                if counter == 1
                    OutFacets = @SVector [jj, OutAreaConvex[1], OutAreaConvex[2], OutAreaConvex[3], OutAreaConvex[4], OutAreaConvex[5]]
                else
                    OutFacets = hcat(OutFacets, [jj, OutAreaConvex[1], OutAreaConvex[2], OutAreaConvex[3], OutAreaConvex[4], OutAreaConvex[5]])
                end
            end


        end

    elseif convexFlag == 0
        OutFacets = areasConcave(dir, rmax, distance, triangles, Ntri)
    end



    #store normal vector components into vector
    normal_v = Vector{SVector{3,Float64}}(undef, size(OutFacets, 2))
    for ii ∈ 1:lastindex(OutFacets[2, :])
        normal_v[ii] = SVector(OutFacets[4, ii], OutFacets[5, ii], OutFacets[6, ii])
    end
    # print(normal_v)

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

