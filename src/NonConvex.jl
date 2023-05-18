
# using StaticArrays, LinearAlgebra, GeometryBasics

# SV3{T} = SVector{3,T}
# include("Origins.jl")
# include("MTalgorithm.jl")

# include("GeometryInputs.jl")

# # triangles = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
# dir = @SVector [1.0, 1.0, 0.0]
# Vrel_v = [5395.145235865259 5395.145235865259 0.0];

# triangles, dir, rmax, distance = GeomInputs(Vrel_v, 0, 0)

# rmax = 3                            #radius of the circular plane from where rays originate
# distance = 10                       #distance at which the circular plane is located (it should be out from the satellite body)
# Ntri = size(triangles, 1)

# # OutTriangles = areasConcave(dir, rmax, distance, MeshVerticesCoords, Ntri)

#--------------------------------------------

""" 
    areasConcave(dir, rmax, distance, triangles, Ntri) 

Obtain the areas of triangular elements that are impinged by the oncoming particles. 
A ray-tracing algorithm is used. The sampler for the generation of the rays can be selected.

#INPUT:
- `triangles::Matrix`   : coordinates of each vertex in the format [x1y1z1x2y2z2x3y3z3;...]
- `dir`                 : vector with the direction to be analyzed
- `rmax`                : radius of the circular plane from where the rays are originated  
- `distance`            : distance of the perpendicular plane from where the rays are originated
#OUTPUT:
- `OutTriangles::Matrix(6,_)` : it stores index (1), area (2), angle (3), triangle's normal (4,5,6)
"""

function areasConcave(dir, rmax, distance, triangles, Ntri)

    #select the sampling method and the density of the sampler [rays/m^2]
    samplerG = GridFilter(50000)
    samplerF = FibonacciSampler(50000)
    samplerMC = MonteCarloSampler(50000)

    O, Norig = generate_ray_origins(samplerG, dir, rmax, distance)        #coordinates of ray origins, number of origins

    #------pre-allocation-------------------
    # index = 0                                          #counter indicating the number of triangles intercepted by the same ray
    intercept_dummy = zeros(Ntri, 7)                   #triangle index, distance from origin to intercept, area of the triangle, angle between velocity vector and triangle's normal
    triIntercept = zeros(6, Norig)                     #triangle index, area of the triangle, angle between velocity vector and triangle's normal
    #---------------------------------------

    for jj ∈ 1:Norig       #iterate over the set of ray origins
        global index
        index = 0
        for ii ∈ 1:Ntri    #iterate over all triangles

            #definition of the triangle vertices
            V1 = [triangles[ii, 1], triangles[ii, 2], triangles[ii, 3]]
            V2 = [triangles[ii, 4], triangles[ii, 5], triangles[ii, 6]]
            V3 = [triangles[ii, 7], triangles[ii, 8], triangles[ii, 9]]


            orig = O[jj, :]
            orig_new = reduce(vcat, orig)
            RTI, modo = MTalgorithm(TriangleFace(SV3(V1), SV3(V2), SV3(V3)), Ray(SV3(orig_new), dir))

            # print(RTI.t)
            # print(RTI.γ_dir)
            t = RTI.t
            gamma = RTI.γ_dir
            area = RTI.area
            u_n = TriangleFace(SV3(V1), SV3(V2), SV3(V3)).normal |> normalize
            # print(u_n)


            if modo == NoIntersection

            elseif modo == BackFaceIntersection || modo == FrontFaceIntersection   #the triangle is intercepted by the ray

                # print(index)
                global index

                index += 1

                intercept_dummy[index, :] = [ii, t, area, gamma, u_n[1], u_n[2], u_n[3]] #last 3 components are normal vector components of the facet

                if index > 1 #if more than one triangle intercepted by the same ray
                    if abs(intercept_dummy[index, 2]) < abs(intercept_dummy[index-1, 2])   #select the most forward triangle to be intercepted

                        triIntercept[:, jj] = [intercept_dummy[index, 1] intercept_dummy[index, 3] intercept_dummy[index, 4] intercept_dummy[index, 5] intercept_dummy[index, 6] intercept_dummy[index, 7]]   #triangle index, area of the triangle, angle between velocity vector and triangle's normal

                    end
                else
                    triIntercept[:, jj] = [intercept_dummy[index, 1] intercept_dummy[index, 3] intercept_dummy[index, 4] intercept_dummy[index, 5] intercept_dummy[index, 6] intercept_dummy[index, 7]]        #triangle index, area of the triangle, angle between velocity vector and triangle's normal
                end
            end
        end


        fill!(intercept_dummy, 0.0)

    end




    #------pre-allocation-------------------
    OutTriangles = zeros(6, Norig)            #matrix storing the triangle index, area, and the angle between the normal and the velocity direction
    triangleCount = 0
    indexList = zeros(Norig, 1)               #list to store the indices of the triangles that are intercepted
    #---------------------------------------

    #check the triangles that have been intercepted and store their properties
    for jj ∈ 1:Norig
        if Int(triIntercept[1, jj]) != 0
            if !(triIntercept[1, jj] ∈ indexList)

                triangleCount += 1
                OutTriangles[:, triangleCount] = triIntercept[:, jj]
                indexList[triangleCount] = triIntercept[1, jj]
            end
        end
    end

    # print(OutTriangles)

    #eliminate non-intercepted triangle entries and size down the output matrix
    OutTriangles = filter(!iszero, OutTriangles)

    # print(OutTriangles)

    OutTriangles = reshape(OutTriangles, (6, Int(lastindex(OutTriangles) / 6)))

    # lastindex(OutTriangles), length(OutTriangles), Int(length(OutTriangles) / 3)

    return OutTriangles

end

#TEST
#---------------------------------

# triangles = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
# dir = @SVector [1.0, 1.0, 0.0]
# Vrel_v = [5395.145235865259 5395.145235865259 0.0];

# # triangles, dir, rmax, distance = GeomInputs(Vrel_v, 0, 0)

# rmax = 3                            #radius of the circular plane from where rays originate
# distance = 10                       #distance at which the circular plane is located (it should be out from the satellite body)
# Ntri = size(triangles, 1)

# OutTriangles = areasConcave(dir, rmax, distance, triangles, Ntri)
