include("culling.jl")

"""
    normal_of_triangle_by_coords(triangle)

Return the normal of a triangle

#INPUT
- `triangle::Vector`
#OUTPUT
- `normal::SVector{3, T}`
"""

function normal_of_triangle_by_coords(triangle)
    V1 = SV3(triangle[1], triangle[2], triangle[3])
    V2 = SV3(triangle[4], triangle[5], triangle[6])
    V3 = SV3(triangle[7], triangle[8], triangle[9])
    TriangleFace(V1, V2, V3).normal
end

"""
    ray_mesh_intersection(triangles::Transpose{Float64,Matrix{Float64}}, ray::Ray)

Identify intercepted triangles by the MT algorithm and filter out back-face intersections. Add index to each intercepted triangle.

#INPUTS
- `triangles::Transpose{Float64,Matrix{Float64}}`
- `ray::Ray`
#OUTPUTS
- `rti :: MTalgorithm(face, ray)`
"""

function ray_mesh_intersection(triangles::Transpose{Float64,Matrix{Float64}}, ray::Ray)
    Ntri = size(triangles, 1)
    f = Filter(rti -> rti.mode ∈ (BackFaceIntersection, FrontFaceIntersection))
    m = Map(ii -> begin   #iterate over all triangles

        #definition of the triangle vertices
        V1 = SV3(triangles[ii, 1], triangles[ii, 2], triangles[ii, 3])
        V2 = SV3(triangles[ii, 4], triangles[ii, 5], triangles[ii, 6])
        V3 = SV3(triangles[ii, 7], triangles[ii, 8], triangles[ii, 9])

        # orig_new = reduce(vcat, orig)
        face = TriangleFace(V1, V2, V3)

        rti = MTalgorithm(face, ray)

        @set rti.face_index = ii
    end)
    foldl(earlier_intersection, 1:Ntri |> m |> f; init=no_intersection(eltype(triangles)))
end

"""
    ray_mesh_intersection(geo::HomogeneousGeometry, ray::Ray{T}) where {T}

Identify intercepted triangles by the MT algorithm and filter out back-face intersections. Add index to each intercepted triangle.

#INPUTS
- `geo::HomogeneousGeometry`
- `ray::Ray`
#OUTPUTS
- `rti :: MTalgorithm(face, ray)`
"""

function ray_mesh_intersection(geo::HomogeneousGeometry, ray::Ray{T}) where {T}
    Ntri = n_faces(geo)

    f = Filter(rti -> rti.mode ∈ (BackFaceIntersection, FrontFaceIntersection))
    m = Map(ii -> begin   #iterate over all triangles
        face = geo.faces[ii]
        rti = MTalgorithm(face, ray)
        @set rti.face_index = ii
    end)
    foldl(earlier_intersection, 1:Ntri |> m |> f; init=no_intersection(T))
end


""" 
    areasConcave(dir, rmax, distance, triangles, Ntri) 

Obtain the areas of triangular elements that are impinged by the oncoming particles. 
A ray-tracing algorithm is used. The sampler for the generation of the rays can be selected.

#INPUT:
- `MeshVertices::Matrix`   : coordinates of each vertex in the format [x1y1z1x2y2z2x3y3z3;...]
- `dir`                 : vector with the direction to be analyzed
- `rmax`                : radius of the circular plane from where the rays are originated  
- `distance`            : distance of the perpendicular plane from where the rays are originated
#OUTPUT:
- `OutTriangles::Matrix(6,_)` : it stores index (1), area (2), angle (3), triangle's normal (4,5,6)
"""
function areasConcave(dir, rmax, distance, MeshVertices, Ntri)

    # dir --> opposite direction to velocity vector
    Vdir = -dir

    Ntri_preculling = Ntri


    println("Ntri_preculling=", Ntri_preculling)
    println("max in MeshVertices (preculling)=", maximum(MeshVertices[:, :]))
    #back-face culling
    triangles = culling(MeshVertices, Vdir)
    #Number of triangles
    Ntri = size(triangles, 1)
    println("culling ratio =", Ntri / Ntri_preculling)
    println("max in MeshVertices=", maximum(MeshVertices[:, :]))


    #select the sampling method and the density of the sampler [rays/m^2]
    samplerG = GridFilter(50000)
    samplerF = FibonacciSampler(1e6)
    samplerMC = MonteCarloSampler(1e6)

    O, Norig = generate_ray_origins(samplerF, dir, rmax, distance)        #coordinates of ray origins, number of origins


    #------pre-allocation-------------------
    # index = 0                                          #counter indicating the number of triangles intercepted by the same ray
    intercept_dummy = zeros(Ntri, 7)                   #triangle index, distance from origin to intercept, area of the triangle, angle between velocity vector and triangle's normal
    triIntercept = zeros(6, Norig)                     #triangle index, area of the triangle, angle between velocity vector and triangle's normal
    #---------------------------------------

    for jj ∈ 1:Norig       #iterate over the set of ray origins
        orig = O[jj]
        ray = Ray(SV3(orig), dir)
        rti = ray_mesh_intersection(triangles, ray)
        if mode(rti) ∈ (BackFaceIntersection, FrontFaceIntersection)   #the triangle is intercepted by the ray
            ii = rti.face_index
            u_n = normal_of_triangle_by_coords(@view triangles[ii, :])
            triIntercept[:, jj] .= [
                ii, rti.area, rti.γ_dir, u_n[1], u_n[2], u_n[3]
            ]
        end
    end

    #------pre-allocation-------------------
    # OutTriangles = zeros(6, Norig)            #matrix storing the triangle index, area, and the angle between the normal and the velocity direction
    triangleCount = 0
    indexList = zeros(Norig, 1)               #list to store the indices of the triangles that are intercepted
    OutFacets = []
    #---------------------------------------

    #store the triangles that have been intercepted --> there could be triangles intercepted by more than 1 ray
    for jj ∈ 1:Norig
        if Int(triIntercept[1, jj]) != 0
            if !(triIntercept[1, jj] ∈ indexList)

                # triangleCount += 1
                # OutTriangles[:, triangleCount] = triIntercept[:, jj]
                # indexList[triangleCount] = triIntercept[1, jj]
                triangleCount += 1
                # print(triangleCount)
                if triangleCount == 1
                    OutFacets = @SVector [triIntercept[1, jj], triIntercept[2, jj], triIntercept[3, jj], triIntercept[4, jj], triIntercept[5, jj], triIntercept[6, jj]]
                    indexList[triangleCount] = triIntercept[1, jj]

                else
                    OutFacets = hcat(OutFacets, [triIntercept[1, jj], triIntercept[2, jj], triIntercept[3, jj], triIntercept[4, jj], triIntercept[5, jj], triIntercept[6, jj]])
                    indexList[triangleCount] = triIntercept[1, jj]
                end
            end
        end
    end

    return OutFacets

end

export areasConcave