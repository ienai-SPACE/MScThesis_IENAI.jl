

function areasConcave(dir, rmax, distance, triangles, Ntri)


    O, Norig = RayOrigins(dir, rmax, distance)        #coordinates of ray origins, number of origins

    #------pre-allocation-------------------
    index = 0                                         #counter indicating the number of triangles intercepted by the same ray
    intercept_dummy = zeros(Ntri, 4)                   #triangle index, distance from origin to intercept, area of the triangle, angle between velocity vector and triangle's normal
    triIntercept = zeros(3, Norig)                     #triangle index, area of the triangle, angle between velocity vector and triangle's normal
    #---------------------------------------

    for jj ∈ 1:Norig       #iterate over the set of ray origins
        for ii ∈ 1:Ntri    #iterate over all triangles

            #definition of the triangle vertices
            V1 = [triangles[ii, 1], triangles[ii, 2], triangles[ii, 3]]
            V2 = [triangles[ii, 4], triangles[ii, 5], triangles[ii, 6]]
            V3 = [triangles[ii, 7], triangles[ii, 8], triangles[ii, 9]]

            area, flag, facing, u, v, t, γ_dir = MTalgorithm(O[jj, :], dir, V1, V2, V3)

            if flag == 0

            elseif flag == 1   #the triangle is intercepted by the ray

                global index
                index += 1

                intercept_dummy[index, :] = [ii, t, area, γ_dir]

                if index > 1 #if more than one triangle intercepted by the same ray
                    if abs(intercept_dummy[index, 2]) < abs(intercept_dummy[index-1, 2])   #select the most forward triangle to be intercepted

                        triIntercept[:, jj] = [intercept_dummy[index, 1] intercept_dummy[index, 3] intercept_dummy[index, 4]]   #triangle index, area of the triangle, angle between velocity vector and triangle's normal

                    end
                else
                    triIntercept[:, jj] = [intercept_dummy[index, 1] intercept_dummy[index, 3] intercept_dummy[index, 4]]           #triangle index, area of the triangle, angle between velocity vector and triangle's normal
                end
            end
        end

        global index
        index = 0
        fill!(intercept_dummy, 0)

    end


    #------pre-allocation-------------------
    OutTriangles = zeros(3, Norig)            #matrix storing the triangle index, area, and the angle between the normal and the velocity direction
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

    #eliminate non-intercepted triangle entries and size down the output matrix
    OutTriangles = filter(!iszero, OutTriangles)
    OutTriangles = reshape(OutTriangles, (3, Int(length(OutTriangles) / 3)))


    return OutTriangles

end