"""
function : GeomInputs(Vrel_v, VdirFlag, convexFlag)

    GOAL:
        - Define whether the object is convex or non-convex
        - Define if the direction of "particle impingement" is set manually or by the velocity vector
        - Define the vertices coordinates of the meshed element
    INPUT:
        - [m/s] Relative velocity vector, flag for direction criterion, flag = 1 for convex shape and flag = 0 for non-convex
    OUTPUT:
        - MeshVerticesCoords : vertices coordinates of the meshed element
        - rmax               : radius of the source for ray-tracing algorithm
        - distance           : distance between source origin and local coordinate system for ray-tracing algorithm
"""

function GeomInputs(Vrel_v, VdirFlag, convexFlag)

    MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]



    if convexFlag != 0
        rmax = 1                            #radius of the circular plane from where rays originate
        distance = 10                       #distance at which the circular plane is located (it should be out from the satellite body)
    end

    if VdirFlag == 1
        #direction vector
        dir = -(Vrel_v / norm(Vrel_v))'
    elseif VdirFlag == 0
        dir = @SVector [1, 1, 0]
    end


    return MeshVerticesCoords, dir, rmax, distance

end