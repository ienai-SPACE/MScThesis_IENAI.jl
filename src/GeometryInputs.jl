using FileIO
# using GeometryBasics

include("inputMesh.jl")

"""
    GeomInputs(Vrel_v, VdirFlag, convexFlag)

Define whether the object is convex or non-convex
Define if the direction of "particle impingement" is set manually or by the velocity vector
Load the the meshed element

#INPUT:
-`Vrel_v`      : [m/s] relative velocity vector
- `VdirFlag`   : flag for direction criterion
- `convexFlag` : flag = 1 for convex shape and flag = 0 for non-convex
#OUTPUT:
- `MeshVerticesCoords`
- `dir`
- `rmax`, distance
        - MeshVerticesCoords : vertices coordinates of the meshed element
        - rmax               : radius of the source for ray-tracing algorithm
        - distance           : distance between source origin and local coordinate system for ray-tracing algorithm
"""

function GeomInputs(Vrel_v, VdirFlag, convexFlag)

    # MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]

    #load the mesh
    mesh = load("C:\\Users\\danie\\Documents\\UC3M\\IENAI internship\\CAD\\sphereMesh4.obj")

    MeshVerticesCoords = finputMesh(mesh)

    if convexFlag == 0
        rmax = maximum(MeshVerticesCoords[:, :])          #radius of the circular plane from where rays originate
        distance = rmax * 10                                                    #distance at which the circular plane is located (it should be out from the satellite body)
    else
        rmax = 0
        distance = 0
    end

    if VdirFlag == 1
        #direction vector
        dir = (Vrel_v / norm(Vrel_v))'
    elseif VdirFlag == 0
        dir = @SVector [1 / sqrt(2), 1 / sqrt(2), 0]
    end


    return MeshVerticesCoords, dir, rmax, distance

end