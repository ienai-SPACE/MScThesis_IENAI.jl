# using StaticArrays
# using FileIO
#-------------------------


"""
     finputMesh()

Transform the input `mesh` into the approapriate format: x1y1z1x2y2z2x3y3z3 such that `inputMesh = Matrix(undef, lastindex(mesh), 9)`    

-`MeshVerticesCoords`
"""

function finputMesh(mesh; scale=1e-3) # assumes mm input by default


    inputMesh = Matrix(undef, lastindex(mesh), 9)
    # s_matrix = SMatrix{lastindex(mesh),9,Float64}(undef)

    for ii ∈ 1:lastindex(mesh)
        #format: x1y1z1x2y2z2x3y3z3

        inputMesh[ii, 1] = mesh[ii].points[1][1]
        inputMesh[ii, 2] = mesh[ii].points[1][2]
        inputMesh[ii, 3] = mesh[ii].points[1][3]
        inputMesh[ii, 4] = mesh[ii].points[2][1]
        inputMesh[ii, 5] = mesh[ii].points[2][2]
        inputMesh[ii, 6] = mesh[ii].points[2][3]
        inputMesh[ii, 7] = mesh[ii].points[3][1]
        inputMesh[ii, 8] = mesh[ii].points[3][2]
        inputMesh[ii, 9] = mesh[ii].points[3][3]

    end
    # s_matrix = map(x -> x, inputMesh)
    Float64.(inputMesh .* scale)                 #pass from [mm] --> [m] & convert from Float32 to Float64
end

#TEST
#-----------------------------------------


# mesh = load("C:\\Users\\danie\\Documents\\UC3M\\IENAI internship\\CAD\\sphereMesh.obj")

# MeshVerticesCoords = finputMesh(mesh)



