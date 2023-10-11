import Gmsh: gmsh
using FilePathsBase
using FilePathsBase: /

# Before using any functions in the Julia API, Gmsh must be initialized:
gmsh.initialize()




#Merge a file. Equivalent to the `File->Merge` menu in the Gmsh app.
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "Tsat.stp")
gmsh.merge(string(path))

#Mesh size
lc = 100 #[mm] mesh size
gmsh.model.mesh.gmsh.model.mesh.field.add("Threshold", 1)  #Create Field[1] 'threshold'(t10.jl used as reference)
gmsh.model.mesh.field.setNumber(1, "SizeMin", lc)
gmsh.model.mesh.field.setNumber(1, "SizeMax", 10 * lc)
# gmsh.model.mesh.field.setNumber(1, "DistMin", 0.15)
# gmsh.model.mesh.field.setNumber(1, "DistMax", 0.5)

#Generate mesh of dimension
gmsh.model.mesh.generate(2) #dim = 2

#Write mesh
gmshOut = gmsh.write("Tsat.msh")





# This should be called when you are done using the Gmsh Julia API:
gmsh.finalize()