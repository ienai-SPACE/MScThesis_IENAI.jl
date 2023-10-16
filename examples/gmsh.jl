import Gmsh: gmsh
using FilePathsBase
using FilePathsBase: /

# --- Before using any functions in the Julia API, Gmsh must be initialized:------
gmsh.initialize()

# --- Merge a file. Equivalent to the `File->Merge` menu in the Gmsh app ---------
#=Input format: STEP, VRML, IGES, and STL =#
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphere_solid.stp")
gmsh.merge(string(path))

# --- Mesh size ------------------------------------------------------------------
lc = 10 #[mm] mesh size
gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
# gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

#= When the element size is fully specified by a mesh size field (as it is in
this example), it is thus often desirable to set. This will prevent over-refinement 
due to small mesh sizes on the boundary.=#

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# --- Algorithm ------------------------------------------------------------------
#= 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay,
6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of
Parallelograms, 11: Quasi-structured Quad) =#

gmsh.option.setNumber("Mesh.Algorithm", 5)

# --- Generate mesh of dimension -------------------------------------------------
gmsh.model.mesh.generate(2) #dim = 2

# --- Write mesh ------------------------------------------------------------------
#= Output mesh format: Mesh output format (1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos
stat, 27: stl, 28: p3d, 30: mesh, 31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39:
inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 49: neu, 50: matlab)=#
format = 27
gmsh.option.setNumber("Mesh.Format", format)

if format == 27
    gmshOut = gmsh.write("geometry_mesh.stl")
elseif 10 || 1
    gmshOut = gmsh.write("geometry_mesh.msh")
elseif 39
    gmshOut = gmsh.write("geometry_mesh.inp")
else
    println("Check format and write the mesh with the correct file extension")
end

# --- This should be called when you are done using the Gmsh Julia API: ----------
gmsh.finalize()