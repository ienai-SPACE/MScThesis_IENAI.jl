# Meshing

## GMSH

Gmsh is a three-dimensional finite element mesh generator with a build-in CAD engine and post-processor. Its design goal is to provide a fast, light and user-friendly meshing tool with parametric input and flexible visualization capabilities.

**Input format:**

For the specific purpose of GMSH in this application, the **input CAD model** must be provided in **\*.stp** format.

**Output format:**

The rendered mesh will be in **\*.stl ASCII** format.

## License info

Source: *C. Geuzaine and J.-F. Remacle, Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering, Volume 79, Issue 11, pages 1309-1331, 2009.*

The precise conditions of the license for Gmsh are found in the General Public License that accompanies the source code (see Appendix F [License], page 409). Further information about this license is available from the GNU Project webpage https://www.gnu.org/copyleft/gpl-faq.html. Detailed copyright information can be found in Appendix E [Copyright and credits], page 405. If you want to integrate parts of Gmsh into a closed-source software, or want to sell a modified closed-source version of Gmsh, you will need to obtain a different license. Please contact us directly for more information.

## Example

Import gmsh API:
```
import Gmsh: gmsh
using FilePathsBase
using FilePathsBase: /
```
Before using any functions in the Julia API, Gmsh must be initialized:
```
gmsh.initialize()
```
Merge a file. Equivalent to the `File->Merge` menu in the Gmsh GUI.
* Input format: STEP (STL does not easily allow to re-mesh)
```
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "T_SatMesh.stp")
gmsh.merge(string(path))
```
**Mesh size:**
```
lc = 100 #[mm] mesh size
gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
# gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
```
When the element size is fully specified by a mesh size field (as it is in this example), it is thus often desirable to set. This will prevent over-refinement due to small mesh sizes on the boundary:
```
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
```
**Algorithm:**
* 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad)

```
gmsh.option.setNumber("Mesh.Algorithm", 5)
```

Generate mesh of dimension:
```
gmsh.model.mesh.generate(2) #dim = 2
```
Write mesh:
*Output mesh format: Mesh output format (1: msh, 2: unv, 10: auto, 16: vtk, 19: vrml, 21: mail, 26: pos
stat, 27: stl, 28: p3d, 30: mesh, 31: bdf, 32: cgns, 33: med, 34: diff, 38: ir3, 39: inp, 40: ply2, 41: celum, 42: su2, 47: tochnog, 49: neu, 50: matlab)
```
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
```
This should be called when you are done using the Gmsh Julia API:
```
gmsh.finalize()
```
