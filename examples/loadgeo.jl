using SatelliteGeometryCalculations

using FilePathsBase
using FilePathsBase: /

# pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
# mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "sphereMesh4.obj")
# load_geometry(mesh_path, SurfaceProps(), true)

mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "T_SatMesh.obj")
load_geometry(mesh_path, SurfaceProps(), false, "mm")