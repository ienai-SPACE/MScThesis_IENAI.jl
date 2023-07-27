using SatelliteGeometryCalculations

using FilePathsBase
using FilePathsBase: /

# pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
# mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphereMesh4.obj")
# load_geometry(mesh_path, SurfaceProps(), true)

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "T_SatMesh.obj")
load_geometry(mesh_path, SurfaceProps(), false, "mm")