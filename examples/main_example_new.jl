using SatelliteGeometryCalculations

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = pkg_path / "test/inputs_models_data/sphereMesh4.obj"
geometry = load_geometry(mesh_path, SurfaceProps(), true)
