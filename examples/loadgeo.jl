using SatelliteGeometryCalculations

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = pkg_path / "test/samples/sphereMesh4.obj"
load_geometry(mesh_path, SurfaceProps(), true)