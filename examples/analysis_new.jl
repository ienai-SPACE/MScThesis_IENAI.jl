using SatelliteGeometryCalculations

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "sphereMesh4.obj")
load_geometry(mesh_path, SurfaceProps(), true)

mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "T_SatMesh.obj")
geo = load_geometry(mesh_path, SurfaceProps(), false)


α = deg2rad(0)
ϕ = deg2rad(90)
v = Viewpoint(geo, α, ϕ)
analyze_areas(geo, v)