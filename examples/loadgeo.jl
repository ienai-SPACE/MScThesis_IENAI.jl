using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "boxMesh.obj")
materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "facetMaterials.json")

#HETEROGENEOUS CASE
geo = load_geometry(mesh_path, materials_path, true, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
# geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters