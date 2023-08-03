# Example 1: loading the geometry and material inputs

Three main options are available to tailor the way the geometry and material properties are introduced, namely the material composition of the geometry, and the units in which it is introduced.

* Firstly, it is necessary to differentiate between a homogeneous (single material) or heterogenous (multi-material) analysis

## Homogeneous case

```
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "boxMesh.obj")

geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
```
## Heterogeneous case

```
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "boxMesh.obj")
materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "facetMaterials.json")

geo = load_geometry(mesh_path, materials_path, true, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

```

* The second parameter that can be modified is the units in which the input geometry is introduced, namely [m] or [mm]. Selected as a `string` in the function [`SatelliteGeometryCalculations.load_geometry`](\ref) 

* The third parameter that must be selected is whether the introduced geometry is **convex** or **non-convex**. This is done by selecting in [`SatelliteGeometryCalculations.load_geometry`](\ref) the `boolean` input as `true` or `false`, respectively.