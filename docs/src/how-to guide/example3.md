# Example 3: full execution of the code

This example shows how to run an analysis for the input geometry [TSAT_coarse_mesh.obj](@ref) with the material proporties shown in [TSAT_coarse_mesh_materials.json](@ref). Therefore, an **heterogeneous** case is run, having selected the **non-convex** method, with input units being [mm]. In this analysis the **viewpoint** is coincident with that of the velocity vector.

```
using SatelliteGeometryCalculations, DelimitedFiles
using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")
materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh_materials.json")

#HETEROGENEOUS CASE
geo = load_geometry(mesh_path, materials_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
# geo = load_geometry(mesh_path, true, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------

outSurfaceProps = SurfaceProps()                                             #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)      #outGasStreamProps.[C, PO, mmean, Ta]

# α = deg2rad(0)  #rotate around z-axis
# ϕ = deg2rad(90) #rotate around x-axis
# v = Viewpoint(geo, α, ϕ)
v = Viewpoint(geo, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Aref, intercept_info, normals, culling = analyze_areas(geo, v)

#---- Aerodynamic coefficients ---------------------------------------------------------------
coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
```