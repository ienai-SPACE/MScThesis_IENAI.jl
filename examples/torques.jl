using SatelliteGeometryCalculations, DelimitedFiles, LinearAlgebra, StaticArrays

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphereMesh.obj")
materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh_materials.json")

#HETEROGENEOUS CASE
# geo = load_geometry(mesh_path, materials_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

α = deg2rad(180)  #rotate around z-axis
ϕ = deg2rad(0) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
# v = Viewpoint(geo, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Atot, intercept_info, normals, culling, barycenters = analyze_areas(geo, v);

println("Aproj = ", Aproj)
println("Aref = ", Atot)


#---- Aerodynamic coefficients ---------------------------------------------------------------

coeffs, Atot, Aproj, coeffs_v = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals);
# #---------------------------------------------------------------------------------------------
println("Cd = ", coeffs[1], coeffs[2])
println("Cl = ", coeffs[3], coeffs[4])
println("Cp = ", coeffs[5], coeffs[6])
println("Ctau = ", coeffs[7], coeffs[8])

#---- Center of pressure ---------------------------------------------------------------
#[Not valid: GSI depend on angle, not only on area contribution] 
# CoP = SatelliteGeometryCalculations.getCoP(Aproj, intercept_info, barycenters);
# println(CoP)

torque_ref = SVector(100.0, 0.0, 0.0)  # point about which moments are calculated
CT = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, barycenters, torque_ref)

chord_plane_n = SVector(0.0, 1.0, 0.0) #chord plane normal
CoP = SatelliteGeometryCalculations.getCoP(T, coeffs, chord_plane_n)

CT2 = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, barycenters, CoP[2])

(CT - CT2)
# #---------------------------------------------------------------------------------------------


SatelliteGeometryCalculations.tock()

