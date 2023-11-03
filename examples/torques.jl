using SatelliteGeometryCalculations, DelimitedFiles, LinearAlgebra, StaticArrays
SGC = SatelliteGeometryCalculations

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

α = deg2rad(90)  #rotate around z-axis
ϕ = deg2rad(0) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
# v = Viewpoint(geo, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Atot, intercept_info, normals, culling, solarCellsGeo, rti_vec, ray_facet_info = analyze_areas(geo, v);

println("Aproj = ", Aproj)
println("Aref = ", Atot)

#---- Aerodynamic coefficients ---------------------------------------------------------------

coeffs, Atot, Aproj, coeffs_v = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals);

println("Cd = ", coeffs[1], coeffs[2])
println("Cl = ", coeffs[3], coeffs[4])
println("Cp = ", coeffs[5], coeffs[6])
println("Ctau = ", coeffs[7], coeffs[8])

#---- Center of pressure ------------------------------------------------------------------------

# For the calculation of the aerodynamic torque, the `torque_ref` should be the CoM 
torque_ref = SVector(0.0, 0.0, 0.0)
CT_A = SGC.getTorques(coeffs_v, intercept_info, torque_ref, ray_facet_info)

CoP = SGC.getCoP(CT_A, coeffs, Aproj)
CT2_A = SGC.getAeroTorque(coeffs, CoP.CoPx, torque_ref, Aproj)

println("CT=", CT_A)
println("CT=", CT2_A)

CoP = SGC.getCoP(CT_A, coeffs, Aproj)

CT3_A = SGC.getTorques(coeffs_v, intercept_info, CoP.CoPy, ray_facet_info)
println("CT=", CT3_A)

# #------------------------------------------------------------------------------------------------


SatelliteGeometryCalculations.tock()


#-------------------------------
coords_v = map(ii-> ray_facet_info.ray_coords[ii].coords, 1:77792)
hit_idx = map(ii-> ray_facet_info.hit_idx[ii], 1:332)
f_hit_idx = map(ii-> ray_facet_info.f_hit_idx[ii], 1:77792)
writedlm("coeffs_v", coeffs_v)
writedlm("coords_v", coords_v)
writedlm("hit_idx", hit_idx)
writedlm("f_hit_idx", f_hit_idx)

