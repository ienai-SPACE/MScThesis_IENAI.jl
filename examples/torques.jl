using SatelliteGeometryCalculations, DelimitedFiles, LinearAlgebra, StaticArrays
SGC = SatelliteGeometryCalculations

# SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphere_centered.obj")
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

α = deg2rad(0)  #rotate around z-axis
ϕ = deg2rad(90) #rotate around x-axis
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

CoP, = SGC.getCoP(CT_A, coeffs, Aproj)
CT2_A = SGC.getAeroTorque(coeffs, CoP.CoPx, torque_ref, Aproj)

println("CT=", CT_A)
println("CT=", CT2_A)

CoP, = SGC.getCoP(CT_A, coeffs, Aproj)

CT3_A = SGC.getTorques(coeffs_v, intercept_info, CoP.CoPz, ray_facet_info)
println("CT=", norm(CT3_A))

# #------------------------------------------------------------------------------------------------


# SatelliteGeometryCalculations.tock()


#..........................................................
# coords = [ray_coords[ii].coords for ii in 1:lastindex(ray_coords)]
# writedlm("coords", coords)
# f_hit_idx = [faces_hit_idx_nonunique[ii] for ii in 1:lastindex(faces_hit_idx_nonunique)]
# writedlm("f_hit_idx", faces_hit_idx_nonunique)
# idx = [hit_idx[ii] for ii in 1:lastindex(hit_idx)]
# writedlm("hit_idx", idx)
# coeffs_v = [coeffs_v[ii] for ii in 1:lastindex(hit_idx)]
# writedlm("coeffs", coeffs_v)

# torque_ref = SVector(0.0, 0.0, 0.0)  # point about which moments are calculated
# 	CT_A = SGC.getTorques(coeffs_v, intercept_info, torque_ref, ray_facet_info)
# 	CoP, = SGC.getCoP(CT_A, coeffs, Aproj)

# 	CT_Ax = SGC.getTorques(coeffs_v, intercept_info, CoP.CoPx, ray_facet_info)
#     CT_Ay = SGC.getTorques(coeffs_v, intercept_info, CoP.CoPy, ray_facet_info)
#     CT_Az = SGC.getTorques(coeffs_v, intercept_info, CoP.CoPz, ray_facet_info)

# 	checkx = norm(CT_Ax) #center of pressure on the chord plane
# 	checky = norm(CT_Ay) #center of pressure on the chord plane
#     checkz = norm(CT_Az) #center of pressure on the chord plane