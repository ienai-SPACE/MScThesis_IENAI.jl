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

α = deg2rad(90)  #rotate around z-axis
ϕ = deg2rad(0) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
# v = Viewpoint(geo, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Atot, intercept_info, normals, culling, solarCellsGeo, rti_vec, _rti, _rpi, faces_hit_idx_nonsorted, hit_idx = analyze_areas(geo, v);

println("Aproj = ", Aproj)
println("Aref = ", Atot)


#---- Aerodynamic coefficients ---------------------------------------------------------------

coeffs, Atot, Aproj, coeffs_v = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals);

println("Cd = ", coeffs[1], coeffs[2])
println("Cl = ", coeffs[3], coeffs[4])
println("Cp = ", coeffs[5], coeffs[6])
println("Ctau = ", coeffs[7], coeffs[8])

#---- Center of pressure ------------------------------------------------------------------------
#[Not valid: GSI depend on angle, not only on area contribution] 
# CoP = SatelliteGeometryCalculations.getCoP(Aproj, intercept_info, barycenters);
# println(CoP)

max_idx = maximum([_rti[ii].idx for ii in 1:lastindex(_rti)])
min_idx = minimum([_rti[ii].idx for ii in 1:lastindex(_rti)])
max_idx2 = maximum([rti_vec[ii].face_index for ii in 1:lastindex(rti_vec)])

# For the calculation of the aerodynamic torque, the `torque_ref` should be the CoM 
torque_ref = SVector(0.0, 0.0, 0.0)  # point about which moments are calculated
CT = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, torque_ref, _rti, _rpi)

chord_plane_n = SVector(0.0, 1.0, 0.0) #chord plane normal
CoP = SatelliteGeometryCalculations.getCoP(CT, coeffs, Aproj, chord_plane_n)
CT2 = SatelliteGeometryCalculations.getAeroTorque(coeffs, CoP[2], torque_ref, Aproj)

println("CT=", CT)
println("CT=", CT2)


chord_plane_n = SVector(0.0, 1.0, 0.0) #chord plane normal
CoP = SatelliteGeometryCalculations.getCoP(CT, coeffs, Aproj, chord_plane_n)

CT3 = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, CoP[2], _rti, _rpi)
println("CT=", CT3)

# #------------------------------------------------------------------------------------------------


SatelliteGeometryCalculations.tock()


