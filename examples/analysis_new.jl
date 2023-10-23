using SatelliteGeometryCalculations, DelimitedFiles, LinearAlgebra

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")
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
ϕ = deg2rad(34) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
# v = Viewpoint(geo, Vrel_v)

#---- DRIA - Sphere ------------------------------------------------------------------------
# CD_sph, cd_j, sumM = SatelliteGeometryCalculations.DRIA_sphere(outSurfaceProps, outGasStreamProps, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Atot, intercept_info, normals, culling, barycenters, solarCellsGeo, rti_vec, _rti, _rpi, faces_hit_idx_nonsorted = analyze_areas(geo, v);

println("Aproj = ", Aproj)
println("Aref = ", Atot)
println("Aproj_sc = ", solarCellsGeo.Aproj)


#---- Aerodynamic coefficients ---------------------------------------------------------------

coeffs, Atot, Aproj, coeffs_v = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals);
# #---------------------------------------------------------------------------------------------
println("Cd = ", coeffs[1], coeffs[2])
println("Cl = ", coeffs[3], coeffs[4])
println("Cp = ", coeffs[5], coeffs[6])
println("Ctau = ", coeffs[7], coeffs[8])

# #---------------------------------------------------------------------------------------------


#=
#h = 200
rho = 5.402465388584428e-7
v_i360 = 7655.309260090936
F_D = coeffs[1] * pi * 0.5^2 * 0.5 * v_i360^2 * rho


#h = 300
rho = 4.4053682747140424e-11
v_i360 = 7751.906068794503
F_D = coeffs[1] * pi * 0.5^2 * 0.5 * v_i360^2 * rho

#h = 400
rho = 2.3928535687949098e-12
v_i360 = 7668.372367380956
F_D = coeffs[1] * pi * 0.5^2 * 0.5 * v_i360^2 * rho
=#