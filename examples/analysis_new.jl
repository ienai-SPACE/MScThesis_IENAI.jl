using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")
materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh_materials.json")

#HETEROGENEOUS CASE
geo = load_geometry(mesh_path, materials_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
# geo = load_geometry(mesh_path, true, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

# VERTICES = [SatelliteGeometryCalculations.face_vertices(geo, idx) for idx in 1:length(geo.faces)]
# writedlm("GRACE_VERTICES.txt", VERTICES)
#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

α = deg2rad(0)  #rotate around z-axis
ϕ = deg2rad(90) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
# v = Viewpoint(geo, Vrel_v)

#DRIA - Sphere
# CD_sph, cd_j, sumM = SatelliteGeometryCalculations.DRIA_sphere(outSurfaceProps, outGasStreamProps, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Aref, intercept_info, normals, culling = analyze_areas(geo, v)

println("Aproj = ", Aproj)
println("Aref = ", Aref)

#---- Aerodynamic coefficients ---------------------------------------------------------------

coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
# #---------------------------------------------------------------------------------------------
println(coeffs)

SatelliteGeometryCalculations.tock()


# println("CD_inhouse = ", coeffs[1])
# println("CD_DRIA = ", CD_sph)

