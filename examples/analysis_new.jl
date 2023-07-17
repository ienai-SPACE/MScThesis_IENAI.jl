using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
# mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "sphereMesh4.obj")
# load_geometry(mesh_path, SurfaceProps(), true)

mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "T_SatMesh.obj")
geo = load_geometry(mesh_path, SurfaceProps(), false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

# VERTICES = [SatelliteGeometryCalculations.face_vertices(geo, idx) for idx in 1:length(geo.faces)]
# writedlm("GRACE_VERTICES.txt", VERTICES)
#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

# α = deg2rad(-110)
# ϕ = deg2rad(-50)
# v = Viewpoint(geo, α, ϕ)
v = Viewpoint(geo, Vrel_v)

#DRIA - Sphere
# CD_sph, cd_j, sumM = SatelliteGeometryCalculations.DRIA_sphere(outSurfaceProps, outGasStreamProps, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Aref, intercept_info, normals = analyze_areas(geo, v)



# writedlm("hit_vertices.txt", hit_vertices)
# writedlm("hit_normals.txt", hit_normals)

# writedlm("_face_vertices_preRT.txt", _face_vertices_preRT)
# writedlm("_face_normals_preRT.txt", _face_normals_preRT)
# writedlm("_face_vertices_preCulling.txt", _face_vertices_preCulling)
# writedlm("_face_normals_preCulling.txt", _face_normals_preCulling)
# writedlm("_face_vertices_preRT2.txt", _face_vertices_preRT2)
# writedlm("_face_normals_preRT2.txt", _face_normals_preRT2)

# writedlm("rti_vec_t.txt", rti_vec.t)


println("Aproj = ", Aproj)
println("Aref = ", Aref)

#---- Aerodynamic coefficients ---------------------------------------------------------------

coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
#---------------------------------------------------------------------------------------------
println(coeffs)

SatelliteGeometryCalculations.tock()


# println("CD_inhouse = ", coeffs[1])
# println("CD_DRIA = ", CD_sph)

