using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "boxMesh.obj")
materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "facetMaterials.json")

#HETEROGENEOUS CASE
# geo = load_geometry(mesh_path, materials_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

# VERTICES = [SatelliteGeometryCalculations.face_vertices(geo, idx) for idx in 1:length(geo.faces)]
# writedlm("GRACE_VERTICES.txt", VERTICES)
#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

α = deg2rad(90)
ϕ = deg2rad(90)
# v = Viewpoint(geo, α, ϕ)
v = Viewpoint(geo, Vrel_v)

#DRIA - Sphere
# CD_sph, cd_j, sumM = SatelliteGeometryCalculations.DRIA_sphere(outSurfaceProps, outGasStreamProps, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Aref, intercept_info, normals, filteredGeo, culling, sampler2, dir2, rmax2, distance2 = analyze_areas(geo, v)



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
# #---------------------------------------------------------------------------------------------
println(coeffs)

SatelliteGeometryCalculations.tock()


# println("CD_inhouse = ", coeffs[1])
# println("CD_DRIA = ", CD_sph)

#CHECK TEST:
# using SatelliteGeometryCalculations, DelimitedFiles

# using FilePathsBase
# using FilePathsBase: /

# pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
# mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")
# geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

# α = deg2rad(90)
# ϕ = deg2rad(90)
# v = Viewpoint(geo, α, ϕ)

# abstract type RaySampler end
# struct FibonacciSampler <: RaySampler
#     ray_density::Float64 # rays / m²  

# end
# samplerF = FibonacciSampler(1e5)
# sampler = samplerF

# O, Norig, Aray = SatelliteGeometryCalculations.generate_ray_origins(sampler, v.direction, v.rmax, v.distance)

# rti_vec, w, filtered_geometry, Aray = SatelliteGeometryCalculations.raytrace(geo, v, sampler)

# rti_vec, w, filtered_geometry, Aray = SatelliteGeometryCalculations.areas_nonconvex(geo, v)