using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
# mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "sphereMesh4.obj")
# load_geometry(mesh_path, SurfaceProps(), true)

mesh_path = FilePathsBase.join(pkg_path, "test", "samples", "T_Sat_fineMesh.obj")
geo = load_geometry(mesh_path, SurfaceProps(), false)

#---------- EVALUATION OF A SINGLE VIEWPOINT DIRECTION--------------------------------------
α = deg2rad(90)
ϕ = deg2rad(180)
v = Viewpoint(geo, α, ϕ)

##################################################################################################################
#---- Non-convex----------------------
Aproj, Aref, hit_indices, hit_areas, hit_vertices, hit_normals, rti_vec, valid_rti, geometry = analyze_areas(geo, v)
#---- Convex--------------------------
# Aproj, Aref = analyze_areas(geo, v)
##################################################################################################################

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

#-------------------------------------------------------------------------------------------

#---------- SWEEP TO GENERATE LOOK-UP TABLE-------------------------------------------------

#Perform sweep
# step = deg2rad(15);
# LookUpTable, AprojLookUpTable = SatelliteGeometryCalculations.sweep_v2(geo, step)

# writedlm("AprojLookUpTable_v2.txt", AprojLookUpTable)
#-------------------------------------------------------------------------------------------



SatelliteGeometryCalculations.tock()
