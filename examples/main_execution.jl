using SatelliteGeometryCalculations, StaticArrays, LinearAlgebra
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs-------------------------------------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
#-----------------------------------------------------------------------------------------------------------------------------

outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

#----Area calculation inputs-------------------------------------------------------------------------------------------------
VdirFlag = 0           # 0: specify direction; 1: direction indicated by velocity vector
convexFlag = 0         # set if the satellite is convex (flag == 1) or non-convex (flag == 0)
MeshVerticesCoords, dir, rmax, distance = SatelliteGeometryCalculations.GeomInputs(Vrel_v, VdirFlag, convexFlag)     #mesh geometry & direction defined inside
# MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
# dir = @SVector([1, 1, 0])
# distance = 10
# rmax = 2
# #-----------------------------------------------------------------------------------------------------------------------------

Aproj, Atot, OutLMNTs, int_geos = SatelliteGeometryCalculations.areas(rmax, distance, dir, MeshVerticesCoords, convexFlag);      #calculation of areas and normals to the impinged surfaces


Vrel_norm = norm(Vrel_v)


# coeffs2, Atot2, Aproj2 = drag_for_orientation_convex(MeshVerticesCoords, outGasStreamProps, outSurfaceProps, α, ϕ, Vrel_norm)

interactions_geometries = InteractionGeometry(OutLMNTs.area[1], OutLMNTs.angle[1])

coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, int_geos, Vrel_v, OutLMNTs.normals)



