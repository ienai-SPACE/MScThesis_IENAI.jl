using SatelliteGeometryCalculations, StaticArrays, LinearAlgebra
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs-------------------------------------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
#-----------------------------------------------------------------------------------------------------------------------------

outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

#----Area calculation inputs-------------------------------------------------------------------------------------------------
VdirFlag = 1          # 0: specify direction; 1: direction indicated by velocity vector
convexFlag = 0        # set if the satellite is convex (flag == 1) or non-convex (flag == 0)



MeshVerticesCoords, dir, rmax, distance = SatelliteGeometryCalculations.GeomInputs(Vrel_v, VdirFlag, convexFlag);     #mesh geometry & direction defined inside

# print(Vrel_v, //, dir, //, rmax, //, distance)


# MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
# dir = @SVector([1, 1, 0])
# distance = 10
# rmax = 2
# #-----------------------------------------------------------------------------------------------------------------------------
# print(dir, //, typeof(dir), //, typeof(rmax), //, typeof(distance))


Aproj, Atot, OutLMNTs, int_geos = SatelliteGeometryCalculations.areas(rmax, distance, dir, MeshVerticesCoords, convexFlag);      #calculation of areas and normals to the impinged surfaces


# α::Float64            : azimuth w.r.t. body frame. (α = 0 aligned with the x-axis)
# ϕ::Float64            : elevation w.r.t. body frame. (ϕ = 0 on the xy-plane) 
α = deg2rad(-142.416116)
ϕ = deg2rad(-71.911394)

Aproj, Atot, OutLMNTs, int_geos = SatelliteGeometryCalculations.areasSpherical(rmax, distance, α, ϕ, MeshVerticesCoords, convexFlag)


step = deg2rad(45);
LookUpTable, AprojLookUpTable = SatelliteGeometryCalculations.sweep(rmax, distance, MeshVerticesCoords, convexFlag, step)



# coeffs2, Atot2, Aproj2 = drag_for_orientation_convex(MeshVerticesCoords, outGasStreamProps, outSurfaceProps, α, ϕ, Vrel_norm)
# interactions_geometries = InteractionGeometry(OutLMNTs.area[1], OutLMNTs.angle[1])


# coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, int_geos, Vrel_v, OutLMNTs.normals)




