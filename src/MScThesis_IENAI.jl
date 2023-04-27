module MScThesis_IENAI

using LinearAlgebra, StaticArrays, SatelliteToolbox

SV3{T} = SVector{3,T}

include("geometry.jl")
include("EnvironmentalConstants.jl")
include("PermanentProperties.jl")
include("SurfaceProps.jl")
include("GasStreamProps.jl")

include("CoefficientCalculations.jl") # -----

include("Areas.jl")
include("GeometryInputs.jl")
include("Orbit&DateInputs.jl")


outSurfaceProps = SurfaceProps()                                                   #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs-------------------------------------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = OrbitandDate()
#-----------------------------------------------------------------------------------------------------------------------------

outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)        #outGasStreamProps.[C, PO, mmean, Ta]

#----Area calculation inputs-------------------------------------------------------------------------------------------------
VdirFlag = 0                                                                            #set whether direction is set by the velocity vector or by a predifined vector inside the function
convexFlag = 1                                                                          #set if the satellite is convex (flag == 1) or non-convex (flag == 0)
MeshVerticesCoords, dir, rmax, distance = GeomInputs(Vrel_v, VdirFlag, convexFlag)      #mesh geometry defined inside
# #-----------------------------------------------------------------------------------------------------------------------------

Aproj, Aref, OutLMNTs, int_geos = areas(rmax, distance, dir, MeshVerticesCoords, convexFlag)      #calculation of areas and normals to the impinged surfaces

# struct InteractionGeometry{T}
#     area::T
#     angle::T
# end

interactions_geometries = InteractionGeometry(OutLMNTs.area[1], OutLMNTs.angle[1])

Vrel_norm = 7000.0
#surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::Vector{<:InteractionGeometry}, Vrel_norm
coeffs, A = compute_coefficients(outSurfaceProps, outGasStreamProps, int_geos, Vrel_norm)

(; Cd, Cl, Cp, Ctau) = coeffs

end
