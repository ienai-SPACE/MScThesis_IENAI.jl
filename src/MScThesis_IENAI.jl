module MScThesis_IENAI

using LinearAlgebra

include("EnvironmentalConstants.jl")
include("PermanentProperties.jl")
include("MutableProps.jl")
include("GasStreamProps.jl")
include("CoefficientCalculations.jl")
include("Areas.jl")

outMutableProps = MutableProperties()        #outMutableProps.[η, Tw, s_cr, s_cd, m_srf]

#-----------NRLMSISE-00 INPUTS: Julian date, geodetic coordinates, and solar/magnetic indices-------------------------------
JD = date_to_jd(2018, 6, 19, 18, 35, 0)    #Julian Day [UTC].
alt = 300e3                                 #Altitude [m].
g_lat = deg2rad(-22)                        #Geodetic latitude [rad].
g_long = deg2rad(-45)                       #Geodetic longitude [rad].
f107A = 73.5                                #81 day average of F10.7 flux (centered on day of year doy).
f107 = 79                                   #Daily F10.7 flux for previous day.
ap = 5.13                                   #Magnetic index.
#----------------------------------------------------------------------------------------------------------------------------
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)        #outGasStreamProps.[C, PO, mmean, Ta]

#----Area calculation inputs-------------------------------------------------------------------------------------------------
triangles = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
dir = @SVector [1, 1, 0]

convexFlag = 1
rmax = 1                            #radius of the circular plane from where rays originate
distance = 10                       #distance at which the circular plane is located (it should be out from the satellite body)

#direction vector
#dir = [-0.6988494037821777, -0.137916655763437, -0.7018464981007773]
#dir = -(Vrel/norm(Vrel))';
#-----------------------------------------------------------------------------------------------------------------------------

Aproj, Aref, OutTriangles = areas(rmax, distance, dir, triangles, convexFlag)



#m_srf = Aluminum.amass                 #[g/mol] atomic mass of the surface atom
#Afacet = 1.5 .* ones(3, 1)             #[m^2]   area of each facet (set as input)
#δ = [0, 0, pi / 2]                     #[rad]   angle between the normal of the surface and the incoming particle   

#CD, CL, CP, CTAU = CoefficientCalculations(m_srf, OutTriangles[2, :], OutTriangles[3, :])
Vrel = 7000.0
CD, CL, CP, CTAU = CoefficientCalculations(outMutableProps, outGasStreamProps, OutTriangles, Vrel)


end