using StaticArrays
using SatelliteToolbox

#include("EnvironmentalConstants.jl")
include("EnvironmentalCalcs.jl")

"""
function : GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)

    GOAL: make environmental calculations based on the atmospheric model and store the results in a struct

    INPUT:
        - Julian date, altitude, g_lat, g_long, f107A, f107, ap, relative velocity vector
    
    OUTPUT:
        - struct : stGasStreamProps(C, PO, mmean, Ta)
            - Gases concentrations, [Pa] partial oxygen pressure, [g/mol] mean molecular mass, [K] ambient temperature
"""

struct stGasStreamProps{T}
    #= a::T           #[km]
     e::T
     r_p::T         #[km]
     r_a::T         #[km]
     h::T           #[km]   mean altitude of the orbit            
     Vatm::T        #[km/s] (h = 200-300 km) velocity of rotation of the atmosphere
     Vsc::T         #[km/s] velocity of the spacecraft
     Vrel::T        #[m/s] relative velocity between spacecraft and atmosphere
     =#
    C::MVector{6,T}           #gas concentratiosn: C = [O H N O2 N2 He]
    PO::T          #[Pa] for 300km (5×10−6;10−6 for 400km and 500km, respectively) computed with NRLMSISE-00
    mmean::T       #[g/mol] mean molecular mass
    Ta::T          #[K] ambient temperature
end


function GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)

    #Environmental calculations derived from the atmospheric model
    C, Ta, PO, mmean = fEnvironmentlCalcs(JD, alt, g_lat, g_long, f107A, f107, ap)  #Output: [-], [K], [Pa], [g/mol]

    return stGasStreamProps(C, PO, mmean, Ta)

end

#TEST
#------------------------------------------------------------------

#=
#-----------NRLMSISE-00 INPUTS: Julian date, geodetic coordinates, and solar/magnetic indices-------------------------------
JD = date_to_jd(2018, 6, 19, 18, 35, 0);    #Julian Day [UTC].
alt = 300e3                                 #Altitude [m].
g_lat = deg2rad(-22)                        #Geodetic latitude [rad].
g_long = deg2rad(-45)                       #Geodetic longitude [rad].
f107A = 73.5                                #81 day average of F10.7 flux (centered on day of year doy).
f107 = 79                                   #Daily F10.7 flux for previous day.
ap = 5.13                                   #Magnetic index.
#----------------------------------------------------------------------------------------------------------------------
=#

#GASOUT = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)
