
#[A. Walker, P. M. Mehta, and J. Koller, “Drag coefficient model using the cercignani-lampis-lord gas- surface interaction model,” Journal of Spacecraft and Rockets, 2014]

using StaticArrays
using SatelliteToolbox

include("PermanentProperties.jl")
include("EnvironmentalConstants.jl")

function OxyPartPress(nrlmsise00_output)

    #den_Total[kg/m^3], the rest of densities[1/m^3], Temperature [K]
    ρ_v = @SVector [nrlmsise00_output.den_He, nrlmsise00_output.den_O, nrlmsise00_output.den_N2, nrlmsise00_output.den_O2, nrlmsise00_output.den_H, nrlmsise00_output.den_N]
    m_v = @SVector [He.m, O.m, N2.m, O2.m, H.m, N.m]


    nD = 0
    m_Total = 0
    for jj ∈ 1:6
        nd = ((ρ_v[jj] * m_v[jj] / NA / 1000) / m_v[jj]) * NA * 1000    #number density
        nD += nd                                                        #total number density
    end
    m_Total = sum(m_v)

    χ_O = ((nrlmsise00_output.den_O * O.m / NA / 1000) / O.m) / (nrlmsise00_output.den_Total / m_Total)
    Ta = nrlmsise00_output.T_alt                                        #ambient temperature

    P0 = nD * χ_O * Ta * kb                                             #[Pa]

    return P0

end

#TEST
#--------------------------------------------------------------------

#=
JD = date_to_jd(2018, 6, 19, 18, 35, 0);          #Julian Day [UTC].
alt = 300e3                                       #Altitude [m].
g_lat = deg2rad(-22)      #Geodetic latitude [rad].
g_long = deg2rad(-45)     #Geodetic longitude [rad].
f107A = 73.5       #81 day average of F10.7 flux (centered on day of year doy).
f107 = 79      #Daily F10.7 flux for previous day.
ap = 5.13        #Magnetic index.


nrlmsise00_output = nrlmsise00(JD, alt, g_lat, g_long, f107A, f107, ap, output_si=true, dversion=true)
OxyPartPress(nrlmsise00_output)
=#