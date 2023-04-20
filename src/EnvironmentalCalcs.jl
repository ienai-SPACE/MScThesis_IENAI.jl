
using StaticArrays
using SatelliteToolbox
include("PermanentProperties.jl")
include("EnvironmentalConstants.jl")
include("OPartPress.jl")

function fEnvironmentlCalcs(JD, alt, g_lat, g_long, f107A, f107, ap)

    nrlmsise00_output = nrlmsise00(JD, alt, g_lat, g_long, f107A, f107, ap, output_si=true, dversion=true) #JD::Number, alt::Number, g_lat::Number, g_long::Number [, f107A::Number, f107::Number, ap::Union{Number,AbstractVector}]; output_si::Bool = true, dversion::Bool = true

    #Concentrations
    C = @MVector [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    C[1] = nrlmsise00_output.den_He / nrlmsise00_output.den_Total * (He.m / NA / 1000)
    C[2] = nrlmsise00_output.den_O / nrlmsise00_output.den_Total * (O.m / NA / 1000)
    C[3] = nrlmsise00_output.den_N2 / nrlmsise00_output.den_Total * (N2.m / NA / 1000)
    C[4] = nrlmsise00_output.den_O2 / nrlmsise00_output.den_Total * (O2.m / NA / 1000)
    C[5] = nrlmsise00_output.den_H / nrlmsise00_output.den_Total * (H.m / NA / 1000)
    C[6] = nrlmsise00_output.den_N / nrlmsise00_output.den_Total * (N.m / NA / 1000)

    #Oxygen partial pressure
    PO = OxyPartPress(nrlmsise00_output)

    #Mean molecular mass
    num = sum([nrlmsise00_output.den_He * He.m, nrlmsise00_output.den_O * O.m, nrlmsise00_output.den_N2 * N2.m, nrlmsise00_output.den_O2 * O2.m, nrlmsise00_output.den_H * H.m, nrlmsise00_output.den_N * N.m])
    den = sum([nrlmsise00_output.den_He, nrlmsise00_output.den_O, nrlmsise00_output.den_N2, nrlmsise00_output.den_O2, nrlmsise00_output.den_H, nrlmsise00_output.den_N])
    mmean = num / den      #[g/mol]

    return C, nrlmsise00_output.T_alt, PO, mmean

end

#TEST
#---------------------------------
#
JD = date_to_jd(2018, 6, 19, 18, 35, 0)          #Julian Day [UTC].
alt = 300e3                                       #Altitude [m].
g_lat = deg2rad(-22)      #Geodetic latitude [rad].
g_long = deg2rad(-45)     #Geodetic longitude [rad].
f107A = 73.5       #81 day average of F10.7 flux (centered on day of year doy).
f107 = 79      #Daily F10.7 flux for previous day.
ap = 5.13        #Magnetic index.

C, nrlmsise00_output.T_alt, PO, mmean = fEnvironmentlCalcs(JD, alt, g_lat, g_long, f107A, f107, ap)

#