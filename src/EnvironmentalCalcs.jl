
using StaticArrays
using SatelliteToolbox
# include("PermanentProperties.jl")
# include("EnvironmentalConstants.jl")
include("OPartPress.jl")

"""
    GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)

Computes the concentration of the 6 gases in the atmosphere, the temperature, the oxygen partial pressure and the mean molecular mass.
"""
function GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)
    nrlmsise00_output = nrlmsise00(JD, alt, g_lat, g_long, f107A, f107, ap, output_si=true, dversion=true) #JD::Number, alt::Number, g_lat::Number, g_long::Number [, f107A::Number, f107::Number, ap::Union{Number,AbstractVector}]; output_si::Bool = true, dversion::Bool = true
    GasStreamProperties(nrlmsise00_output)
end

function GasStreamProperties(nrlmsise00_output::SatelliteToolbox.NRLMSISE00_Output)
    out = nrlmsise00_output

    #Concentrations
    m_v = @SVector [He.m, O.m, N2.m, O2.m, H.m, N.m]
    den_v = @SVector [out.den_He, out.den_O, out.den_N2, out.den_O2, out.den_H, out.den_N]
    C = @. den_v * m_v / (NA * 1000 * nrlmsise00_output.den_Total)

    #Oxygen partial pressure
    PO = oxygen_partial_pressure(nrlmsise00_output)

    #Mean molecular mass
    num = den_v'm_v
    den = sum(den_v)
    mmean = num / den      #[g/mol]
    GasStreamProperties(C, PO, mmean, nrlmsise00_output.T_alt)
end

export GasStreamProperties


#TEST
#---------------------------------
#=
JD = date_to_jd(2018, 6, 19, 18, 35, 0)          #Julian Day [UTC].
alt = 300e3                                       #Altitude [m].
g_lat = deg2rad(-22)      #Geodetic latitude [rad].
g_long = deg2rad(-45)     #Geodetic longitude [rad].
f107A = 73.5       #81 day average of F10.7 flux (centered on day of year doy).
f107 = 79      #Daily F10.7 flux for previous day.
ap = 5.13        #Magnetic index.

C, Talt, PO, mmean = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)

=#
