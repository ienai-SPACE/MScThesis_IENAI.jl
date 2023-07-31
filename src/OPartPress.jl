
"""
    oxygen_partial_pressure(nrlmsise00_output::SatelliteToolbox.NRLMSISE00Output)

Calculate the partial pressure of oxygen in the atmosphere. Source: A. Walker, P. M. Mehta, and J. Koller, “Drag coefficient model using the cercignani-lampis-lord gas- surface interaction model,” Journal of Spacecraft and Rockets, 2014

# Input
-`nrlmsise00_output::SatelliteToolbox.NRLMSISE00Output`     : atmospheric model
# Output
- `P0`                                                      : oxygen partial pressure
"""
function oxygen_partial_pressure(nrlmsise00_output)

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

export oxygen_partial_pressure