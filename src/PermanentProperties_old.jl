# Permanent properties parameters are defined in this function: properties of the surface/material, gases

"""
    SurfaceAtomProps{T}

mass::T      :mass of an atom [g]
amass::T     :atomic mass [u]
"""
struct SurfaceAtomProps{T}

    mass::T
    amass::T

end

"""
    GasProps{T}

mass::T         :molecular mass [g/mol]
"""
struct GasProps{T}
    mass::T
end

"""
permProps(Gindex, Mindex, Sindex)

#INPUTS:
- `Gindex`        : gas element (mass) identifier
- `Mindex`        : material (solar radiation) identifier
- `Sindex`        : material surface (masses) identifier
#OUTPUTS
- `gasProperties`
"""

function permProps(Gindex, Mindex, Sindex)

    eta = 0   #generator efficiency 

    Tw = 300 #[K] temperature at the wall

    #= Solar radiation pressure coefficients (for calculation of force of the satellite due to
    solar radiation according to [Luthcke et al., 1997]) =#
    cr = 0.15 #specular reflectivity component
    cd = 0.25 #diffuse reflectivity component

    gasProperties = [GasProps(15.9994),
    ]
    push!(gasProperties, GasProps(15.9994)) #"O"
    push!(gasProperties, GasProps(1.0079))  #"H"
    push!(gasProperties, GasProps(14.0067)) #"N"
    push!(gasProperties, GasProps(31.9988)) #"O2"
    push!(gasProperties, GasProps(28.0134)) #"N2"
    push!(gasProperties, GasProps(2.0))       #"He"



    #Output parameters:

    return gasProperties[Gindex].mass


end