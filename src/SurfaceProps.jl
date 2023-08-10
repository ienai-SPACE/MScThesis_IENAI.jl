""" 
    SurfaceProps{T}

Store properties of the satellite surface into a struct

- `η::T       generator efficiency`
- `Tw::T      [K] temperature at the wall `
- `s_cr::T    specular reflectivity component`
- `s_cd::T    diffuse reflectivity component`
- `m_srf::T   [g/mol] atomic mass of the surface atom`
"""
struct SurfaceProps{T}
    η::T       #generator efficiency
    Tw::T      #[K] temperature at the wall 
    s_cr::T    #specular reflectivity component
    s_cd::T    #diffuse reflectivity component
    m_srf::T   #[g/mol] atomic mass of the surface atom
end

"Default constructor - assumes aluminum at 300 K"
SurfaceProps() = SurfaceProps(Aluminum)


""" 
    SurfaceProps(sfatom::SurfaceAtomProperties)

Store properties of the satellite surface into a struct based on input properties

- `η::T       generator efficiency`
- `Tw::T      [K] temperature at the wall `
- `s_cr::T    specular reflectivity component`
- `s_cd::T    diffuse reflectivity component`
- `m_srf::T   [g/mol] atomic mass of the surface atom`
"""
function SurfaceProps(sfatom::SurfaceAtomProperties)
    m_srf = sfatom.amass                 #[g/mol] atomic mass of the surface atom
    η = 0.0   #generator efficiency  
    Tw = 300.0 #[K] temperature at the wall        
    #Solar radiation pressure coefficients (for calculation of force of the satellite due to solar radiation according to [Luthcke et al., 1997])
    s_cr = 0.15 #specular reflectivity component
    s_cd = 0.25 #diffuse reflectivity component    
    SurfaceProps(η, Tw, s_cr, s_cd, m_srf)
end

export SurfaceProps


"""
    struct Face{G,T}

- `geometry::G`
- `properties::SurfaceProps{T}`
"""
struct Face{G,T}
    geometry::G
    properties::SurfaceProps{T}
end