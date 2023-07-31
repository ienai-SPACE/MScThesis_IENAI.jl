"""
    Gas{T}

Gas atomic weights

# Fields
- `m::T`    # molecular mass [g/mol]
"""
struct Gas{T}
    m::T
end

const O = Gas(15.9994)
const H = Gas(1.0079)
const N = Gas(14.0067)
const O2 = Gas(31.9988)
const N2 = Gas(28.0134)
const He = Gas(2.0)
const atmospheric_gases = @SVector[O, H, N, O2, N2, He]

"""
    SurfaceAtomProperties{T}

Surface material atomic masses

- `mass::T`      #mass of an atom [g]
- `amass::T`     #atomic mass [u] 
"""
struct SurfaceAtomProperties{T}
    mass::T      #mass of an atom [g]
    amass::T     #atomic mass [g/mol]
end

SurfaceAtomProperties(; atomic_mass=amass) = SurfaceAtomProperties(atomic_mass / 6.02214076e23, atomic_mass)

const Aluminum = SurfaceAtomProperties(4.48e-23, 26.9815)
const Silicon = SurfaceAtomProperties(4.664e-23, 28.0855)
const Monocrystalline_Silicon = SurfaceAtomProperties(5.333651448335449e-23, 32.12)

"""
    SolarMatProps{T} 

Solar radiation material properties

- `ads::T`     adsorptance (adimensional)
- `emitt::T`   emittance (adimensional)
"""
struct SolarMatProps{T}
    ads::T     #adsorptance
    emitt::T    #emittance
end

# TODO: comment source
const Polished_beryllium = SolarMatProps(0.44, 0.01)
const Goldized_kapton = SolarMatProps(0.25, 0.02)
const Aluminum_tape = SolarMatProps(0.21, 0.04)
const Solar_cells_Si = SolarMatProps(0.75, 0.82)

export O, H, N, O2, N2, He, atmospheric_gases, Aluminum, Silicon, Polished_beryllium, Goldized_kapton, Aluminum_tape, Solar_cells_Si
export eta, Tw, s_cr, s_cd