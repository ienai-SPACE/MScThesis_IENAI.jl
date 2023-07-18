"""
    Viewpoint{T}
- `rmax::T`
- `distance::T`
- `direction::SV3{T}`  
"""

struct Viewpoint{T}
    rmax::T
    distance::T
    direction::SV3{T}
end

"""
    Viewpoint(rmax::T, distance::T, azimuth::T, elevation::T) where {T}

Storage of: viewpoint information (azimuth and elevation wrt to local reference frame), maximum Euclidean distance (`rmax`), and `distance` to a perpendicular plane

#INPUTS
- `rmax::T`
- `distance::T`
- `azimuth::T`
- `elevation::T`
#OUTPUTS
- `Viewpoint{T}`  : fields are `rmax`, `distance`, `dir`
"""

function Viewpoint(rmax::T, distance::T, azimuth::T, elevation::T) where {T}
    α = azimuth
    ϕ = elevation
    ux = cos(ϕ) * cos(α)
    uy = cos(ϕ) * sin(α)
    uz = sin(ϕ)
    dir = SV3([-ux, -uy, -uz])
    Viewpoint{T}(rmax, distance, dir)
end

"""
    Viewpoint(rmax::T, distance::T, azimuth::T, elevation::T) where {T}

Storage of: viewpoint information (azimuth and elevation wrt to local reference frame), maximum Euclidean distance (`rmax`), and `distance` to a perpendicular plane

#INPUTS
- `rmax::T`
- `distance::T`
- `Vrel::Vector{T}`
#OUTPUTS
- `Viewpoint{T}`        : fields are `rmax`, `distance`, `dir`
"""

function Viewpoint(rmax::T, distance::T, Vrel::Vector{T}) where {T}
    Vdir = normalize(Vrel)
    dir = SV3(-Vdir)
    Viewpoint{T}(rmax, distance, dir)
end

export Viewpoint