struct Viewpoint{T}
    rmax::T
    distance::T
    direction::SV3{T}
end

function Viewpoint(rmax::T, distance::T, azimuth::T, elevation::T) where {T}
    α = azimuth
    ϕ = elevation
    ux = cos(ϕ) * cos(α)
    uy = cos(ϕ) * sin(α)
    uz = sin(ϕ)
    dir = SV3([-ux, -uy, -uz])
    Viewpoint{T}(rmax, distance, dir)
end

export Viewpoint