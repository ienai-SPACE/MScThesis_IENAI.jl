import Base: +, -, /, *, zero

include("GSICalculations.jl")
include("SurfaceProps.jl")
include("IlluminationConvex.jl")


struct ElementInteractionProps{T}
    δ::T    # angle between oncoming direction vector and normal to the surface
    SRF::T  # molecular mass of the surface atom
    Tw::T   # temperature at the wall
end

#I DONT COMPLETELY UNDERSTAND WHAT THIS DOES! 
ElementInteractionProps(surfprops::SurfaceProps, angle) = ElementInteractionProps(angle, surfprops.m_srf, surfprops.Tw)

mutable struct stATMnDYN{T}
    Ta::T
    Vrel::T
    PO::T
    C::SVector{6,T}
end

struct InteractionGeometry{T}
    area::T
    angle::T
end

struct AerodynamicCoefficients{T}
    Cd::T
    Cl::T
    Cp::T
    Ctau::T
end


function drag_for_orientation_convex(geometry::Vector{<:Face}, gsp::GasStreamProperties, α, ϕ)

    Aproj, Aref, OutFacets = fIlluminationConvex(geometry, α, ϕ)
    InteractionGeometry(OutFacets[2, :], OutFacets[3, :])


    *(a::Real, c::AerodynamicCoefficients) = AerodynamicCoefficients(a * c.Cd, a * c.Cl, a * c.Cp, a * c.Ctau)
    *(c::AerodynamicCoefficients, a::Real) = a * c
    +(a::AerodynamicCoefficients, b::AerodynamicCoefficients) = AerodynamicCoefficients(a.Cd + b.Cd, a.Cl + b.Cl, a.Cp + b.Cp, a.Ctau + b.Ctau)
    /(a::AerodynamicCoefficients, d::Real) = AerodynamicCoefficients(a.Cd / d, a.Cl / d, a.Cp / d, a.Ctau / d)


    function compute_coefficients(surfprops::SurfaceProps, gasprops::gsp, intgeo::InteractionGeometry, Vrel_norm)
        element_interaction = ElementInteractionProps(surfprops, intgeo.angle)
        Cd, Cl, Cp, Ctau = DRIA_GSI(element_interaction, gasprops, Vrel_norm)
        AerodynamicCoefficients(Cd, Cl, Cp, Ctau)
    end

    function compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::Vector{<:InteractionGeometry}, Vrel_norm)
        aero_coeffs = map(intgeo) do intgeo
            compute_coefficients(surfprops, gasprops, intgeo, Vrel_norm)
        end
        areas = map(intgeo) do intgeo
            intgeo.area
        end
        num = sum(aero_coeffs .* areas)
        num / sum(areas), sum(areas)
    end

    return compute_coefficients

end


export InteractionGeometry, drag_for_orientation_convex