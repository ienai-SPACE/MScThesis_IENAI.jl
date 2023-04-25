import Base: +, -, /, *, zero

include("GSICalculations.jl")

""" function: CoefficientCalculations

GOAL:
    - Calculate the total drag, lift, pressure, and shear coefficients 
        - Contribution of each area element are area-weighted and added up for all contributing areas
INPUT:
    - outSurfaceProps       : struct with mutable properties
    - outGasStreamProps     : struct with gas stream properties
    - OutLMNTs              : [index of the element, [m^2] area of the element, [rad] angle between oncoming direction vector and normal to the surface]
    - Vrel_norm             : [m/s] magnitude of the relatice velocity
OUTPUT:
    - Total CD, CL, CP, CTAU
    """

struct ElementInteractionProps{T}  #check this --> ::T returns the error: !Matched::T
    δ::T    # angle between oncoming direction vector and normal to the surface
    SRF::T  # molecular mass of the surface atom
    Tw::T   # temperature at the wall
end

ElementInteractionProps(surfprops::SurfaceProps, angle) = ElementInteractionProps(angle, surfprops.m_srf, surfprops.Tw)

mutable struct stATMnDYN{T}  #check this --> ::T returns the error: !Matched::T
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

*(a::Real, c::AerodynamicCoefficients) = AerodynamicCoefficients(a * c.Cd, a * c.Cl, a * c.Cp, a * c.Ctau)
*(c::AerodynamicCoefficients, a::Real) = a * c
+(a::AerodynamicCoefficients, b::AerodynamicCoefficients) = AerodynamicCoefficients(a.Cd + b.Cd, a.Cl + b.Cl, a.Cp + b.Cp, a.Ctau + b.Ctau)
/(a::AerodynamicCoefficients, d::Real) = AerodynamicCoefficients(a.Cd / d, a.Cl / d, a.Cp / d, a.Ctau / d)

function compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::InteractionGeometry, Vrel_norm)
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


function CoefficientCalculations(outSurfaceProps, outGasStreamProps, OutLMNTs, Vrel_norm)

    #---Pre-allocation--------------------------
    #m_srf = outSurfaceProps.m_srf
    #C = outGasStreamProps.C
    Afacet = OutLMNTs[2, :]
    #δ = OutTriangles[3, :]

    Cd_facet = @MMatrix zeros(length(Afacet), 1)                                  #drag coefficient at each facet
    Cl_facet = @MMatrix zeros(length(Afacet), 1)                                  #lift coefficient at each facet
    Cp_facet = @MMatrix zeros(length(Afacet), 1)                                  #pressure coefficient at each facet
    Ctau_facet = @MMatrix zeros(length(Afacet), 1)                                #shear coefficient at each facet

    ATMnDYN = stATMnDYN(outGasStreamProps.Ta, Vrel_norm, outGasStreamProps.PO, outGasStreamProps.C)
    #--------------------------------------------

    A = sum(Afacet)

    for jj ∈ 1:Int(length(Afacet))

        LMNT = ElementInteractionProps(OutLMNTs[3, jj], outSurfaceProps.m_srf, outSurfaceProps.Tw)

        #Cd, Cl, Cp, Ctau = DRIA_GSI(δ[jj], C, m_srf)
        Cd, Cl, Cp, Ctau = DRIA_GSI(LMNT, ATMnDYN)

        Cd_facet[jj] = Cd * Afacet[jj] / A
        Cl_facet[jj] = Cl * Afacet[jj] / A
        Cp_facet[jj] = Cp * Afacet[jj] / A
        Ctau_facet[jj] = Ctau * Afacet[jj] / A

    end

    CD = sum(Cd_facet)
    CL = sum(Cl_facet)
    CP = sum(Cp_facet)
    CTAU = sum(Ctau_facet)

    return CD, CL, CP, CTAU

end

export InteractionGeometry, compute_coefficients