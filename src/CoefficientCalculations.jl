import Base: +, -, /, *, zero

include("GSICalculations.jl")

""" 
    `compute_coefficients`
    GOAL: compute the toal drag, lift, pressure, and shear coefficients evaluated for a single area element
    INPUT:
    - surfprops::SurfaceProps
    - gasprops::GasStreamProperties
    - intgeo::InteractionGeometry
    - Vrel_norm : T     [m/s]

    `function compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::Vector{<:InteractionGeometry}, Vrel_norm)`
    GOAL: compute the toal drag, lift, pressure, and shear coefficients evaluated for all areas of the element
    INPUT:
    - surfprops::SurfaceProps
    - gasprops::GasStreamProperties
    - intgeo::Vector{<:InteractionGeometry}
    - Vrel_norm : T
    OUTPUT:
    - Total CD, CL, CP, CTAU, Aref
    """

struct ElementInteractionProps{T}
    δ::T    # angle between oncoming direction vector and normal to the surface
    SRF::T  # molecular mass of the surface atom
    Tw::T   # temperature at the wall
end

#check potential mistake here: ElementInteractionProps(angle, surfprops::SurfaceProps) 
ElementInteractionProps(surfprops::SurfaceProps, angle) = ElementInteractionProps(angle, surfprops.m_srf, surfprops.Tw)

# mutable struct stATMnDYN{T}
#     Ta::T
#     Vrel::T
#     PO::T
#     C::SVector{6,T}
# end

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
    areas = map(intgeo) do x
        x.area
    end
    angles = map(intgeo) do x
        x.angle
    end
    num = sum(aero_coeffs .* areas)

    # return scaled_coefficients = sum(C*A_i)/Aref, A_tot
    num / sum(areas .* cos.(angles)), sum(areas), sum(areas .* cos.(angles))


end


#EXAMPLES
# julia> map(x->x^2, 1:6)
# 6-element Vector{Int64}:
#   1
#   4
#   9
#  16
#  25
#  36

# julia> squares = map(1:6) do x
#            x^2
#        end
# 6-element Vector{Int64}:
#   1
#   4
#   9
#  16
#  25
#  36



#=
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
=#
export InteractionGeometry, compute_coefficients