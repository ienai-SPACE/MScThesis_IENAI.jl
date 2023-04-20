include("GSICalculations.jl")

""" function: CoefficientCalculations

GOAL:
    - Calculate the total drag, lift, pressure, and shear coefficients 
        - Contribution of each area element are area-weighted and added up for all contributing areas

INPUT:
    - m_srf  :: [g/mol]    molar mass of the satellite surface material
    - Afacet :: [m^2] area of the contributing element
    - δ      :: [rad] angle between the normal to the surface and the velocity vector

OUTPUT:
    - Total CD, CL, CP, CTAU
    """

mutable struct stLMNT{T}  #check this --> ::T returns the error: !Matched::T
    C::MVector{T}
    δ::Float64
    SRF::Float64
    PO::Float64
    Tw::Float64
end
struct stATMnDYN{T}  #check this --> ::T returns the error: !Matched::T
    Ta::T
    Vrel::T
end



function CoefficientCalculations(outMutableProps, outGasStreamProps, OutTriangles, Vrel)

    #---Pre-allocation--------------------------
    #m_srf = outMutableProps.m_srf
    #C = outGasStreamProps.C
    Afacet = OutTriangles[2, :]
    #δ = OutTriangles[3, :]

    Cd_facet = @MMatrix zeros(length(Afacet), 1)                                  #drag coefficient at each facet
    Cl_facet = @MMatrix zeros(length(Afacet), 1)                                  #lift coefficient at each facet
    Cp_facet = @MMatrix zeros(length(Afacet), 1)                                  #pressure coefficient at each facet
    Ctau_facet = @MMatrix zeros(length(Afacet), 1)                                #shear coefficient at each facet

    ATMnDYN = stATMnDYN(outGasStreamProps.Ta, Vrel)
    #--------------------------------------------

    A = sum(Afacet)

    for jj ∈ 1:Int(length(Afacet))

        LMNT = stLMNT(outGasStreamProps.C, OutTriangles[3, jj], outMutableProps.m_srf, outGasStreamProps.PO, outMutableProps.Tw)

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