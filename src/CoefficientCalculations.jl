include("GSICalculations.jl")

""" function: CoefficientCalculations

GOAL:
    - Calculate the total drag, lift, pressure, and shear coefficients 
        - Contribution of each area element are area-weighted and added up for all contributing areas
INPUT:
    - outMutableProps       : struct with mutable properties
    - outGasStreamProps     : struct with gas stream properties
    - OutLMNTs              : [index of the element, [m^2] area of the element, [rad] angle between oncoming direction vector and normal to the surface]
    - Vrel_norm             : [m/s] magnitude of the relatice velocity
OUTPUT:
    - Total CD, CL, CP, CTAU
    """

struct stLMNT{T}  #check this --> ::T returns the error: !Matched::T
    δ::T
    SRF::T
    Tw::T
end
mutable struct stATMnDYN{T}  #check this --> ::T returns the error: !Matched::T
    Ta::T
    Vrel::T
    PO::T
    C::MVector{6,T}
end



function CoefficientCalculations(outMutableProps, outGasStreamProps, OutLMNTs, Vrel_norm)

    #---Pre-allocation--------------------------
    #m_srf = outMutableProps.m_srf
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

        LMNT = stLMNT(OutLMNTs[3, jj], outMutableProps.m_srf, outMutableProps.Tw)

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