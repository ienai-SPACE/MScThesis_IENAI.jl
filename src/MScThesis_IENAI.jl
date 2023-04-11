module MScThesis_IENAI

using LinearAlgebra
# Write your package code here.



include("PermanentPropertiesNew.jl")
include("EnvironmentalInputs.jl")
include("GSICalculations.jl")


m_srf = Aluminum.amass;              #[g/mol] atomic mass of the surface atom
Afacet = 1.5.*ones(3,1);             #[m^2]   area of each facet (set as input)
δ = [0, 0, pi/2];                    #[rad]   angle between the normal of the surface and the incoming particle   


Cd, Cl, Cp, Ctau = DRIA_GSI(Afacet,δ,C,m_srf)

#using PermanentPropertiesNew
end