""" Different Gas-Surface-Interaction models """


using SpecialFunctions


function DRIA_GSI(Afacet, δ, C, m_srf)

    #DRIA model#

    #=[DOORNBOS 2012, Thermospheric Density and Wind Determination from Satellite Dynamics] + modification based on atomic oxygen abosption=#

    #INPUTS:
    #----------------------------------------------------------------------
    #m_srf  :: mass of the atoms of the surface particles
    #C      :: gas concentrations
    #Afacet :: area of each facet
    # δ     :: incoming particle angle wrt surface's normal
    #----------------------------------------------------------------------
    #OUTPUTS:
    #----------------------------------------------------------------------
    #CD     :: Total Drag Coefficient
    #CL     :: Total Lift Coefficient
    #CP     :: Total Pressure Coefficient
    #CTAU   :: Total Shear Coefficient
    #----------------------------------------------------------------------


    K = 1.44e6                                  #Best-fit Langmuir adsorbate constants for DRIA GSI model
    θ = K * P0 / (1 + K * P0)                   #Fraction of the surface contaminated by atomic oxygen


    A = sum(Afacet)

    γ = cos.(δ)
    ell = sin.(δ)



    #Pre-allocation:
    s = zeros(6, 2)
    cd_j = zeros(length(δ), 6, 2)
    cl_j = zeros(length(δ), 6, 2)
    mgas = [m[1].m, m[2].m, m[3].m, m[4].m, m[5].m, m[6].m] #[g/mol] atomic mass of gas constituents


    for jj ∈ 1:2     # clean and contaminated surface
        for j ∈ 1:6 # Species specific mass concentration

            μ_srf = mgas[j] / m_srf                 #ratio between the mass of the atoms of the incoming gas with the mass of the surface particles
            Ks = 3.6                              #substrate coefficient (6 < s < 11 the use of Ks = 3.6 is appropriate)
            α_c = Ks * μ_srf / (1 + μ_srf)^2            #accomodation coefficient for clean surface

            α_vec = [α_c, 1]

            α = α_vec[jj]



            s[j, jj] = Vrel / sqrt(2 * (kb / (mgas[j] / NA / 1000) * Ta)) #thermal speed
            P = exp.(-γ .^ 2 .* s[j, jj] .^ 2) ./ s[j, jj]
            G = 1 / (2 * s[j, jj] .^ 2)
            Q = 1 + G
            Z = ones(size(γ)) + erf.(γ .* s[j, jj])
            RR = R / mgas[j] * 1000
            Vratio = sqrt((1 / 2) * (1 + α * ((4 * RR * Tw) / Vrel^2 - 1))) # [Doornbos 2012]
            cd_j[:, j, jj] = (P ./ sqrt(π) .+ γ .* Q .* Z .+ (0.5 * γ) .* Vratio .* (γ .* sqrt(pi) .* Z .+ P)) .* C[j] * mgas[j]
            cl_j[:, j, jj] = (ell * G .* Z .+ (0.5 .* ell) * Vratio .* (γ .* sqrt(pi) .* Z .+ P)) .* C[j] * mgas[j]
        end
    end



    Cd_weighted = sum(cd_j, dims=2) / sum(C .* mgas)         # mass weighted
    Cl_weighted = sum(cl_j, dims=2) / sum(C .* mgas)         # mass weighted

    Cd_weighted_clean = Cd_weighted[:, 1, 1]                 # clean (α ~= 1)
    Cl_weighted_clean = Cl_weighted[:, 1, 1]                 # clean (α ~= 1)
    Cd_weighted_cont = Cd_weighted[:, 1, 2]                  # contaminated (α = 1)
    Cl_weighted_cont = Cl_weighted[:, 1, 2]                  # contaminated (α = 1)


    Cd_facet = (Cd_weighted_clean * (1 - θ) + Cd_weighted_cont * θ) .* Afacet ./ A
    Cl_facet = (Cl_weighted_clean * (1 - θ) + Cl_weighted_cont * θ) .* Afacet ./ A

    Cp_facet = Cd_facet .* cos.(δ) + Cl_facet .* sin.(δ)
    Ctau_facet = Cd_facet .* sin.(δ) - Cl_facet .* cos.(δ)

    CD = sum(Cd_facet)
    CL = sum(Cl_facet)
    CP = sum(Cp_facet)
    CTAU = sum(Ctau_facet)

    return CD, CL, CP, CTAU

end



























#= FOR A FLAT PLATE: NO Cd AND Cl DECOMPOSITION
-------------------------------------------------------------------------------------------------------------------
Ks = 3.6;                            #substrate coefficient (6 < s < 11 the use of Ks = 3.6 is appropriate)
theta = K*P0/(1+K*P0);               #fraction of the surface contaminated by atomic oxygen

#--- loop from 1:6 (gas species)
mu_srf = gas_i/m_srf;                #ratio between the mass of the atoms of the incoming gas with the mass of the surface particles
alpha = Ks*mu_srf/(1+mu_srf)^2;      #for clean surface
alpha_d = 1;                         #for contaminated surface (completely diffuse) 
erf_s_d = 0 ;                #pending: implement error function

Tki_D = (2/3)*Ta*s_i^2;
Tkr_d = Tki*(1-alpha) + alpha*Tw;

cd_d_i = ((2+1/s_i^2)*erf_s + 2/(sqrt(pi)*s_i)*exp(-s_i^2) + sqrt(Tkr/Ta)*sqrt(pi)/s);   #perpendicular drag for each gas species for a contaminated surface
cd_i = (2+1/s_i^2)*erf_s + 2/(sqrt(pi)*s_i)*exp(-s_i^2) + sqrt(Tkr/Ta)*sqrt(pi)/s;     #perpendicular drag for each gas species for a clean surface
#--- end of loop

cd_d = dot(cd_d_i,m_i)/m_tot;        #weighted averaged drag (contaminated)
cd = dot(cd_i,m_i)/m_tot;            #weighted averaged drag (clean)


CD_k = theta*cd_d + (1-theta)*cd;     #for each facet

CD = dot(CD_k*Afacet)/Nfacet;         #total drag         
=#