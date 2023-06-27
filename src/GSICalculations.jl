using SpecialFunctions
using StaticArrays
"""
    DRIA_GSI(element_interaction, gasprops, Vrel)

Reference: [DOORNBOS 2012, Thermospheric Density and Wind Determination from Satellite Dynamics] + modification based on atomic oxygen abosption=#

#INPUTS:
- `element_interaction::ElementInteractionProps{T}`       : struct with `δ`, `SRF`, and `Tw` as fields
- `gasprops::GasStreamProperties`                         : mass of the atoms of the surface particles
- `Vrel:Vector`    

#OUTPUTS:
-`Cd_facet`     : Drag Coefficient
-`Cl_facet`     : Lift Coefficient
-`Cp_facet`     : Pressure Coefficient
-`Ctau_facet`   : Shear Coefficient
"""

function DRIA_GSI(element_interaction, gasprops, Vrel, normals)

    δ = element_interaction.δ
    C = gasprops.C
    m_srf = element_interaction.SRF
    P0 = gasprops.PO
    Tw = element_interaction.Tw
    Ta = gasprops.Ta
    Vrel_norm = norm(Vrel)



    K = 1.44e6                                  #Best-fit Langmuir adsorbate constants for DRIA GSI model
    θ = K * P0 / (1 + K * P0)                   #Fraction of the surface contaminated by atomic oxygen


    γ = cos(δ)
    ell = sin(δ)



    #-------Pre-allocation---------------------------------------------
    mgas = @SVector [O.m, H.m, N.m, O2.m, N2.m, He.m]                   #[g/mol] atomic mass of gas constituents
    #mgas = @SVector [m[1].m, m[2].m, m[3].m, m[4].m, m[5].m, m[6].m] #[g/mol] atomic mass of gas constituents
    cd_j = @MMatrix zeros(length(mgas), 2)                               #drag coefficient at each facet
    cl_j = @MMatrix zeros(length(mgas), 2)                               #lift coefficient at each facet
    #s = @MArray [zeros(length(mgas), 2)] 
    s = @MMatrix zeros(length(mgas), 2)                                  #thermal speed: two values per species, i.e. contaminated and clean surface
    #------------------------------------------------------------------


    for jj ∈ 1:2     #clean and contaminated surface
        for j ∈ 1:6  #species specific mass concentration

            μ_srf = mgas[j] / m_srf                 #ratio between the mass of the atoms of the incoming gas with the mass of the surface particles
            Ks = 3.6                                #substrate coefficient (6 < s < 11 the use of Ks = 3.6 is appropriate)
            α_c = Ks * μ_srf / (1 + μ_srf)^2        #accomodation coefficient for clean surface

            α_vec = [α_c, 1]

            α = α_vec[jj]



            s[j, jj] = Vrel_norm / sqrt(2 * (kb / (mgas[j] / NA / 1000) * Ta))    #thermal speed
            P = exp.(-γ .^ 2 .* s[j, jj] .^ 2) ./ s[j, jj]
            G = 1 / (2 * s[j, jj] .^ 2)
            Q = 1 + G
            Z = 1 + erf.(γ .* s[j, jj])
            RR = R / mgas[j] * 1000
            Vratio = sqrt((1 / 2) * (1 + α * ((4 * RR * Tw) / Vrel_norm^2 - 1)))  #[Doornbos 2012]
            cd_j[j, jj] = (P ./ sqrt(π) .+ γ .* Q .* Z .+ (0.5 * γ) .* Vratio .* (γ .* sqrt(pi) .* Z .+ P)) .* C[j] * mgas[j]
            cl_j[j, jj] = (ell * G .* Z .+ (0.5 .* ell) * Vratio .* (γ .* sqrt(pi) .* Z .+ P)) .* C[j] * mgas[j]
        end
    end



    Cd_weighted = sum(cd_j, dims=1) / sum(C .* mgas)         # mass weighted
    Cl_weighted = sum(cl_j, dims=1) / sum(C .* mgas)         # mass weighted

    Cd_weighted_clean = Cd_weighted[1, 1]                    # clean (α ~= 1)
    Cl_weighted_clean = Cl_weighted[1, 1]                    # clean (α ~= 1)
    Cd_weighted_cont = Cd_weighted[1, 2]                     # contaminated (α = 1)
    Cl_weighted_cont = Cl_weighted[1, 2]                     # contaminated (α = 1)


    Cd_facet = (Cd_weighted_clean * (1 - θ) + Cd_weighted_cont * θ)
    Cl_facet = (Cl_weighted_clean * (1 - θ) + Cl_weighted_cont * θ)

    Cp_facet = Cd_facet * cos(δ) + Cl_facet * sin(δ)
    Ctau_facet = Cd_facet * sin(δ) - Cl_facet * cos(δ)





    return Cd_facet, Cl_facet, Cp_facet, Ctau_facet

end
