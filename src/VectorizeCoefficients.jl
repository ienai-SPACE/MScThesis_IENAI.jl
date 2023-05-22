include("unitaryDirections.jl")

"""
    vectorizeCoeffs(aero_coeffs, normals, Vrel_v)

Provide a vectorial direction to all `coeffs_vec = coefficients.*area` and find the resultant direction
"""

function vectorizeCoeffs(aero_coeffs, normals, Vrel_v)

    coeffs_vec = [zeros(3), zeros(3), zeros(3), zeros(3)]

    for ii ∈ 1:Int(size(aero_coeffs, 1))

        u_D, u_L, u_P, u_tau, = unitaryDirections(Vrel_v, normals[ii])

        coeffs_vec_ii = [aero_coeffs[ii].Cd * u_D, aero_coeffs[ii].Cl * u_L, aero_coeffs[ii].Cp * u_P, aero_coeffs[ii].Ctau * u_tau]
        coeffs_vec += coeffs_vec_ii

    end

    return coeffs_vec
end