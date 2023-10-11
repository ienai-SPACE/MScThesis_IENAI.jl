using LinearAlgebra

"""
    getTorques(coeffs_v, bc)

Aerodynamic torques wrt input reference frame
"""
function getTorques(coeffs_v, bc, ref_point::Vector)
    Cp_v = map(ii -> coeffs_v[ii][3], 1:lastindex(bc))

    Cp_x = [Cp_v[ii][1] for ii in 1:length(bc)]
    Cp_y = [Cp_v[ii][2] for ii in 1:length(bc)]
    Cp_z = [Cp_v[ii][3] for ii in 1:length(bc)]
    bc_x = [bc[ii][1] - ref_point[1] for ii in 1:lastindex(bc)]
    bc_y = [bc[ii][2] - ref_point[2] for ii in 1:lastindex(bc)]
    bc_z = [bc[ii][3] - ref_point[3] for ii in 1:lastindex(bc)]

    Mx = -sum((Cp_z .* bc_y - Cp_y .* bc_z))
    My = -sum((Cp_x .* bc_z - Cp_z .* bc_x))
    Mz = -sum((Cp_y .* bc_x - Cp_x .* bc_y))

    [Mx, My, Mz]
end

# """
# [Not valid] It returns a singular matrix because the CoP in 3D is a line, not a point
# """
# function getCoP(T, coeffs)
#     F = coeffs[1] * coeffs[2] + coeffs[3] * coeffs[4]
#     mat = [0 -F[3] F[2]; F[3] 0 -F[1]; -F[2] F[1] 0]
#     println("det(mat) = ", det(mat))
#     CoP = mat \ T

#     return CoP
# end

"""
    getCoP(T, coeffs, cpn::Vector)

This function calculates de center of pressure by applying the opposite of a cross product (F x u = -T) such that:
for a x b = c, b is not uniquely determined by a and c. Moreover, there is no solution unless a
and c are orthogonal. If a and c are orthogonal, then the solutions are (c×a)/(a.a)+ta
for arbitrary scalars t.
Source: https://www.goengineer.com/  and https://www.youtube.com/watch?v=GARurpX-VXE 

# Inputs
- `T`
- `coeffs`
- `cpn::Vector` : chord plane normal

# Outputs 
- `u1`
- `u2`  : CoP based on line of action and chord plane
"""
function getCoP(T, coeffs, cpn::Vector)

    F = coeffs[1] * coeffs[2] + coeffs[3] * coeffs[4]
    orthogonal_check = dot(F, T)
    println("orthogonal_check = ", orthogonal_check)
    if abs(orthogonal_check) < 1e-3
        FcrossT = cross(F, T)
        FdotF = dot(F, F)
        u1 = FcrossT / FdotF
        t = dot(u1, cpn) / dot(F, cpn)
        u2 = u1 - t * F
        [u1, u2]
    else
        println("F and T are not orthogonal, hece, no solution found")
    end
end
