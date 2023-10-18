using LinearAlgebra

"""
    getTorques(coeffs_v, bc)

Non-dimensional aerodynamic torques wrt input reference frame. Normalized with aerodynamic pressure.

# Inputs
- `coeffs_v::Vector{Vector{AbstractVector{Float64}}}`
- `bc::Vector{Vector}`
- `ref_point::SVector{3,Float64}`
# Outputs
- `::SVector{3, Float64}`   : coefficients of torques (M = CM/(0.5*rho*v^2)) i.e. CMx, CMy, CMz
"""
function getTorques(coeffs_v::Vector{Vector{AbstractVector{Float64}}}, intercept_info::Vector{<:InteractionGeometry}, bc::Vector{Vector}, ref_point::SVector{3,Float64})
    Cp_v = map(ii -> coeffs_v[ii][3], 1:lastindex(bc))

    Cp_x = [Cp_v[ii][1] for ii in 1:length(bc)]
    Cp_y = [Cp_v[ii][2] for ii in 1:length(bc)]
    Cp_z = [Cp_v[ii][3] for ii in 1:length(bc)]
    bc_x = [bc[ii][1] - ref_point[1] for ii in 1:lastindex(bc)]
    bc_y = [bc[ii][2] - ref_point[2] for ii in 1:lastindex(bc)]
    bc_z = [bc[ii][3] - ref_point[3] for ii in 1:lastindex(bc)]
    areas = [intercept_info[ii].area for ii in 1:lastindex(bc)]

    CMx = -sum((Cp_z .* bc_y - Cp_y .* bc_z) .* areas)
    CMy = -sum((Cp_x .* bc_z - Cp_z .* bc_x) .* areas)
    CMz = -sum((Cp_y .* bc_x - Cp_x .* bc_y) .* areas)

    SVector(CMx, CMy, CMz)
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
Source: https://www.goengineer.com/blog/calculating-center-of-pressure-solidworks-flow-simulation#:~:text=To%20do%20this%2C%20create%20a%203D%20sketch%2C%20add,the%20chord%20plane%20is%20the%20center%20of%20pressure.
  and https://www.youtube.com/watch?v=GARurpX-VXE 

# Inputs
- `T`               : coefficient of torque (non-dimensionalized with the dynamic pressure)
- `coeffs::Tuple`   : aerodynamic coefficients
- `A`               : projected cross-sectional area
- `cpn::Vector`     : chord plane normal

# Outputs 
- `u1::SVector{3, Float64}`     : [m]
- `u2::SVector{3, Float64}`     : [m] CoP based on line of action and chord plane
"""
function getCoP(T, coeffs::Tuple, A, cpn::SVector)

    F = (coeffs[1] * coeffs[2] + coeffs[3] * coeffs[4]) * A #(Cd + Cl)*Aproj
    orthogonal_check = dot(F, T)
    println("orthogonal_check: dot(F, T) = ", orthogonal_check)
    if abs(orthogonal_check) < 1e-3  # F and T are orthogonal
        FcrossT = cross(F, T)
        FdotF = dot(F, F)
        u1 = SVector{3}(FcrossT / FdotF)
        t = dot(u1, cpn) / dot(F, cpn)
        u2 = SVector{3}(u1 - t * F) #CoP
        [u1, u2]
    else
        println("F and T are not orthogonal, hece, no solution found")
    end
end


"""
    getCoP(PCSA, intercept_info, barycenters)

Calculation of the center of pressure wrt the local ref.frame at which the vertices of the geometry are defined
Note: cannot be used for GSI because the aerodyn. coeff. also depends on the angle of incidence, not just the area

# Inputs
- `PCSA` : projected cross-sectional area = Aproj
- `intercept_info::Vector{InteractionGeometryHomo{Float64}}` or `Vector{InteractionGeometry{Float64}}`
- `barycenters::Vector{Vectors}`

# Outputs
- `CoP::SVector{3,Float64}`
"""
function getCoP(PCSA, intercept_info::Vector{<:InteractionGeometry}, barycenters)
    println("length(barycenters) = ", length(barycenters))
    println("length(intercept_info) = ", length(intercept_info))
    CoP = (1 / PCSA) * sum([intercept_info[ii].area * cos(intercept_info[ii].angle) for ii in 1:lastindex(barycenters)] .* barycenters, dims=1)
    SVector(CoP[1][1], CoP[1][2], CoP[1][3])
end

"""
    getTa(coeffs::Tuple, CoP::Vector, CoM::Vector, A, rho, V)

Calculation of the aerodynamic torque
"""
function getTa(coeffs::Tuple, CoP::Vector, CoM::Vector, A, rho, V)
    u_v = V / norm(V)
    d = CoM - CoP
    0.5 * rho * coeffs[1] * A * norm(V)^2 * cross(u_v, d)
end