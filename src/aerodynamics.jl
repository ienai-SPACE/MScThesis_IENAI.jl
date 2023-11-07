using LinearAlgebra

"""
	GetCoP{T}

Center of pressure with respect to y-z, x-z, and x-y planes, respectively. Units: [m]
"""
struct GetCoP{T}
	CoPx::T
	CoPy::T
	CoPz::T
end

"""
	getCoP(M::SVector, coeffs::Tuple, A)

Calculation of center of pressure at y-z, x-z, and x-y planes

# Inputs:
- `M::SVector`
- `coeffs::Tuple`
- `A`			: projected cross-sectional area

# Outputs:
- `GetCoP`		: [m] (CoPx, CoPy, CoPz)	: CoP
- `GetCoP`		: [m] (CoPx, CoPy, CoPz)	: point on the line of action of forces
"""
function getCoP(M::SVector, coeffs::Tuple, A)
	cpn_x = SVector(1.0, 0.0, 0.0)  # chord plane normal
	cpn_y = SVector(0.0, 1.0, 0.0)
	cpn_z = SVector(0.0, 0.0, 1.0)
	CoP_u1u2 = [_getCoP(M, coeffs, A, cpn_x), _getCoP(M, coeffs, A, cpn_y), _getCoP(M, coeffs, A, cpn_z)]
	GetCoP(SVector(CoP_u1u2[1][2]), SVector(CoP_u1u2[2][2]), SVector(CoP_u1u2[3][2])), GetCoP(SVector(CoP_u1u2[1][1]), SVector(CoP_u1u2[2][1]), SVector(CoP_u1u2[3][1]))
end

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
function _getCoP(T, coeffs::Tuple, A, cpn::SVector)

	# F = (coeffs[1] * coeffs[2] + coeffs[3] * coeffs[4]) * A #(Cd + Cl)*Aproj = Fnet
	F = (coeffs[5] * coeffs[6]) * A

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
	CoP = (1 / PCSA) * sum([intercept_info[ii].area * cos(intercept_info[ii].angle) for ii in 1:lastindex(barycenters)] .* barycenters, dims = 1)
	SVector(CoP[1][1], CoP[1][2], CoP[1][3])
end


"""
	getTorques(coeffs_v::Vector{Vector{AbstractVector{Float64}}}, intercept_info::Vector{<:InteractionGeometry}, ref_point::SVector{3,Float64}, ray_coords::Vector{SatelliteGeometryCalculations.IntersectCoords{Float64,Int64}}, Nray_idx, hit_idx, f_hit_idx)

Non-dimensional aerodynamic torques wrt input reference frame. Normalized with aerodynamic pressure.

# Inputs
- `coeffs_v::Vector{Vector{AbstractVector{Float64}}}`
- `intercept_info::Vector{<:InteractionGeometry}`
- `ref_point::SVector{3,Float64}`
- `ray_coords::Vector{SatelliteGeometryCalculations.IntersectCoords{Float64,Int64}}`
-  `Nray_idx::Vector{Int64}`	: number of rays striking each impinged facet
- `hit_idx::Vector{Int64}`		: sorted vector which stores the indices of all impinged facets
- `f_hit_idx::Vector{Int64}`	: vector of length equal to the number of striking rays in which each position stores the index of the impinged facet

# Outputs
- `::SVector{3, Float64}`               : coefficients of torques (M = CM*A/(0.5*rho*v^2)) i.e. CMx, CMy, CMz
"""
function getTorques(
	coeffs_v::Vector{Vector{AbstractVector{Float64}}},
	intercept_info::Vector{<:InteractionGeometry},
	ref_point::SVector{3, Float64},
	rf::SatelliteGeometryCalculations.Ray_Facet_info,
)

	areas = [intercept_info[ii].area * cos(intercept_info[ii].angle) for ii in 1:lastindex(intercept_info)]
	Cnet_v = map(ii -> coeffs_v[ii][1] + coeffs_v[ii][2], 1:lastindex(intercept_info)) .* areas #(Cd + Cl)*Aproj

	#Aerodynamic coefficient * projected area of each of the intercepted facets
	Cnet_x = [Cnet_v[ii][1] for ii in 1:lastindex(intercept_info)]
	Cnet_y = [Cnet_v[ii][2] for ii in 1:lastindex(intercept_info)]
	Cnet_z = [Cnet_v[ii][3] for ii in 1:lastindex(intercept_info)]

	#distance between striking points of the rays and the reference point
	bc_x = [rf.ray_coords[ii].coords[1] - ref_point[1] for ii in 1:lastindex(rf.ray_coords)]
	bc_y = [rf.ray_coords[ii].coords[2] - ref_point[2] for ii in 1:lastindex(rf.ray_coords)]
	bc_z = [rf.ray_coords[ii].coords[3] - ref_point[3] for ii in 1:lastindex(rf.ray_coords)]

	CMx = 0
	CMy = 0
	CMz = 0
	# jj_v = Vector{Int64}(undef,lastindex(ray_coords))
	for ii in 1:lastindex(rf.ray_coords)
		jj = findfirst(x -> x == rf.f_hit_idx[ii], rf.hit_idx)

		CMx_ = (Cnet_z[jj] / rf.Nray_idx[jj]) * bc_y[ii] - (Cnet_y[jj] / rf.Nray_idx[jj]) * bc_z[ii]
		CMy_ = (Cnet_x[jj] / rf.Nray_idx[jj]) * bc_z[ii] - (Cnet_z[jj] / rf.Nray_idx[jj]) * bc_x[ii]
		CMz_ = (Cnet_y[jj] / rf.Nray_idx[jj]) * bc_x[ii] - (Cnet_x[jj] / rf.Nray_idx[jj]) * bc_y[ii]
		CMx = CMx - CMx_
		CMy = CMy - CMy_
		CMz = CMz - CMz_
		# jj_v[ii] = jj
	end

	SVector(CMx, CMy, CMz)
end

"""
	getAeroTorque(coeffs::Tuple, CoP::Vector, CoM::Vector, A, rho, V)

Calculation of the aerodynamic torque
"""
function getAeroTorque(coeffs::Tuple, CoP::SVector, CoM::SVector, A, rho, V)
	u_v = V / norm(V)
	d = CoM - CoP
	0.5 * rho * coeffs[1] * A * norm(V)^2 * cross(u_v, d)
end

"""
	getAeroTorque(coeffs::Tuple, CoP::Vector, CoM::Vector, A)

Calculation of the aerodynamic torque, normalized with the dynamic pressure
"""
function getAeroTorque(coeffs::Tuple, CoP::SVector, CoM::SVector, A)
	u_v = -coeffs[2]
	println("u_v = ", u_v)
	println("typeof(CoP)= ", typeof(CoP))
	println("typeof(CoM) = ", typeof(CoM))
	d = CoM - CoP
	coeffs[1] * A * cross(u_v, d)
end
