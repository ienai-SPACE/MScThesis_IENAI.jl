# Torques and center of pressure calculation

The calculation of the center of pressure and aerodynamic torques is done from the aerodynamic coefficient information. [`SatelliteGeometryCalculations.getCoP`](@ref) and [`SatelliteGeometryCalculations.getTorques`](@ref) implement these calculations.

## Torque from aerodynamic forces

The output of [`SatelliteGeometryCalculations.getTorques`](@ref) is the coefficients of torques, normalized with the dynamic pressure: ``CM = M/(0.5*\rho*v^2)``. The balance of moments is evaluated for each contributing triangular facet, and the sum of all contributions provides the total torque.

$CMx = -\sum(Cp_z * bc_y - Cp_y * bc_z) * A$
$CMy = -\sum(Cp_x * bc_z - Cp_z * bc_x) * A$
$CMz = -\sum(Cp_y * bc_x - Cp_x * bc_y) * A$

where ``Cp_k`` corresponds to the pressure coefficient component in the ``k-``direction, and ``bc_k`` stands for the distance between the barycenter of a facet to the `ref_point` in the ``k-``direction.

The torque is calculated with respect to a reference point `ref_point`, which should be zero when computed at the center of pressure (CoP). And when computed about the center of mass (CoM), it provides the aerodynamic moment.

```
getTorques(coeffs_v::Vector{Vector{AbstractVector{Float64}}}, intercept_info::Vector{<:InteractionGeometry}, bc::Vector{Vector}, ref_point::SVector{3,Float64})
```

## Center of pressure

In 3D, the **center of pressure** corresponds to a line. This is found by applying the opposite of a cross product (F x u = -T) where F is the force, u the distance to the reference point, and T the resulting torque, such that: for a x b = c, b is not uniquely determined by a and c. Moreover, there is no solution unless a and c are orthogonal. If a and c are orthogonal, then the solutions are ``(c×a)/(a.a)+ta`` for arbitrary scalars t. 

Source: [Calculating the Center of Pressure in SOLIDWORKS Flow Simulation](https://www.goengineer.com/blog/calculating-center-of-pressure-solidworks-flow-simulation#:~:text=To%20do%20this%2C%20create%20a%203D%20sketch%2C%20add,the%20chord%20plane%20is%20the%20center%20of%20pressure.)
and [GoEngineering Youtube channel](https://www.youtube.com/watch?v=GARurpX-VXE) 

[`SatelliteGeometryCalculations.getCoP`](@ref) retrieves the CoP with respect to the `ref_point` about which moments are calculated. The function requires the torques (normalized with the dynamic pressure), the aerodynamic coefficients, the projected cross-sectional area, and the normal of the chord plane. The output are two points (u1 and u2) provided in [m] which correspond to points on the line of action of the force. 

```
CoP = SatelliteGeometryCalculations.getCoP(CT, coeffs, Aproj, chord_plane_n)
```

u2 (CoP[2]) is the CoP based on line of action and chord plane.