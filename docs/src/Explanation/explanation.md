# General overview

## Functionalities of the code:

* Calculation of PCSA (projected cross-sectional area).
* Calculation of aerodynamic coefficients.
* 360-sweep covering both of the above funtionalities. The output will be a look-up table that can be accessed with an interpolator.

## Structure of the code:

The code may be divided into 4 core sections:
1. The environment
2. Interaction with the environment
3. Area calculation
4. Numerical models

Firstly, the environmental conditions where the satellite intends to operate should be addressed. For this, the  [`SatelliteToolbox`](@ref) package is used, utilizing the MSISE00 model [`SatelliteToolbox.NRLMSISE00_Output`](@ref). 

Secondly, the interactions taking place between the satellite and the gas particles of the environment; gas-surface interactions
are calculated with the implementation of the DRIA model (diffuse reflection with incomplete accommodation) [`SatelliteGeometryCalculations.compute_coefficients`](@ref) [`SatelliteGeometryCalculations.DRIA_GSI`](@ref). 

The third and fourth sections are inherently related because while the calculation of areas can be done in many different ways
and without prior inputs, the selection of a specific numerical model sets preferences over the way the areas are calculated. For the most general case, the **non-convex case** [`SatelliteGeometryCalculations.areas_nonconvex`](@ref) case, a ray-tracing method is implemented based on a Möller-Trumbore algorithm ([`SatelliteGeometryCalculations.MTalgorithm`](@ref)). It allows the identification of the intercepted facets and the counting of the total amount of rays that strike; each ray has an equivalent area, which after projection based on the viewpoint yields the **projected cross-sectional area** and **total reference area**. For the **convex case** [`SatelliteGeometryCalculations.areas_convex`](@ref), the forward facing areas are identified and projected onto the desired viewpoint [`SatelliteGeometryCalculations.Viewpoint`](@ref).