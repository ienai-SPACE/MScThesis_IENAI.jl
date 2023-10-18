# General overview

## Functionalities of the code:

* Calculation of PCSA (projected cross-sectional area).
* Calculation of aerodynamic coefficients.
* 360-sweep covering both of the above funtionalities. The output will be a look-up table that can be accessed with an interpolator.

## Structure of the code:

The code may be divided into some well-differentiated sections:
1. The environment
2. Interaction with the environment
3. Area calculation
4. Numerical models
5. Post-processing of aerodynamic coefficients
6. Meshing capabilities

Firstly, the environmental conditions where the satellite intends to operate should be addressed. For this, the  [`SatelliteToolbox`](@ref) package is used, utilizing the MSISE00 model [`SatelliteToolbox.NRLMSISE00_Output`](@ref). 

Secondly, the interactions taking place between the satellite and the gas particles of the environment; gas-surface interactions
are calculated with the implementation of the DRIA model (diffuse reflection with incomplete accommodation) [`SatelliteGeometryCalculations.compute_coefficients`](@ref) [`SatelliteGeometryCalculations.DRIA_GSI`](@ref). 

The third and fourth sections are inherently related because while the calculation of areas can be done in many different ways and without prior inputs, the selection of a specific numerical model sets preferences over the way the areas are calculated. For the most general case, the **non-convex case** [`SatelliteGeometryCalculations.areas_nonconvex`](@ref) case, a ray-tracing method is implemented based on a Möller-Trumbore algorithm ([`SatelliteGeometryCalculations.MTalgorithm`](@ref)). It allows the identification of the intercepted facets and the counting of the total amount of rays that strike; each ray has an equivalent area, which after projection based on the viewpoint yields the **projected cross-sectional area** and **total reference area**. For the **convex case** [`SatelliteGeometryCalculations.areas_convex`](@ref), the forward facing areas are identified and projected onto the desired viewpoint [`SatelliteGeometryCalculations.Viewpoint`](@ref).

The fifth section allows the calculation of the center of pressure and aerodynamic torques from the aerodynamic coefficient calculation. [`SatelliteGeometryCalculations.getCoP`](@ref) and [`SatelliteGeometryCalculations.getTorques`](@ref) implement these calculations.

Lastly, and additional functionality has been implemented to further increase the input capabilities and user's ease of use. Given that that a triangular mesh is necessary for the ray-tracer panel method (RTPM) to perform calculations, the user may directly provide the mesh in \*.obj or \*.stl ASCII format (with specific material properties assigned to each triangular facet in \*.json format) or a CAD model in \*.stp format. In the latter case, GMSH API is called to mesh the geometry, giving the output mesh in \*.stl ASCII format.

## Functional flow

In order to provide a first-glance understanding of how the code works, the following functional flowchart is provided:

![alt](C:\Users\danie\Documents\UC3M\IENAI\Work\S10\flowchart\architecture.png)