# Inputs/Outputs discussion

## Outputs
The outputs the software offers are:
* The projected cross-sectional area (Aproj)
* The total area of impinged surfaces (Atot): the addition of all mesh elements that have been intercepted by a gas particle
* The aerodynamic coefficients in vectorial form (CD, CL, CP, Cτ). With the prime goal of quick and reliable calculations, two main functionalities are created, each providing the aforementioned outputs:
    * **F1 - Single viewpoint**: the calculations are performed for a specified direction, defined by a local spherical reference frame (α, ϕ), or by the direction of the velocity vector.
    * **F2 - Sweep viewpoint**: a look-up table is created for a sweep of directions such that: α = [−180º : step : 180º] and ϕ = [−90º : step : 90º]. This operation, even though computationally expensive for small values of step, is performed only once. In order to obtain the outputs for any direction, values are interpolated from the look-up table.

## Inputs

For both functionalities, information about the environmental and dynamic conditions, and geometry needs to be set.

* **Orbital information**: orbital elements are used to determine the relative velocity and altitude.
* **Julian date, geodetic coordinates, and solar/magnetic indices**: as inputs of the atmospheric model.
* **Triangular mesh of the 3D model to be evaluated in ’.obj’ format**. Note that depending on the wind-up scheme of the triangles, the definition of the normal to the surface of that element may vary; in order to adjust to a unique convention, triangles need to be numbered in a counter-clockwise fashion. To ensure this condition is satisfied, one may use free software such as Blender to set all normals pointing outward from the body. Triangular elements have been selected as the preferred mesh element due to their mesh adaptability to complex shapes and their geometrical properties to calculate features -such as areas and normals- in 3D space.
* **Information on the wall temperature and composition**: they are needed for the GSI calculations. Regarding the composition of the satellite surface, facet-specific material properties can be assigned by loading a ’.json’ file where this information is stored. See examples in the *How-to guide*.