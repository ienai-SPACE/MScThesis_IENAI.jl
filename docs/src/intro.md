# Summary of the tool functionalities

A tool that allows the calculation of aerodynamic coefficients of satellites in LEO is presented. Two main functionalities may be differentiated: the calculation of the aerodynamic coefficients per se, and the calculation of the projected cross-sectional area (PCSA) for a specific viewpoint. A ray-tracing panel method (RTPM) has been selected as the most optimal numerical
solution to address the problem of the area and aerodynamic coefficients calculation of
a satellite. The RTPM provides a reasonably mild computational cost and complexity,
allowing in-time calculations for different scenarios. These features are the main advantages in opposition to DSMC or TPMC, which are more computationally expensive.

The DRIA model
is selected to calculate the aerodynamic coefficients. With regard to the calculation of
the area and later aerodynamic coefficient computations, a ray-tracing particle method is
implemented.

The PCSA is validated and verified against CROC (a module from ESA’s DRAMA
tool), and the computation of aerodynamic coefficients against DRIA’s closed-form solution for a sphere.

Additionally, the tool offers the possibility of differentiating between convex and nonconvex shapes (including self-shading effects) and between a single (homogeneous) or
several (heterogeneous) surface material definition.
