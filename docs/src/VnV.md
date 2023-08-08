# Validation and Verification campaign

The validation and verification phase has been performed at two different levels: area
level and GSI level. The former addresses the performance of the ray-tracer area calculation through ray-counting for non-convex shapes, and the area projection method for convex; whereas the latter assesses the quality of the GSI calculations.

## Area calculation: convex shapes

For the convex case an area-projection method is implemented, because of this reason, the limiting factor is the fidelity with which the meshed geometry corresponds to the real geometry. This scenario has been assessed by evaluating the degree to which the code identifies the area of the surface of a sphere as well as its projected area.

![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\VnV_convex_table.jpg)
![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\mesh_sensitivity_analysis.jpg)

Table 5.1 and the plots underneath show the asymptotic behavior of the projected area, a key parameter involved in the aerodynamic coefficient calculations. Regarding the calculation of the total area, numerical errors that arise in the calculation of the orientation of the facets are assumed to be responsible for the over or underestimation of the area without any specific trend. As a final comment, it is important to remark on the fact that the more elements the more costly the computations are, and thus a compromise may need to be found.


## Area calculation: non-convex shapes

For non-convex shapes, a commercial tool has been used: CROC. CROC (CROss Section of Complex Bodies) is a module of ESA’s Debris Risk Assessment and Mitigation Analysis software. CROC provides the cross-sectional area of the satellite under a specified aspect angle, which is the surface area under which perturbations interact, such as SRP or drag. The CROC's *randomly tumbling satellite functionality* performs a sweep around the whole satellite. This is equivalent to the look-up table functionality of the in-house RTPM code, and thus the justification to compare them. The figure below shows the comparison between the two, and the large level of almost-identical correspondence between both, as shown in table 5.2.

![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\CROC_VV2.jpg)
![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\VnV_nonconvex_table.jpg)

### Look-up table interpolation

Recalling that an interpolating method is implemented to increase the speed of the calculations, the interpolated areas also need to be verified. Four random directions have been compared to CROC measurements. Results are shown in the following table:

![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\VnV_interpolation_table.jpg)

## Gas-surface interactions

The other leg of this software that needs to undergo a validation and verification process is the calculation of the gas-surface interactions. Remember that the code implements a panel method algorithm, finding the aerodynamic coefficients at each facet (control area) and performing a weighted addition of these contributions to obtain the overall aerodynamic coefficients. In essence, when the code performs the calculations on the facets, it is calculating the GSIs over a flat plate, which is one of the closed-form solutions that can be acquired analytically. Another known closed-form solution is that of a sphere. Because of this reason, in order to perform the V&V of the GSI calculations in the in-house code, the aerodynamic coefficients over a meshed sphere have been obtained and compared to the DRIA closed-form solution:

$C_{D, \mathrm{sphere}}=\frac{4 s^4+4 s^2-1}{2 s^4} \operatorname{erf}(s)+\frac{2 s^2+1}{\sqrt{\pi} s^3} e^{-s^2}+\frac{2 \sqrt{\pi}}{3 s} \sqrt{\frac{T_{k, r}}{T_{\infty}}}$


Due to the symmetry of a sphere, only the coefficient of drag is analyzed since the coefficients of lift and shear should be zero.

![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\DRIA_sphere.jpg)

In the analysis, the CD is calculated for different altitudes and two different wall temperatures. Three main observations are drawn from the results shown in the figure above: The first one indicates that the in-house code provides a higher correlation to closed-form results for altitudes up to 350 km (2% difference), although for altitudes from 800-900km the percentage difference between both approaches a 2% again. The second one is that for altitudes from 350-700km the CD shows the maximum values of drag coefficient, reaching its maximum at 500km. Lastly, the increase of the wall temperature causes a higher CD and has no significant difference in the percentage error between the in-house measurements and the closed-form solution. 

Overall, both methods show close results, and the percentage difference between both is much smaller than, or at least equal to deviations from the ’right value’ that could arise from other uncertainties. These uncertainties could be originated by not knowing the exact wall temperature, the partial pressure or concentration of oxygen or other gas species, or even the velocity of rotation of the atmosphere.