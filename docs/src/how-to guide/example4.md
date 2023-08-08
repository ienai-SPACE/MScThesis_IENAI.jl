# Example 4: execution of sweep and interpolation

## Sweep
This example shows how to run a sweep analysis for the input geometry [TSAT_coarse_mesh.obj](@ref) with homogeneous material proporties. Therefore, an **homogeneous** case is run, having selected the **non-convex** method, with input units being [mm]. 

In this analysis the **viewpoint** will vary such that α = [−180º : step : 180º] and ϕ = [−90º : step : 90º], being possible to modify the value of the `step` in `grid = SatelliteGeometryCalculations.Grid(step)` (the value of the step must be in radians).

*The flag `sweepFlag` allows turning ON/OFF this functionality.*

## Interpolation

Once a look-up table is build with the sweep, this one can be accessed.

`alpha_in` and `phi_in` denote the angles, in degrees, that provide the viewpoint to be assessed; these correspond to the rotation angles on the x-y and y-z planes, respectively.  

*The flag `interpFlag` allows turning ON/OFF this functionality.*
```
using SatelliteGeometryCalculations, DelimitedFiles
using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")

#HETEROGENEOUS CASE
# geo = load_geometry(mesh_path, materials_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
outSurfaceProps = SurfaceProps()                                               #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)     #outGasStreamProps.[C, PO, mmean, Ta]

#---------- # SWEEP TO GENERATE LOOK-UP TABLE # ----------------------------------------------
sweepFlag = 1  #Flag to activate the sweep
if sweepFlag == 1
    step = deg2rad(15)
    grid = SatelliteGeometryCalculations.Grid(step)
    LookUpTable, AprojLookUpTable, CdLookUp, ClLookUp, CpLookUp, CtauLookUp, culling = SatelliteGeometryCalculations.sweep_v2(geo, grid, outSurfaceProps, outGasStreamProps, Vrel_v)

    writedlm("homo_TSAT_coarse_AprojLookUpTable.txt", AprojLookUpTable)
    writedlm("homo_TSAT_coarse_CDlookup.txt", CdLookUp)

end

##---------- Interpolation --------------------------------------------------------------------
interpFlag = 0  #Flag to activate the interpolation
if interpFlag == 1
    _A_look_up_table = readdlm("homo_TSAT_coarse_AprojLookUpTable.txt")
    _Cd_look_up_table = readdlm("homo_TSAT_coarse_CDlookup.txt")

    alpha_in = 142
    phi_in = -20

    Aproj_int = SatelliteGeometryCalculations.interpolator(_A_look_up_table, alpha_in, phi_in)
    CD_int = SatelliteGeometryCalculations.interpolator(_Cd_look_up_table, alpha_in, phi_in)

end
```