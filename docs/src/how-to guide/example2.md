# Example 2: calculation of areas

Two geometric scenarios are possible: **convex** and **non-convex** geometries, as indicated in [Example 1: loading the geometry and material inputs](@ref). This example is applicable to both cases, except for the [Selection of the sampler/filter](@ref) section.


## Selection of the sampler/filter

***Only applicable to non-convex geometries!***

The generation of the circular source from which the rays are projected onto the satellite can be modified, having 3 options: **Grid Filter**, **Fibonacci Sampler**, and **Monte Carlo Sampler**. Moreover, the density of the sampler [rays/m^3] may also be modified `samplerX = nameSampler(density)`.

![alt](C:\Users\danie\Documents\UC3M\IENAI internship\Overleaf\figures\FibonacciSampling1000.jpg)

These options may be modified in the function [`SatelliteGeometryCalculations.areas_nonconvex`](@ref)

```
function areas_nonconvex(geometry::AbstractGeometry, viewpoint::Viewpoint)

    samplerG = GridFilter(1e5)
    samplerF = FibonacciSampler(1e5)
    samplerMC = MonteCarloSampler(1e5)
    sampler = samplerF

    rti_vec, w, filtered_geometry, Aray= raytrace(geometry, viewpoint, sampler) #culling +  ray tracing
    
    ...
```

## Selecting the viewpoint

Two options are available:

### a) Definition of viewpoint based on spherical coordinates in body reference frame
```
α = deg2rad(0)  #rotate around z-axis
ϕ = deg2rad(90) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
```
### b) Viewpoint coincident with the velocity vector (in cartesian form) in body reference frame

```
v = Viewpoint(geo, Vrel_v)
```
Where `Vrel_v` is defined inside of [`SatelliteGeometryCalculations.OrbitandDate()`](@ref).

## Example: running the complete code

```
using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent

mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")

#HOMOGENEOUS CASE
geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters


#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # ---------------------------------
outSurfaceProps = SurfaceProps()             #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

# --- Viewpoint ----------------------------------------------------------------------------
α = deg2rad(0)  #rotate around z-axis
ϕ = deg2rad(90) #rotate around x-axis
v = Viewpoint(geo, α, ϕ)
# v = Viewpoint(geo, Vrel_v)

#---- Area calculations --------------------------------------------------------------------
Aproj, Aref, intercept_info, normals, culling = analyze_areas(geo, v)
```