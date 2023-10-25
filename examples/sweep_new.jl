
using SatelliteGeometryCalculations, DelimitedFiles

# SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")
#HETEROGENEOUS CASE
# geo = load_geometry(mesh_path, materials_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters
#HOMOGENEOUS CASE
geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

#---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

#----Orbit and date inputs------------------------------------------------------------------
JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]


#---------- # SWEEP TO GENERATE LOOK-UP TABLE # ----------------------------------------------
sweepFlag = 1
if sweepFlag == 1
    step = deg2rad(15)
    grid = SatelliteGeometryCalculations.Grid(step)
    # LookUpTable, AprojLookUpTable = SatelliteGeometryCalculations.sweep_v2(geo, grid)
    LookUpTable, AprojLookUpTable, CdLookUp, ClLookUp, CpLookUp, CtauLookUp, culling = SatelliteGeometryCalculations.sweep_v2(geo, grid, outSurfaceProps, outGasStreamProps, Vrel_v)

    writedlm(string(FilePathsBase.join(pkg_path, "test", "inputs_models_data", "results", "homo_TSAT_coarse_AprojLookUpTable.txt")), AprojLookUpTable)
    writedlm(string(FilePathsBase.join(pkg_path, "test", "inputs_models_data", "results", "homo_TSAT_coarse_CDlookup.txt")), CdLookUp)
end

##---------- Interpolation --------------------------------------------------------------------
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
full_path = string(FilePathsBase.join(pkg_path, "test", "inputs_models_data", "results", "homo_TSAT_coarse_CDlookup.txt"))

α = 90  #[deg]
ϕ = 0   #[deg]
out_interpolant = SatelliteGeometryCalculations.interpolant_RTPM(full_path, α, ϕ)


# SatelliteGeometryCalculations.tock()
