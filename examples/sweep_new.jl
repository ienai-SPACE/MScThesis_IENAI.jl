
using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

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

    writedlm("homo_TSAT_coarse_AprojLookUpTable.txt", AprojLookUpTable)
    writedlm("homo_TSAT_coarse_CDlookup.txt", CdLookUp)

end

##---------- Interpolation --------------------------------------------------------------------
interpFlag = 1
if interpFlag == 1
    _A_look_up_table = readdlm("homo_TSAT_coarse_AprojLookUpTable.txt")
    _Cd_look_up_table = readdlm("homo_TSAT_coarse_CDlookup.txt")
    alpha_in = 0
    phi_in = 90

    Aproj_int = SatelliteGeometryCalculations.interpolator(_A_look_up_table, alpha_in, phi_in)
    CD_int = SatelliteGeometryCalculations.interpolator(_Cd_look_up_table, alpha_in, phi_in)

end


SatelliteGeometryCalculations.tock()
