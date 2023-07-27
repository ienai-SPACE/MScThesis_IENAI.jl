
using SatelliteGeometryCalculations, DelimitedFiles

SatelliteGeometryCalculations.tick()

using FilePathsBase
using FilePathsBase: /

pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "GRACE.obj")
geo = load_geometry(mesh_path, SurfaceProps(), false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

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

    writedlm("GRACE_AprojLookUpTable.txt", AprojLookUpTable)
    writedlm("GRACE_CDlookup.txt", CdLookUp)
    writedlm("GRACE_CLlookup.txt", ClLookUp)
    writedlm("GRACE_CPlookup.txt", CpLookUp)
    writedlm("GRACE_CTAUlookup.txt", CtauLookUp)
    writedlm("GRACE_cullingRatio", culling)
    # writedlm("TSAT_AprojLookUpTable.txt", AprojLookUpTable)

end

##---------- Interpolation --------------------------------------------------------------------
interpFlag = 0
if interpFlag == 1
    _A_look_up_table = readdlm("GRACE_AprojLookUpTable.txt")
    _Cd_look_up_table = readdlm("GRACE_CDlookup.txt")
    alpha_in = 142
    phi_in = -20

    Aproj_int = SatelliteGeometryCalculations.interpolator(_A_look_up_table, alpha_in, phi_in)
    CD_int = SatelliteGeometryCalculations.interpolator(_Cd_look_up_table, alpha_in, phi_in)

end


SatelliteGeometryCalculations.tock()