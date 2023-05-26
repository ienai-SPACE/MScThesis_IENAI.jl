using Test, SatelliteToolbox, StaticArrays, SatelliteGeometryCalculations


@testset "Oxygen partial pressure other environmental calculations" begin
    JD = date_to_jd(2018, 6, 19, 18, 35, 0)           #Julian Day [UTC].
    alt = 300e3                                       #Altitude [m].
    g_lat = deg2rad(-22)      #Geodetic latitude [rad].
    g_long = deg2rad(-45)     #Geodetic longitude [rad].
    f107A = 73.5       #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79      #Daily F10.7 flux for previous day.
    ap = 5.13        #Magnetic index.

    nrlmsise00_output = nrlmsise00(JD, alt, g_lat, g_long, f107A, f107, ap, output_si=true, dversion=true)
    @test oxygen_partial_pressure(nrlmsise00_output) == 1.8056322369751164e-5
    gsp = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)

    @test gsp.C ≈ [0.0014285330240918048, 0.7731601312894321, 0.20321083664911124, 0.009497284958350873, 3.4239449619055476e-5, 0.011599903369159519]
    @test gsp.Ta ≈ 833.3451338126531
    @test gsp.PO ≈ 1.8056322369751164e-5
    @test gsp.mmean ≈ 17.387351734564074

    # outGSP = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)
    # @test outGSP.C == [0.0014285330240918048, 0.7731601312894321, 0.20321083664911124, 0.009497284958350873, 3.4239449619055476e-5, 0.011599903369159519]
    # @test outGSP.Ta == 833.3451338126531
    # @test outGSP.PO == 1.8056322369751164e-5
    # @test outGSP.mmean == 17.387351734564074

end

@testset "AreasConvex" begin
    vertices3 = [0, 0, 0, 0, 1, 1, 0, -1, 1] #x1y1z1x2y2z2x3y3z3
    dir = [-1 -1 0]
    vertices4 = [0 0 0 0 1 0 0 1 1 0 0 1]

    @test areasConvex(vertices3, dir) == [1.0, 0.7853981633974484, 1.0, 0.0, 0.0]
    @test areasConvex(vertices4, dir) == [1.0, 0.7853981633974484, 1.0, 0.0, 0.0]
end

@testset "area function with convex shape" begin
    triangles = [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
    dir = [-1, -1, 0]
    convexFlag = 1
    rmax = 2
    distance = 10

    Aproj, Atot, OutLMNTs, InteractionGeometry_v = SatelliteGeometryCalculations.areas(rmax, distance, dir, triangles, convexFlag)

    @test Aproj ≈ 1.4142135623730947
    @test Atot ≈ 1.7320508075688772
    # @test OutFacets ≈ [1.0 2.0; 0.8660254037844386 0.8660254037844386; 0.6154797086703871 0.6154797086703871]
    @test OutLMNTs.area ≈ [0.8660254037844386; 0.8660254037844386]
    @test OutLMNTs.angle ≈ [0.6154797086703871; 0.6154797086703871]
    @test InteractionGeometry_v[1].area ≈ 0.8660254037844386
    @test InteractionGeometry_v[2].area ≈ 0.8660254037844386
    @test InteractionGeometry_v[1].angle ≈ 0.6154797086703875
    @test InteractionGeometry_v[2].angle ≈ 0.6154797086703875
end


@testset "Coefficients" begin
    JD = date_to_jd(2018, 6, 19, 18, 35, 0)          #Julian Day [UTC].
    alt = 300e3                                       #Altitude [m].
    g_lat = deg2rad(-22)      #Geodetic latitude [rad].
    g_long = deg2rad(-45)     #Geodetic longitude [rad].
    f107A = 73.5       #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79      #Daily F10.7 flux for previous day.
    ap = 5.13        #Magnetic index.


    outSurfaceProps = SurfaceProps()
    outGSP = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)
    #test a single struct as input: intgeo::InteractionGeometry
    interactions_geometries = InteractionGeometry(0.8660254037844386, 0.6154797086703875)
    Vrel_v = [1 / sqrt(2), 1 / sqrt(2), 0] * 7000.0
    normals = [@SVector[0.5773502691896258, 0.5773502691896258, 0.5773502691896258]]

    coeffs = compute_coefficients(outSurfaceProps, outGasStreamProps, interactions_geometries, Vrel_v, normals)
    @test coeffs.Cd ≈ 1.855272320022949
    @test coeffs.Cl ≈ 0.15663734326247042
    @test coeffs.Cp ≈ 1.6052581182863233
    @test coeffs.Ctau ≈ 0.9432481181659013

    #test vector of structs as input: intgeo::Vector{<:InteractionGeometry}
    int_geos = [InteractionGeometry(0.8660254037844386, 0.6154797086703875), InteractionGeometry(0.8660254037844386, 0.6154797086703875)]

    normals2 = [@SVector[0.5773502691896258, 0.5773502691896258, 0.5773502691896258], @SVector[-0.5773502691896258, -0.5773502691896258, 0.5773502691896258]]

    coeffs2, areas, Aproj = compute_coefficients(outSurfaceProps, outGSP, int_geos, Vrel_v, normals2)


    @test coeffs2[1] ≈ 2.2722352589828807
    @test coeffs2[3] ≈ 0.19184078282911463
    @test coeffs2[5] ≈ 1.1350889009950158
    @test coeffs2[7] ≈ 0.6669771406965588
    @test areas ≈ 1.7320508075688772
    @test Aproj ≈ 1.4142135623730947

end


@testset "main" begin
    # using SatelliteGeometryCalculations
    # outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

    # #----Orbit and date inputs-------------------------------------------------------------------------------------------------
    # JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
    # #-----------------------------------------------------------------------------------------------------------------------------

    # outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

    # #----Area calculation inputs-------------------------------------------------------------------------------------------------
    # VdirFlag = 0           # 0: specify direction; 1: direction indicated by velocity vector
    # convexFlag = 0         # set if the satellite is convex (flag == 1) or non-convex (flag == 0)
    # MeshVerticesCoords, dir, rmax, distance = SatelliteGeometryCalculations.GeomInputs(Vrel_v, VdirFlag, convexFlag)     #mesh geometry & direction defined inside
    # # #-----------------------------------------------------------------------------------------------------------------------------

    # Aproj, Atot, OutLMNTs, int_geos = SatelliteGeometryCalculations.areas(rmax, distance, dir, MeshVerticesCoords, convexFlag)      #calculation of areas and normals to the impinged surfaces

    # # # # print(int_geos)

    # # α = deg2rad(45)
    # # ϕ = 0
    # Vrel_norm = 7000.0


    # # #coeffs2, Atot2, Aproj2 = drag_for_orientation_convex(MeshVerticesCoords, outGasStreamProps, outSurfaceProps, α, ϕ, Vrel_norm)

    # # #interactions_geometries = InteractionGeometry(OutLMNTs.area[1], OutLMNTs.angle[1])

    # coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, int_geos, Vrel_norm)
    # # #coeffs = compute_coefficients(outSurfaceProps, outGasStreamProps, interactions_geometries, Vrel_norm)
    # # #(; Cd, Cl, Cp, Ctau) = coeffs
    # # #print(coeffs)

    # print(size(MeshVerticesCoords, 1), //)
    # print(Atot, //)
    # print(Aproj, //)
    # print(coeffs, //)

    using SatelliteGeometryCalculations, StaticArrays, LinearAlgebra
    outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

    #----Orbit and date inputs-------------------------------------------------------------------------------------------------
    JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
    #-----------------------------------------------------------------------------------------------------------------------------

    outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

    #----Area calculation inputs-------------------------------------------------------------------------------------------------
    VdirFlag = 0           # 0: specify direction; 1: direction indicated by velocity vector
    convexFlag = 0         # set if the satellite is convex (flag == 1) or non-convex (flag == 0)
    MeshVerticesCoords, dir, rmax, distance = SatelliteGeometryCalculations.GeomInputs(Vrel_v, VdirFlag, convexFlag)     #mesh geometry & direction defined inside
    # MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
    # dir = @SVector([1, 1, 0])
    # distance = 10
    # rmax = 2
    # #-----------------------------------------------------------------------------------------------------------------------------

    Aproj, Atot, OutLMNTs, int_geos = SatelliteGeometryCalculations.areas(rmax, distance, dir, MeshVerticesCoords, convexFlag)      #calculation of areas and normals to the impinged surfaces


    Vrel_norm = norm(Vrel_v)


    # coeffs2, Atot2, Aproj2 = drag_for_orientation_convex(MeshVerticesCoords, outGasStreamProps, outSurfaceProps, α, ϕ, Vrel_norm)

    interactions_geometries = InteractionGeometry(OutLMNTs.area[1], OutLMNTs.angle[1])

    coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, int_geos, Vrel_v, OutLMNTs.normals)

end

# using LinearAlgebra
# v1 = [1, 1, 0]
# v2 = [0, 1, -1]
# v3 = [1, 0, -1]
# # v2 = [1, 1, 0, 1, 0, -1, 0, 1, -1]
# # v3 = [0.5, 0.5, 0, 1, 1, 1, 0, 0, 1]
# edge1 = v2 - v1
# edge2 = v3 - v2
# crossProd = cross(edge1, edge2) |> LinearAlgebra.normalize