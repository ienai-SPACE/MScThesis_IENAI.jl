using Test, SatelliteToolbox, StaticArrays, SatelliteGeometryCalculations, SpecialFunctions
using FilePathsBase
using FilePathsBase: /

@testset "analyzeAreas and compute_coefficients for CONVEX" begin

    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphereMesh.obj")
    geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

    #---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
    outSurfaceProps = SurfaceProps(0.0, 300.0, 0.0, 0.0, 26.9815)                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

    #----Orbit and date inputs------------------------------------------------------------------
    JD = date_to_jd(2018, 6, 19, 18, 35, 0)     #Julian Day [UTC].
    h = 300
    alt = h * 1000                                #Altitude [m].
    g_lat = deg2rad(-22)                        #Geodetic latitude [rad].
    g_long = deg2rad(-45)                       #Geodetic longitude [rad].
    f107A = 73.5                                #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79                                   #Daily F10.7 flux for previous day.
    ap = 5.13                                   #Magnetic index.
    Vrel = 7629.8875635332
    Vrel_v = [1.0, 0.0, 0.0] * Vrel
    # JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
    outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

    α = deg2rad(0)
    ϕ = deg2rad(0)
    v = Viewpoint(geo, α, ϕ)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Aref, intercept_info, normals = analyze_areas(geo, v)

    @test Aproj ≈ 0.7773842779770798
    @test Aref ≈ 1.552052015831979
    @test length(intercept_info) == 308
    @test length(normals) == 308

    coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
    @test coeffs[1] ≈ 2.2393444265693843
    @test coeffs[3] ≈ 0.0004955217635818692
    @test coeffs[5] ≈ 1.210758832881363
    @test coeffs[7] ≈ 1.0285862895743416
    @test coeffs[2] ≈ [-1.0, 0.0, 0.0]
    @test coeffs[4] ≈ [0.0, -0.4242967712738129, 0.9055231912472577]
    @test coeffs[6] ≈ [-0.9999996646571998, 0.00034397188560096874, 0.0007432151973150171]
    @test coeffs[8] ≈ [-0.9999997181888721, -0.0006092977220194402, -0.00043860969250936286]
end

@testset "analyzeAreas and compute_coefficients for NON-CONVEX" begin

    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "T_SatMesh.obj")
    geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

    #---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
    outSurfaceProps = SurfaceProps(0.0, 300.0, 0.0, 0.0, 26.9815)                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]

    #----Orbit and date inputs------------------------------------------------------------------
    JD = date_to_jd(2018, 6, 19, 18, 35, 0)     #Julian Day [UTC].
    h = 300
    alt = h * 1000                                #Altitude [m].
    g_lat = deg2rad(-22)                        #Geodetic latitude [rad].
    g_long = deg2rad(-45)                       #Geodetic longitude [rad].
    f107A = 73.5                                #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79                                   #Daily F10.7 flux for previous day.
    ap = 5.13                                   #Magnetic index.
    Vrel = 7629.8875635332
    Vrel_v = [1.0, 0.0, 0.0] * Vrel
    # JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
    outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

    α = deg2rad(0)
    ϕ = deg2rad(0)
    v = Viewpoint(geo, α, ϕ)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Aref, intercept_info, normals = analyze_areas(geo, v)

    @test Aproj ≈ 0.2099296370954245
    @test Aref ≈ 0.2099296371694464
    @test length(intercept_info) == 80
    @test length(normals) == 80

    coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
    @test coeffs[1] ≈ 2.287456085225732
    @test coeffs[3] ≈ 1.024395736336198e-6
    @test coeffs[5] ≈ 2.287456083828845
    @test coeffs[7] ≈ 7.129319419221911e-6
    @test coeffs[2] ≈ [-1.0, 0.0, 0.0]
    @test coeffs[4] ≈ [0.0, 0.9940389486206442, 0.10902554115969573]
    @test coeffs[6] ≈ [-0.9999999999940649, 3.436514327611419e-6, 2.461974684398811e-7]
    @test coeffs[8] ≈ [-0.00019754849590484602, -0.9940389292201997, -0.10902553906929591]
end




@testset "DRIA_GSI" begin


    JD = date_to_jd(2018, 6, 19, 18, 35, 0)     #Julian Day [UTC].
    h = 300
    alt = h * 1000                                #Altitude [m].
    g_lat = deg2rad(-22)                        #Geodetic latitude [rad].
    g_long = deg2rad(-45)                       #Geodetic longitude [rad].
    f107A = 73.5                                #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79                                   #Daily F10.7 flux for previous day.
    ap = 5.13                                   #Magnetic index.
    Vrel = 7629.8875635332
    Vrel_v = [1.0, 0.0, 0.0] * Vrel
    # JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
    outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]




    δ = deg2rad(20)
    C = outGasStreamProps.C
    m_srf = 26.9815
    P0 = outGasStreamProps.PO
    Tw = 300.0
    Ta = outGasStreamProps.Ta
    Vrel_norm = 7629.8875635332



    K = 1.44e6                                  #Best-fit Langmuir adsorbate constants for DRIA GSI model
    θ = K * P0 / (1 + K * P0)                   #Fraction of the surface contaminated by atomic oxygen


    γ = cos(δ)
    ell = sin(δ)



    # #-------Pre-allocation---------------------------------------------
    mgas = @SVector [O.m, H.m, N.m, O2.m, N2.m, He.m]                   #[g/mol] atomic mass of gas constituents
    #mgas = @SVector [m[1].m, m[2].m, m[3].m, m[4].m, m[5].m, m[6].m] #[g/mol] atomic mass of gas constituents
    cd_j = @MMatrix zeros(length(mgas), 2)                               #drag coefficient at each facet
    cl_j = @MMatrix zeros(length(mgas), 2)                               #lift coefficient at each facet
    #s = @MArray [zeros(length(mgas), 2)] 
    s = @MMatrix zeros(length(mgas), 2)                                  #thermal speed: two values per species, i.e. contaminated and clean surface
    # #------------------------------------------------------------------


    for jj ∈ 1:2     #clean and contaminated surface
        for j ∈ 1:6  #species specific mass concentration

            μ_srf = mgas[j] / m_srf                 #ratio between the mass of the atoms of the incoming gas with the mass of the surface particles
            Ks = 3.6                                #substrate coefficient (6 < s < 11 the use of Ks = 3.6 is appropriate)
            α_c = Ks * μ_srf / (1 + μ_srf)^2        #accomodation coefficient for clean surface

            α_vec = [α_c, 1]

            α = α_vec[jj]



            s[j, jj] = Vrel_norm / sqrt(2 * (kb / (mgas[j] / NA / 1000) * Ta))    #thermal speed
            P = exp.(-γ .^ 2 .* s[j, jj] .^ 2) ./ s[j, jj]
            G = 1 / (2 * s[j, jj] .^ 2)
            Q = 1 + G
            Z = 1 + erf.(γ .* s[j, jj])
            RR = R / mgas[j] * 1000
            Vratio = sqrt((1 / 2) * (1 + α * ((4 * RR * Tw) / Vrel_norm^2 - 1)))  #[Doornbos 2012]
            cd_j[j, jj] = (P ./ sqrt(π) .+ γ .* Q .* Z .+ (0.5 * γ) .* Vratio .* (γ .* sqrt(pi) .* Z .+ P)) .* C[j] * mgas[j]
            cl_j[j, jj] = (ell * G .* Z .+ (0.5 .* ell) * Vratio .* (γ .* sqrt(pi) .* Z .+ P)) .* C[j] * mgas[j]
        end
    end



    Cd_weighted = sum(cd_j, dims=1) / sum(C .* mgas)         # mass weighted
    Cl_weighted = sum(cl_j, dims=1) / sum(C .* mgas)         # mass weighted

    Cd_weighted_clean = Cd_weighted[1, 1]                    # clean (α ~= 1)
    Cl_weighted_clean = Cl_weighted[1, 1]                    # clean (α ~= 1)
    Cd_weighted_cont = Cd_weighted[1, 2]                     # contaminated (α = 1)
    Cl_weighted_cont = Cl_weighted[1, 2]                     # contaminated (α = 1)


    Cd_facet = (Cd_weighted_clean * (1 - θ) + Cd_weighted_cont * θ)
    Cl_facet = (Cl_weighted_clean * (1 - θ) + Cl_weighted_cont * θ)

    Cp_facet = Cd_facet * cos(δ) + Cl_facet * sin(δ)
    Ctau_facet = Cd_facet * sin(δ) - Cl_facet * cos(δ)

    @test Cd_facet ≈ 2.1366229680840956
    @test Cl_facet ≈ 0.09358046743714629
    @test Cp_facet ≈ 2.0397752413956463
    @test Ctau_facet ≈ 0.6428312190766557
end

@testset "Euclidean norm" begin

    MeshVerticesCoords = [0 0 0 1 1 1 5 5 5; 0 0 0 0.1 0.1 0.1 0.5 0.5 0.5]
    max = SatelliteGeometryCalculations.euclidean_norm(MeshVerticesCoords::Matrix{Float64})

    @test max ≈ 8.660254037844387

end

@testset "Check forward facing" begin

    Vrel_v = [1.0, 0.0, 0.0] #only direction is necessary for this test
    viewpoint1 = SatelliteGeometryCalculations.Viewpoint(0.0, 0.0, Vrel_v)
    Vrel_v2 = [-1.0, 0.0, 0.0] #only direction is necessary for this test
    viewpoint2 = SatelliteGeometryCalculations.Viewpoint(0.0, 0.0, Vrel_v2)
    Vrel_v3 = [0.0, 1.0, 0.0] #only direction is necessary for this test
    viewpoint3 = SatelliteGeometryCalculations.Viewpoint(0.0, 0.0, Vrel_v3)
    v1 = SVector(0.0, 0.0, 0.0)
    v2 = SVector(0.0, 1.0, 0.0)
    v3 = SVector(0.0, 0.0, 1.0)
    face = SatelliteGeometryCalculations.TriangleFace(v1, v2, v3)

    @test SatelliteGeometryCalculations.is_visible(face, viewpoint1) == true
    @test SatelliteGeometryCalculations.is_visible(face, viewpoint2) == false
    @test SatelliteGeometryCalculations.is_visible(face, viewpoint3) == false  #90º wrt normal is not included to avoid numerical errors

end

@testset "MT algorithm" begin
    v1 = SVector(0.0, 0.0, 0.0)
    v2 = SVector(0.0, 1.0, 0.0)
    v3 = SVector(0.0, 0.0, 1.0)
    triangle = SatelliteGeometryCalculations.TriangleFace(v1, v2, v3)

    origin = SVector(10.0, 0.0, 0.0)
    direction = SVector(-1.0, 0.0, 0.0)
    ray = SatelliteGeometryCalculations.Ray(origin, direction)

    MT_result = SatelliteGeometryCalculations.MTalgorithm(triangle, ray; ϵ=sqrt(eps(Float64)))

    # @test (MT_result.mode == FrontFaceIntersection) == true       #doubts here
    @test MT_result.t ≈ 10.0
    @test MT_result.γ_dir ≈ 0.0
    @test MT_result.area ≈ 0.5
    @test MT_result.face_index ≈ 0

end

@testset "projection" begin

    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphereMesh.obj")
    geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

    α = deg2rad(0)
    ϕ = deg2rad(0)
    viewpoint = Viewpoint(geo, α, ϕ)

    #back-face culling
    filtered_geometry = SatelliteGeometryCalculations.filter_backfaces(geo, viewpoint)
    #Number of triangles
    Ntri = SatelliteGeometryCalculations.n_faces(filtered_geometry)
    #Calculate areas

    Aproj = sum(SatelliteGeometryCalculations.projection(filtered_geometry, viewpoint, Ntri))
    @test Aproj ≈ 0.7773512745188897

end

@testset "sampler origins" begin
    using SatelliteGeometryCalculations, DelimitedFiles

    using FilePathsBase
    using FilePathsBase: /

    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh.obj")
    geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

    α = deg2rad(90)
    ϕ = deg2rad(90)
    v = Viewpoint(geo, α, ϕ)

    samplerF = FibonacciSampler(1e5)
    sampler = samplerF

    O, Norig, Aray = SatelliteGeometryCalculations.generate_ray_origins(sampler, v.direction, v.rmax, v.distance)

    @test Norig == 144515
    @test Aray ≈ 9.999932351268763e-6
    @test Aray * Norig / (pi * v.rmax^2) < 1.00001
    @test Aray * Norig / (pi * v.rmax^2) > 0.99999
end

@testset "raytrace" begin
    using SatelliteGeometryCalculations, DelimitedFiles

    using FilePathsBase
    using FilePathsBase: /

    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "boxMesh.obj")
    geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

    α = deg2rad(90)
    ϕ = deg2rad(90)
    v = Viewpoint(geo, α, ϕ)

    samplerF = FibonacciSampler(1e5)
    sampler = samplerF

    rti_vec, filtered_geometry, Aray = SatelliteGeometryCalculations.raytrace(geo, v, sampler)

    @test rti_vec[1].area ≈ 0.12499849701503256
    @test Aray ≈ 9.999976615704741e-6
end

@testset "barycenter_triangle" begin
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "singleTriangle.obj")
    #HOMOGENEOUS CASE
    geo = load_geometry(mesh_path, true, "mm") #homogeneous, convex case


    b = SatelliteGeometryCalculations.face_barycenters(geo, 1)
    b_theo = [0, 0.2 / 3, 0]
    eps = 1e-5
    @test abs(b[1] - b_theo[1]) < eps
    @test abs(b[2] - b_theo[2]) < eps
    @test abs(b[3] - b_theo[3]) < eps

end

@testset "CoP_Circumference_convex" begin
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "flatCircumference.obj")
    #HOMOGENEOUS CASE
    geo = load_geometry(mesh_path, true, "mm") #homogeneous, convex case

    α = deg2rad(0)  #rotate around z-axis
    ϕ = deg2rad(90) #rotate around x-axis
    v = Viewpoint(geo, α, ϕ)
    # v = Viewpoint(geo, Vrel_v)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Atot, intercept_info, normals, culling, barycenters = analyze_areas(geo, v)

    CoP = SatelliteGeometryCalculations.getCoP(Aproj, intercept_info, barycenters)

    CoP_theo = [0, 0, 0]
    eps = 1e-5
    @test abs(CoP[1] - CoP_theo[1]) < eps
    @test abs(CoP[2] - CoP_theo[2]) < eps
    @test abs(CoP[3] - CoP_theo[3]) < eps
end

@testset "CoP_Circumference_nonconvex" begin
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "flatCircumference.obj")
    #HOMOGENEOUS CASE
    geo = load_geometry(mesh_path, false, "mm") #homogeneous, convex case

    α = deg2rad(0)  #rotate around z-axis
    ϕ = deg2rad(90) #rotate around x-axis
    v = Viewpoint(geo, α, ϕ)
    # v = Viewpoint(geo, Vrel_v)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Atot, intercept_info, normals, culling, barycenters = analyze_areas(geo, v)

    CoP = SatelliteGeometryCalculations.getCoP(Aproj, intercept_info, barycenters)

    CoP_theo = [0, 0, 0]
    eps = 1e-5
    @test abs(CoP[1] - CoP_theo[1]) < eps
    @test abs(CoP[2] - CoP_theo[2]) < eps
    @test abs(CoP[3] - CoP_theo[3]) < eps
end

@testset "CoP_Sphere_convex" begin
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphereMesh.obj")
    #HOMOGENEOUS CASE
    geo = load_geometry(mesh_path, true, "mm") #homogeneous, convex case

    α = deg2rad(0)  #rotate around z-axis
    ϕ = deg2rad(90) #rotate around x-axis
    v = Viewpoint(geo, α, ϕ)
    # v = Viewpoint(geo, Vrel_v)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Atot, intercept_info, normals, culling, barycenters = analyze_areas(geo, v)

    CoP = SatelliteGeometryCalculations.getCoP(Aproj, intercept_info, barycenters)

    CoP_theo = [0, 0, 0.3298]
    eps = 1e-3
    @test abs(CoP[1] - CoP_theo[1]) < eps
    @test abs(CoP[2] - CoP_theo[2]) < eps
    @test abs(CoP[3] - CoP_theo[3]) < eps
end

@testset "CoP_Sphere_nonconvex" begin
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "sphereMesh.obj")
    #HOMOGENEOUS CASE
    geo = load_geometry(mesh_path, false, "mm") #homogeneous, convex case

    α = deg2rad(0)  #rotate around z-axis
    ϕ = deg2rad(90) #rotate around x-axis
    v = Viewpoint(geo, α, ϕ)
    # v = Viewpoint(geo, Vrel_v)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Atot, intercept_info, normals, culling, barycenters = analyze_areas(geo, v)

    CoP = SatelliteGeometryCalculations.getCoP(Aproj, intercept_info, barycenters)

    CoP_theo = [0, 0, 0.3298]
    eps = 1e-3
    @test abs(CoP[1] - CoP_theo[1]) < eps
    @test abs(CoP[2] - CoP_theo[2]) < eps
    @test abs(CoP[3] - CoP_theo[3]) < eps
end

@testset "Torque and CoP" begin
    #The sum of moments about the CoP shall be zero

    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "T_SatMesh.obj")
    materials_path = FilePathsBase.join(pkg_path, "test", "inputs_models_data", "TSAT_coarse_mesh_materials.json")

    #HOMOGENEOUS CASE
    geo = load_geometry(mesh_path, false, "mm") # UNITS: "m" -> meters and "mm" -> milimiters

    #---------- # EVALUATION OF A SINGLE VIEWPOINT DIRECTION # --------------------------------------
    outSurfaceProps = SurfaceProps()                                                       #outSurfaceProps.[η, Tw, s_cr, s_cd, m_srf]
    #----Orbit and date inputs------------------------------------------------------------------
    JD, alt, g_lat, g_long, f107A, f107, ap, Vrel_v = SatelliteGeometryCalculations.OrbitandDate()
    outGasStreamProps = GasStreamProperties(JD, alt, g_lat, g_long, f107A, f107, ap)       #outGasStreamProps.[C, PO, mmean, Ta]

    α = deg2rad(180)  #rotate around z-axis
    ϕ = deg2rad(0) #rotate around x-axis
    v = Viewpoint(geo, α, ϕ)
    # v = Viewpoint(geo, Vrel_v)

    #---- Area calculations --------------------------------------------------------------------
    Aproj, Atot, intercept_info, normals, culling, barycenters = analyze_areas(geo, v)
    #---- Aerodynamic coefficients ---------------------------------------------------------------
    coeffs, Atot, Aproj, coeffs_v = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
    #----- Torques and CoP ---------------------------------------------------------------------

    torque_ref = SVector(100.0, 0.0, 0.0)  # point about which moments are calculated
    CT = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, barycenters, torque_ref)

    chord_plane_n = SVector(0.0, 1.0, 0.0) #chord plane normal
    CoP = SatelliteGeometryCalculations.getCoP(T, coeffs, chord_plane_n)

    CT2 = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, barycenters, CoP[2])
    CT3 = SatelliteGeometryCalculations.getTorques(coeffs_v, intercept_info, barycenters, CoP[1])

    check = CT2 #center of pressure on the chord plane
    check2 = CT3 #center of pressure on a point along the line of application
    eps = 5e-3

    @test abs(check[1]) < eps
    @test abs(check[2]) < eps
    @test abs(check[3]) < eps
    @test abs(check2[1]) < eps
    @test abs(check2[2]) < eps
    @test abs(check2[3]) < eps
end