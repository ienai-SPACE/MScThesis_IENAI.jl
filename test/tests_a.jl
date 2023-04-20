using Test

@testset "AreasCovex" begin
    vertices3 = [0 0 0 0 1 1 0 -1 1] #x1y1z1x2y2z2x3y3z3
    dir = [1 1 0]
    vertices4 = [0 0 0 0 1 0 0 1 1 0 0 1]

    @test areasConvex(vertices3, dir) == [1.0, 0.7853981633974484]
    @test areasConvex(vertices4, dir) == [1.0, 0.7853981633974484]
end

@testset "area function with convex shape" begin
    triangles = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
    dir = @SVector [1, 1, 0]
    convexFlag = 1

    Aproj, Aref, OutFacets = areas(rmax, distance, dir, triangles, convexFlag)

    Aproj == 1.4142135623730947
    Aref == 1.7320508075688772
    @test OutFacets == [1.0 2.0; 0.8660254037844386 0.8660254037844386; 0.6154797086703875 0.6154797086703875]
end

@testset "Oxygen partial pressure other environmental calculations" begin
    JD = date_to_jd(2018, 6, 19, 18, 35, 0)          #Julian Day [UTC].
    alt = 300e3                                       #Altitude [m].
    g_lat = deg2rad(-22)      #Geodetic latitude [rad].
    g_long = deg2rad(-45)     #Geodetic longitude [rad].
    f107A = 73.5       #81 day average of F10.7 flux (centered on day of year doy).
    f107 = 79      #Daily F10.7 flux for previous day.
    ap = 5.13        #Magnetic index.

    nrlmsise00_output = nrlmsise00(JD, alt, g_lat, g_long, f107A, f107, ap, output_si=true, dversion=true)
    @test OxyPartPress(nrlmsise00_output) == 1.8056322369751164e-5

    @test fEnvironmentlCalcs(JD, alt, g_lat, g_long, f107A, f107, ap) == 
end

