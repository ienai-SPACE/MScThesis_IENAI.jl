
function areasSpherical(rmax, distance, α, ϕ, MeshVerticesCoords, convexFlag)

    #from spherical to cartesian coordinates
    ux = cos(ϕ) * cos(α)
    uy = cos(ϕ) * sin(α)
    uz = sin(ϕ)
    dir = SV3([ux, uy, uz]) #unitary directional vector (opposite sense to oncoming ray beam)

    Aproj, Atot, OutLMNTs, int_geos = SatelliteGeometryCalculations.areas(rmax, distance, dir, MeshVerticesCoords, convexFlag)      #calculation of areas and normals to the impinged surfaces


    return Aproj, Atot, OutLMNTs, int_geos

end