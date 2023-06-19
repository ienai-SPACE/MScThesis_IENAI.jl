
"""
    areasSpherical(rmax, distance, α, ϕ, MeshVerticesCoords, convexFlag)

#INPUT:
- `rmax`            : radius of the circular plane from where rays originate
- `distance`        : distance at which the circular plane is located (it should be out from the satellite body)
- `α`               : azimuth w.r.t. body frame. (α = 0 aligned with the x-axis)
- `ϕ`               : elevation w.r.t. body frame. (ϕ = 0 on the xy-plane) 
- `MeshVerticesCoords`       : vertices coordinates of the triangular/quad mesh element
- `convexFlag`      :1 for convex, 0 for non-convex
#OUTPUT:
- `OutLMNTs:: OutGeometry{T}`                                   : struct of 3 fields (`area`, `angle`, and `normals`) each containing a vector with the respective magnitude of all intercepted surfaces
- `int_geos::Vector{InteractionGeometry{T}}`                    : vector of struct storing the areas and angles of all intercepted surfaces       
- `Aproj`                                                       : projection of the intercepted triangular areas onto the selected direction 
- `Atot`                                                        : sum of all intercepted triangular areas
"""

function areasSpherical(rmax, distance, α, ϕ, MeshVerticesCoords, convexFlag)

    #from spherical to cartesian coordinates
    ux = cos(ϕ) * cos(α)
    uy = cos(ϕ) * sin(α)
    uz = sin(ϕ)
    dir = SV3([-ux, -uy, -uz]) #unitary directional vector (same sense as oncoming ray beam)

    Aproj, Atot, OutLMNTs, int_geos = areas(rmax, distance, dir, MeshVerticesCoords, convexFlag)      #calculation of areas and normals to the impinged surfaces


    return Aproj, Atot, OutLMNTs, int_geos

end

"""
    areasSpherical(rmax, distance, Vdir, MeshVerticesCoords, convexFlag)

#INPUT:
- `rmax`            : radius of the circular plane from where rays originate
- `distance`        : distance at which the circular plane is located (it should be out from the satellite body)
- `Vrel`            : relative velocity vector
- `MeshVerticesCoords`       : vertices coordinates of the triangular/quad mesh element
- `convexFlag`      :1 for convex, 0 for non-convex
#OUTPUT:
- `OutLMNTs:: OutGeometry{T}`                                   : struct of 3 fields (`area`, `angle`, and `normals`) each containing a vector with the respective magnitude of all intercepted surfaces
- `int_geos::Vector{InteractionGeometry{T}}`                    : vector of struct storing the areas and angles of all intercepted surfaces       
- `Aproj`                                                       : projection of the intercepted triangular areas onto the selected direction 
- `Atot`                                                        : sum of all intercepted triangular areas
"""

function areasSpherical(rmax, distance, Vrel, MeshVerticesCoords, convexFlag)

    Vdir = normalize(Vrel)
    dir = SV3(-Vdir)
    Aproj, Atot, OutLMNTs, int_geos = areas(rmax, distance, dir, MeshVerticesCoords, convexFlag)      #calculation of areas and normals to the impinged surfaces

    return Aproj, Atot, OutLMNTs, int_geos

end