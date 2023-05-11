function drag_for_orientation_convex(geometry::SMatrix, gsp::GasStreamProperties, surfprops::SurfaceProps, α, ϕ, Vrel_norm)

    Aproj, Aref, OutLMNTs, areas_and_angles = fIlluminationConvex(geometry, α, ϕ)

    coeffs, A = compute_coefficients(surfprops, gsp, areas_and_angles, Vrel_norm)

    return coeffs, A

end

