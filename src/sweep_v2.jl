
"""
    SweepStorage_v2{Float64}

-`azimuth::Float64`
-`altitude::Float64`
-`Aproj::Float64`
-`Atot::Float64`
"""

struct SweepStorage_v2{Float64}
    azimuth::Float64
    altitude::Float64
    Aproj::Float64
    Atot::Float64
end
#_eltype(::SweepStorage{T}) where {T} = T


"""
    sweep_v2(geo, grid, outSurfaceProps, outGasStreamProps, Vrel_v)

Create a 2D matrix as function of azimuth and elevation storing the projected `Aproj` and total area `Aref`, and the aerodynamic coefficients

#INPUTS:
- `geo::AbstractGeometry`
- `grid::Grid{T}`       : with fields `alpha::Vector` and `phi::Vector`
- `outSurfaceProps::SurfaceProps{T}`
- `outGasStreamProps::GasStreamProperties{T}`
- `Vrel_v::Vector`

#OUTPUTS:
- `AlphaPhiStorage:: Matrix(undef, lastindex(ϕ), lastindex(α))`     : it contains `SweepStorage{Float64}` with fields `azimuth`, `altitude`, `Aproj`, `Aref`, `OutLMNTs`
- `AprojLookUp::Matrix(undef, lastindex(ϕ), lastindex(α))`          : matrix storing the projected areas for each α, ϕ pair 
- `CdLookUp::Matrix(undef, lastindex(ϕ), lastindex(α))`
- `ClLookUp::Matrix(undef, lastindex(ϕ), lastindex(α))`
- `CpLookUp::Matrix(undef, lastindex(ϕ), lastindex(α))`
- `CtauLookUp::Matrix(undef, lastindex(ϕ), lastindex(α))`
- `culling::Matrix(undef, lastindex(ϕ), lastindex(α))`
"""


function sweep_v2(geo, grid, outSurfaceProps, outGasStreamProps, Vrel_v)
    α = grid.alpha
    ϕ = grid.phi

    #pre-allocation of the vector to be populated by structs
    AlphaPhiStorage = Matrix(undef, length(ϕ), length(α))
    AprojLookUp = Matrix(undef, length(ϕ), length(α))
    CdLookUp = Matrix(undef, length(ϕ), length(α))
    ClLookUp = Matrix(undef, length(ϕ), length(α))
    CpLookUp = Matrix(undef, length(ϕ), length(α))
    CtauLookUp = Matrix(undef, length(ϕ), length(α))
    culling = Matrix(undef, length(ϕ), length(α))

    counter = 0
    counter2 = 0
    counter3 = 0

    for aa ∈ α
        counter2 += 1
        for pp ∈ ϕ
            v = Viewpoint(geo, aa, pp)
            Aproj, Aref, intercept_info, normals, culling_ratio = analyze_areas(geo, v)
            coeffs, Atot, Aproj = compute_coefficients(outSurfaceProps, outGasStreamProps, intercept_info, Vrel_v, normals)
            counter += 1
            storageValues = SweepStorage_v2(rad2deg(aa), rad2deg(pp), Aproj, Aref)
            AlphaPhiStorage[counter, counter2] = storageValues
            AprojLookUp[counter, counter2] = Aproj
            CdLookUp[counter, counter2] = coeffs[1]
            ClLookUp[counter, counter2] = coeffs[3]
            CpLookUp[counter, counter2] = coeffs[5]
            CtauLookUp[counter, counter2] = coeffs[7]
            culling[counter, counter2] = culling_ratio

            counter3 += 1
            println("iter num=", counter3)
        end
        counter = 0
    end


    # return AlphaPhiStorage, AprojLookUp
    return AlphaPhiStorage, AprojLookUp, CdLookUp, ClLookUp, CpLookUp, CtauLookUp, culling

end



