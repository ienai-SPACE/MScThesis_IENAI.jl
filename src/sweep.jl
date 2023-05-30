struct SweepStorage{Float64}
    azimuth::Float64
    altitude::Float64
    Aproj::Float64
    Atot::Float64
    OutLMNTs::OutGeometry{Float64}
end
#_eltype(::SweepStorage{T}) where {T} = T


function sweep(rmax, distance, MeshVerticesCoords, convexFlag, step)

    #step = pi / 20
    α = 0:step:pi/2
    ϕ = 0:step:pi

    #pre-allocation of the vector to be populated by structs
    #T = _eltype(SweepStorage{T})
    AlphaPhiStorage = Matrix(undef, length(ϕ), length(α))
    AprojLookUp = Matrix(undef, length(ϕ), length(α))

    counter = 0
    counter2 = 0

    for aa ∈ α
        counter2 += 1
        for pp ∈ ϕ

            # counter += 1
            #AlphaPhiStorage[counter,counter2] = fIlluminationConvex(triangles, aa, pp) |> SweepStorage

            Aproj, Aref, OutLMNTs, areas_and_angles = areasSpherical(rmax, distance, aa, pp, MeshVerticesCoords, convexFlag)
            counter += 1
            storageValues = SweepStorage(rad2deg(aa), rad2deg(pp), Aproj, Aref, OutLMNTs)
            AlphaPhiStorage[counter, counter2] = storageValues
            AprojLookUp[counter, counter2] = Aproj

        end
        counter = 0
    end

    return AlphaPhiStorage, AprojLookUp

end



