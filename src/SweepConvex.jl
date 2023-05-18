#---------------------------------------------------------------------
using SatelliteGeometryCalculations, StaticArrays, LinearAlgebra

include("IlluminationConvex.jl")
#---------------------------------------------------------------------

struct SweepStorage{Float64}
    azimuth::Float64
    altitude::Float64
    Aproj::Float64
    Atot::Float64
    OutLMNTs::auxOutGeometry{Float64} #previously: OutGeometry{Float64}
end
#_eltype(::SweepStorage{T}) where {T} = T


function sweep(triangles::Matrix, step)

    #step = pi / 20
    α = 0:step:pi/2
    ϕ = 0:step:pi

    #pre-allocation of the vector to be populated by structs
    #T = _eltype(SweepStorage{T})
    AlphaPhiStorage = Matrix(undef, length(ϕ), length(α))

    counter = 0
    counter2 = 0

    for aa ∈ α
        counter2 += 1
        for pp ∈ ϕ

            # counter += 1
            #AlphaPhiStorage[counter,counter2] = fIlluminationConvex(triangles, aa, pp) |> SweepStorage

            Aproj, Aref, OutLMNTs, areas_and_angles = fIlluminationConvex(triangles, aa, pp)
            counter += 1
            storageValues = SweepStorage(rad2deg(aa), rad2deg(pp), Aproj, Aref, OutLMNTs)
            AlphaPhiStorage[counter, counter2] = storageValues

        end
        counter = 0
    end

    return AlphaPhiStorage

end


# #piping example
# #---------------------------------
# function f(x, y)
#     x * y
#     x / 5
# end
# function g(z, t)
#     z^2
#     t * 5
# end

# result = f(1, 2) |> g

#TEST
#-------------------------------------
# MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]



# mesh = load("C:\\Users\\danie\\Documents\\UC3M\\IENAI internship\\CAD\\sphereMesh.obj")

#load the mesh
using FileIO
using FilePathsBase
using FilePathsBase: /
pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
mesh = load(pkg_path / "test" / "samples" / "sphereMesh.obj")


step = pi / 20
# aa = pi
# pp = pi / 2

MeshVerticesCoords = SatelliteGeometryCalculations.finputMesh(mesh; scale=1e-3)

AlphaPhiStorage = sweep(MeshVerticesCoords, step)


