#---------------------------------------------------------------------
using SatelliteGeometryCalculations, StaticArrays, LinearAlgebra
SV3{T} = SVector{3,T}
include("geometry.jl")
include("IlluminationNonConvex.jl")
include("NonConvex.jl")
include("Origins.jl")
include("MTalgorithm.jl")
#---------------------------------------------------------------------



"""
    SweepStorage{Float64}

-`azimuth::Float64`
-`altitude::Float64`
-`Aproj::Float64`
-`Aref::Float64`
-`OutLMNTs::auxOutGeometry{Float64}`
"""
struct SweepStorage{Float64}
    azimuth::Float64
    altitude::Float64
    Aproj::Float64
    Aref::Float64
    OutLMNTs::auxOutGeometry{Float64}
end


"""
    sweep(triangles::Matrix, step)

Create a 2D matrix as function of azimuth and elevation storing the projected `Aproj` and total area `Aref`, and the specific data of each impinged element `OutLMNTs`

#INPUTS:
- `triangles::SMatrix` : coordinates of the vertices of all area elements
- `step::Float64`      : step used in the definition of `α::StepRangeLen{Float64}` `ϕ::StepRangeLen{Float64}`

#OUTPUTS:
- AlphaPhiStorage:: Matrix(undef, lastindex(ϕ), lastindex(α))    : fields `azimuth`, `altitude`, `Aproj`, `Aref`, `OutLMNTs`
"""


function sweep(triangles::Matrix, step)

    #step = pi / 20
    α = 0:step:pi/2
    ϕ = 0:step:pi

    #pre-allocation of the vector to be populated by structs
    AlphaPhiStorage = Matrix(undef, lastindex(ϕ), lastindex(α))

    counter = 0
    counter2 = 0

    for aa ∈ α
        counter2 += 1
        for pp ∈ ϕ

            Aproj, Aref, OutLMNTs, areas_and_angles = fIlluminationNonConvex(triangles, aa, pp)
            counter += 1
            storageValues = SweepStorage(rad2deg(aa), rad2deg(pp), Aproj, Aref, OutLMNTs)
            AlphaPhiStorage[counter, counter2] = storageValues

        end
        counter = 0
    end

    return AlphaPhiStorage

end


#TEST
#-------------------------------------
# MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]
# step = pi / 20
# aa = pi
# pp = pi / 2


# AlphaPhiStorage = sweep(MeshVerticesCoords, step)


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