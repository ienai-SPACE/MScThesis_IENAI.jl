using LinearAlgebra
using StaticArrays
include("trialGeometry.jl")
include("Convex.jl")


"""
function : fIlluminationConvex(α, ϕ)

    INPUT:
        - gas properties??
        - include      : Geometry.jl
            - Nfacets  : number of facets
            - MeshVerticesCoords : coordinates of facets' vertices
        - α            : azimuth w.r.t. body frame. (α = 0 aligned with the x-axis)
        - ϕ            : elevation w.r.t. body frame. (ϕ = 0 on the xy-plane)       
    OUTPUT:
        - OutFacets    :: matrix storing the triangle index, area, and the angle between the normal and the velocity direction
            - Matrix{Float64} -> size: (3, number of intercepted triangles)
        - Aproj           :: projection of the intercepted triangular areas onto the velocity direction 
        - Aref            :: sum of all intercepted triangular areas
"""


function fIlluminationConvex(triangles, α, ϕ)



    #from spherical to cartesian coordinates
    ux = cos(ϕ) * cos(α)
    uy = cos(ϕ) * sin(α)
    uz = sin(ϕ)
    dir = [ux, uy, uz] #unitary directional vector (opposite sense to oncoming ray beam)


    #Number of triangles/facets
    Nfacet = size(triangles, 1)

    #loop to iterate over all triangles
    for jj ∈ 1:Nfacet
        MeshVerticesCoords = triangles[jj, 1:9]
        OutAreaConvex = areasConvex(MeshVerticesCoords, dir)
        if jj == 1
            if OutAreaConvex[1] == 0
                OutFacets = @SVector [0, 0, 0]
            else
                OutFacets = @SVector [jj, OutAreaConvex[1], OutAreaConvex[2]]
            end
        else
            if OutAreaConvex[1] == 0
                OutFacets = hcat(OutFacets, [0, 0, 0])
            else
                OutFacets = hcat(OutFacets, [jj, OutAreaConvex[1], OutAreaConvex[2]])
            end
        end
    end

    #eliminate non-intercepted triangle entries and size down the output matrix
    OutFacets = filter(!iszero, OutFacets)
    OutFacets = reshape(OutFacets, (3, Int(length(OutFacets) / 3)))



    #sum of all intercepted triangular areas
    Aref = sum(OutFacets[2, :])


    #Project all intercepted triangular areas and find projected total area

    #------pre-allocation-------------------
    Aproj = zeros(size(OutFacets, 2), 1)
    #---------------------------------------

    for ii ∈ 1:Int(size(OutFacets, 2))
        Aproj[ii] = abs(OutFacets[2, ii] * cos(OutFacets[3, ii]))
    end
    Aproj = sum(Aproj)                 #sum of all intercepted triangular projected areas


    return Aproj, Aref, OutFacets


end


#TESTING
#-------------------------------------------------------
α = deg2rad(25)
ϕ = deg2rad(250)

Aproj, Aref, OutFacets = fIlluminationConvex(α, ϕ)






