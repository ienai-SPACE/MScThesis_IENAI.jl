using LinearAlgebra
using StaticArrays


"""
    fIlluminationConvex(triangles, α, ϕ)

#INPUT:
- `triangles::SMatrix{number of triangles, 9, Float64, 9*number of triangles}` : coordinates of facets' vertices
- `α::Float64`            : azimuth w.r.t. body frame. (α = 0 aligned with the x-axis)
- `ϕ::Float64`            : elevation w.r.t. body frame. (ϕ = 0 on the xy-plane)       
#OUTPUT:
- `OutLMNTs::OutGeometry{T}`    : struct with fields `area::Vector{T}` and `angles::Vector{T}`
- `Aproj::Float64`              : projection of the intercepted triangular areas onto the velocity direction 
- `Aref::Float64`               : sum of all intercepted triangular areas
- `areas_and_angles::Vector{impinged_geometries{T}}` : vector storing struct with fields 
     
"""

struct impinged_geometries{T}
    area::T
    angle::T
end

struct auxOutGeometry{T}
    area::Vector{T}
    angle::Vector{T}
end


_eltype(::impinged_geometries{T}) where {T} = T

function fIlluminationConvex(triangles, α, ϕ)


    #from spherical to cartesian coordinates
    ux = cos(ϕ) * cos(α)
    uy = cos(ϕ) * sin(α)
    uz = sin(ϕ)
    dir = [ux, uy, uz] #unitary directional vector (opposite sense to oncoming ray beam)
    # print(dir)
    # print("illumination")
    #dir = [1 1 0]
    #Number of triangles/facets
    Nfacet = size(triangles, 1)


    for jj ∈ 1:Nfacet #for loop to iterate over all triangles/quads
        vertices = triangles[jj, 1:9]
        OutAreaConvex = areasConvex(vertices, dir)
        if jj == 1
            if OutAreaConvex[1] == 0
                OutFacets = @SVector [0.0, 0.0, 0.0]
            else
                OutFacets = @SVector [jj, OutAreaConvex[1], OutAreaConvex[2]]
            end
        else
            if OutAreaConvex[1] == 0
                OutFacets = hcat(OutFacets, [0.0, 0.0, 0.0])
            else
                OutFacets = hcat(OutFacets, [jj, OutAreaConvex[1], OutAreaConvex[2]])
            end
        end
    end

    #eliminate non-intercepted triangle entries and size down the output matrix
    OutFacets = filter(!iszero, OutFacets)
    OutFacets = reshape(OutFacets, (3, Int(length(OutFacets) / 3)))

    # if iszero(OutFacets)
    #     OutFacets = [0.0; 0.0; 0.0]
    # end

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


    OutLMNTs = auxOutGeometry(OutFacets[2, :], OutFacets[3, :]) #previously: outGeometry(OutFacets[2, :], OutFacets[3, :])


    imp_geo = impinged_geometries(OutFacets[2, :], OutFacets[3, :])
    T = _eltype(imp_geo)

    #area_and_angle = Vector{InteractionGeometry{T}}(undef, length(OutFacets[2, :]))
    areas_and_angles = map(ii -> impinged_geometries(OutFacets[2, ii], OutFacets[3, ii]), 1:length(OutFacets[2, :]))



    if Aproj == 0
        Aproj = 0.0
    end
    if Aref == 0
        Aref = 0.0
    end



    return Aproj, Aref, OutLMNTs, areas_and_angles


end


#TESTING
#-------------------------------------------------------
#=
α = deg2rad(25)
ϕ = deg2rad(250)

Aproj, Aref, OutFacets = fIlluminationConvex(α, ϕ)
=#






