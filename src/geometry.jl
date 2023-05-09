abstract type FaceGeometry end

"""
    TriangleFace{T} <: FaceGeometry 
    
Calculate the area and the normal of a 3D triangular surface defined by three static vectors (one per vertex)

- `vertices::SVector{3,SV3{T}}`
- `area::T`
- `normal::SVector{3,T}`
"""

struct TriangleFace{T} <: FaceGeometry
    vertices::SVector{3,SV3{T}}
    area::T
    normal::SVector{3,T}
    function TriangleFace(v1::SV3{T}, v2::SV3{T}, v3::SV3{T}) where {T}
        vertices = SVector(v1, v2, v3)
        edge1 = v2 - v1
        edge2 = v3 - v2
        crossProd = cross(edge1, edge2)
        area = norm(crossProd) / 2
        if area == 0
            normal = zeros(3)
        elseif area != 0
            normal = crossProd / (2area)
        end

        new{T}(vertices, area, normal)
    end
end

#=
v3 = @SVector[1, 2, 5]
v2 = @SVector[1, 2, 3]
v1 = @SVector[1, 2, 4]

TriangleFace(v1, v2, v3)
=#