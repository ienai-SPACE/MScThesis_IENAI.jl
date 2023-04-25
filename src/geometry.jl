abstract type FaceGeometry end

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
        normal = crossProd / (2area)
        new{T}(vertices, area, normal)
    end
end
