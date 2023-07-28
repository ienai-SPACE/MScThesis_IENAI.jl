"""
    Ray{T}

- `origin::SV3{T}`
- `direction::SV3{T}`
"""

struct Ray{T}
    origin::SV3{T}
    direction::SV3{T}
    function Ray(origin::SV3{T}, direction::SV3{T}) where {T}
        new{T}(origin, normalize(direction))
    end
end

@enum IntersectionMode begin
    NoIntersection
    FrontFaceIntersection
    BackFaceIntersection
end

"""
    RayTriangleIntersection{T}

- `mode::IntersectionMode`
- `t::T`
- `γ_dir::T`
- `area::T`
- `face_index::Int64`
"""

struct RayTriangleIntersection{T}
    mode::IntersectionMode
    t::T
    γ_dir::T
    area::T
    face_index::Int64
end

"""
    earlier_intersection(rti1::RayTriangleIntersection, rti2::RayTriangleIntersection)

Select which of the two input magnitudes is the biggest one

#INPUTS
- `rti1::RayTriangleIntersection`
- `rti2::RayTriangleIntersection`
#OUTPUTS
- `rti1::RayTriangleIntersection`
- `rti2::RayTriangleIntersection`
"""

function earlier_intersection(rti1::RayTriangleIntersection, rti2::RayTriangleIntersection)
    if rti1.mode == NoIntersection
        return rti2
    elseif rti2.mode == NoIntersection
        return rti1
    else
        return rti1.t < rti2.t ? rti1 : rti2
    end
end

mode(rti::RayTriangleIntersection) = rti.mode

RayTriangleIntersection(mode::IntersectionMode, t::T, γ_dir::T, area::T, i::Integer) where {T} = RayTriangleIntersection{T}(mode, t, γ_dir, area, i)
no_intersection(::Type{T}) where {T} = RayTriangleIntersection(NoIntersection, zero(T), zero(T), zero(T), 0)

"""
    MTalgorithm(triangle::TriangleFace{T}, ray::Ray{T}; ϵ=sqrt(eps(T))) where {T}

Identify the triangular mesh elements that are intercepted by the direction vector through a ray-tracing algorithm (Möller-Trumbore algorithm)
Also, calculate the area of the triangular element and angle between the triangle's normal and the direction (velocity) vector is also provided

#INPUT
- `triangle::TriangleFace{T}`       : it stores `vertices::SVector{3, SV3{T}}`, `area::T` and `normal::SVector{3, T}` fields
- `ray::Ray{T}`                     : it stores the ray information with fields `direction` and `origin`
#OUTPUT
- `RayTriangleIntersection(which_face, t, γ_dir, triangle.area, 0)`     : it contains the fields `mode::IntersectionMode`, `t::T`, `γ_dir::T` [rad], `area::T` [m^2], `face_index::Int64`
"""

function MTalgorithm(triangle::TriangleFace{T}, ray::Ray{T}; ϵ=sqrt(eps(T))) where {T}

    Vdir = -ray.direction
    # println("dir in MT", dir)

    dot_prod = dot(triangle.normal, Vdir)

    γ_dir = acos(dot_prod / (norm(triangle.normal) * norm(Vdir)))
    # println(dot_prod, rad2deg(γ_dir))
    if abs(dot_prod) < ϵ
        return no_intersection(T)
    end
    which_face = dot_prod > 1e-5 ? FrontFaceIntersection : BackFaceIntersection

    V = triangle.vertices
    edge1 = V[2] - V[1]
    edge2 = V[3] - V[1]
    pvec = cross(ray.direction, edge2)
    det = dot(edge1, pvec)
    invDet = 1 / det
    tvec = ray.origin - V[1]
    u = invDet * dot(tvec, pvec)
    #calculation of triangle's area

    (u < 0.0 || u > 1) && return no_intersection(T)

    qvec = cross(tvec, edge1)
    v = invDet * dot(ray.direction, qvec)

    (v < 0.0 || u + v > 1.0) && return no_intersection(T)
    t = invDet * dot(edge2, qvec) # distance from the ray origin to P 
    RayTriangleIntersection(which_face, t, γ_dir, triangle.area, 0)
end

MTalgorithm(face::Face, ray) = MTalgorithm(face.geometry, ray)

