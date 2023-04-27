
"""        function: Möller-Trumbore algorithm

GOAL:
    - Identify the triangular mesh elements that are intercepted by the direction vector through a ray-tracing algorithm (Möller-Trumbore algorithm)
    - Area of the triangular element and angle between the triangle's normal and the direction (velocity) vector is also provided

INPUTS:
     - V1, V2, V3 :: vertices of the triangle
     - dir        :: direction of the ray
     - O          :: origing of the ray
OUTPUTS:
      - flag      :: 1 for intercept, 0 for NO intercept
      - facing    :: 1 for frontfacing, 0 for backfacing
      - u, v      :: barycentric coordinates
      - t         :: distance from the ray origin to P
      - area      :: [m^2] area of the triangular element
      - γ_dir     :: [rad] angle between the triangle's normal and the direction (velocity) vector is also provided
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

struct RayTriangleIntersection{M,T}
    t::T
    γ_dir::T
end

RayTriangleIntersection{M}(t::T, γ_dir::T) where {M,T} = RayTriangleIntersection{M,T}(t, γ_dir)
no_intersection(::Type{T}) where {T} = RayTriangleIntersection{NoIntersection}(zero(T), zero(T))


function MTalgorithm(triangle::TriangleFace{T}, ray::Ray{T}; ϵ=sqrt(eps(T))) where {T}
    dot_prod = dot(triangle.normal, ray.direction)
    γ_dir = acos(dot_prod)
    if abs(dot_prod) < ϵ
        return no_intersection(T) # TODO: refine 
    end
    which_face = dot_prod > 0.0 ? BackFaceIntersection : FrontFaceIntersection

    V = triangle.vertices
    edge1 = V[2] - V[1]
    edge2 = V[3] - V[1]
    pvec = cross(ray.direction, edge2)
    det = dot(edge1, pvec)
    invDet = 1 / det
    tvec = ray.origin - V[1]
    u = invDet * dot(tvec, pvec)

    (u < 0.0 || u > 1) && return no_intersection(T)

    qvec = cross(tvec, edge1)
    v = invDet * dot(dir, qvec)

    (v < 0.0 || u + v > 1.0) && return no_intersection(T)
    t = invDet * dot(edge2, qvec) # distance from the ray origin to P 
    RayTriangleIntersection{which_face}(t, γ_dir)
end


#=

function MTalgorithm(O, dir, V1, V2, V3)

    #defined for tolerance purposes
    dummyFlag = 0                   #flag to indicate: 0-> intercept; 1-> no intercept

    edge1 = V2 - V1
    edge2 = V3 - V1
    pvec = cross(dir, edge2)        #perpendicular vector
    crossProd = cross(edge1, edge2)
    area = norm(crossProd) / 2        #area of the triangle
    det = dot(edge1, pvec)          #determinant of the matrix M

    dot_prod = dot(crossProd, dir)

    #Calculate the angle between triangle's normal and velocity direction 
    γ_dir = acos(dot_prod / (norm(crossProd) * norm(dir)))    #[radians]


    #determine side from which ray incides on the triangle: CULLING 

    if det > 0
        facing = 1 #front-facing
    else
        facing = 0 #back-facing
    end


    if (det > -ϵ && det < ϵ)
        # the vector is parallel to the plane (the intersection is at infinity)
        flag = 0
        u = 0
        v = 0
        t = 0
        dummyFlag = 1
    end

    if dummyFlag == 0

        invDet = 1 / det
        tvec = O - V1
        u = invDet * dot(tvec, pvec)


        if (u < 0.0 || u > 1)
            # the ray does not intersect the triangle
            flag = 0
            u = 0
            v = 0
            t = 0
            dummyFlag = 1
        end


        if dummyFlag == 0
            qvec = cross(tvec, edge1)
            v = invDet * dot(dir, qvec)

            if (v < 0.0 || u + v > 1.0)
                # the ray does not intersect the triangle
                flag = 0
                u = 0
                v = 0
                t = 0
                dummyFlag = 1
            end

            if dummyFlag == 0
                t = invDet * dot(edge2, qvec) # distance from the ray origin to P 
                flag = 1                   # the ray intersects the triangle
            end

        end
    end


    return area, flag, facing, u, v, t, γ_dir

end
=#