using FileIO
import FilePathsBase
using FilePathsBase: /
# using GeometryBasics

@enum GeometryType Convex NonConvex

struct Geometry{GT,F<:Face}
    faces::Vector{F}
    rmax::Float64
    function Geometry{GT}(faces::Vector{F}) where {GT,F}
        rmax = maximum([_max_coord(face) for face in faces])
        new{GT,F}(faces, rmax)
    end
end

struct HomogeneousGeometry{GT,T,F<:FaceGeometry}
    faces::Vector{F}
    surface_props::SurfaceProps{T}
    rmax::Float64
    function HomogeneousGeometry{GT}(faces::Vector{F}, surface_props::SurfaceProps{T}) where {GT,F,T}
        rmax = maximum([_max_coord(face) for face in faces])
        new{GT,T,F}(geometry, surface_props, rmax)
    end
end

_max_coord(face::Face) = _max_coord(face.geometry)
_max_coord(face::TriangleFace) = maximum([maximum(v) for v ∈ face.vertices])

function load_geometry(path, surface_props::SurfaceProps, is_convex::Bool)
    mesh = load(path)
    faces = Face{TriangleFace{Float64},Float64}[]
    for triangle in mesh
        points = map(i -> triangle.points[i].data |> SVector{3,Float64}, 1:3)
        face_geometry = TriangleFace(points...)
        push!(faces, face_geometry)
    end
    gtype = is_convex ? Convex : NonConvex
    HomogeneousGeometry{gtype}(faces, surface_props)
end

# regular Geometry
# function load_geometry(path, surface_props::SurfaceProps, is_convex::Bool)
#     mesh = load(path)
#     faces = Face{TriangleFace{Float64},Float64}[]
#     for triangle in mesh
#         points = map(i -> triangle.points[i].data |> SVector{3,Float64}, 1:3)
#         geometry = TriangleFace(points...)
#         push!(faces, Face(geometry, surface_props))
#     end
#     gtype = is_convex ? Convex : NonConvex
#     Geometry{gtype}(faces)
# end

export load_geometry

include("inputMesh.jl")

"""
    GeomInputs(Vrel_v, VdirFlag, convexFlag)

Define whether the object is convex or non-convex
Define if the direction of "particle impingement" is set manually or by the velocity vector
Load the the meshed element

#INPUT:
-`Vrel_v`      : [m/s] relative velocity vector
- `VdirFlag`   : flag for direction criterion
- `convexFlag` : flag = 1 for convex shape and flag = 0 for non-convex
#OUTPUT:
- `MeshVerticesCoords`
- `dir`
- `rmax`, distance
        - MeshVerticesCoords : vertices coordinates of the meshed element
        - rmax               : radius of the source for ray-tracing algorithm
        - distance           : distance between source origin and local coordinate system for ray-tracing algorithm
"""

function GeomInputs(Vrel_v, VdirFlag, convexFlag)

    # MeshVerticesCoords = @SMatrix [1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1]

    #load the mesh
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh = load(pkg_path / "test" / "samples" / "sphereMesh4.obj")

    MeshVerticesCoords = finputMesh(mesh)

    if convexFlag == 0
        rmax = maximum(MeshVerticesCoords[:, :])          #radius of the circular plane from where rays originate
        distance = rmax * 10                                                    #distance at which the circular plane is located (it should be out from the satellite body)
    else
        rmax = 0
        distance = 0
    end

    if VdirFlag == 1
        #direction vector
        dir = (Vrel_v / norm(Vrel_v))'
    elseif VdirFlag == 0
        dir = @SVector [1 / sqrt(2), 1 / sqrt(2), 0]
    end


    return MeshVerticesCoords, dir, rmax, distance

end