using FileIO
import FilePathsBase
using FilePathsBase: /
# using GeometryBasics

@enum GeometryType Convex NonConvex

abstract type AbstractGeometry{GT} end

struct Geometry{GT,F<:Face} <: AbstractGeometry{GT}
    faces::Vector{F}
    rmax::Float64
    function Geometry{GT}(faces::Vector{F}) where {GT,F}
        rmax = maximum([_max_coord(face) for face in faces])
        new{GT,F}(faces, rmax)
    end
end

struct HomogeneousGeometry{GT,T,F<:FaceGeometry} <: AbstractGeometry{GT}
    faces::Vector{F}
    surface_props::SurfaceProps{T}
    rmax::Float64
    function HomogeneousGeometry{GT}(faces::Vector{F}, surface_props::SurfaceProps{T}) where {GT,F,T}
        rmax = maximum([_max_coord(face) for face in faces])
        new{GT,T,F}(faces, surface_props, rmax)
    end
end

is_convex(geo::AbstractGeometry{Convex}) = true
is_convex(geo::AbstractGeometry{NonConvex}) = false

_max_coord(face::Face) = _max_coord(face.geometry)
_max_coord(face::TriangleFace) = maximum([maximum(v) for v ∈ face.vertices])

get_point_data(point) = point.main.data
get_point_data(point::Point) = point.data

n_faces(geo::AbstractGeometry) = length(geo.faces)


function load_geometry(path, surface_props::SurfaceProps, is_convex::Bool, scale_factor=1e-3)
    mesh = load(path)
    faces = TriangleFace{Float64}[]
    for triangle in mesh
        points = map(i -> get_point_data(triangle.points[i]) |> SVector{3,Float64} |> p -> p * scale_factor, 1:3)
        face_geometry = TriangleFace(points...)
        push!(faces, face_geometry)
    end
    gtype = is_convex ? Convex : NonConvex
    HomogeneousGeometry{gtype}(faces, surface_props)
end

function Viewpoint(geo::AbstractGeometry, azimuth, elevation)
    rmax = get_rmax(geo)
    distance = 100 * rmax
    Viewpoint(rmax, distance, azimuth, elevation)
end

get_rmax(geo::AbstractGeometry) = geo.rmax

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
- `Vrel_v`      : [m/s] relative velocity vector
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

function GeomInputs(Vrel_v::Vector{Float64}, VdirFlag::Int64, convexFlag::Int64)



    #load the mesh
    pkg_path = FilePathsBase.@__FILEPATH__() |> parent |> parent
    mesh = load(pkg_path / "test" / "samples" / "T_SatMesh.obj") #T_SatMesh  sphereMesh cubeMesh coneMesh T_Sat_fineMesh


    MeshVerticesCoords = finputMesh(mesh)
    # MeshVerticesCoords = Float64.([1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1])

    if convexFlag == 0
        rmax = maximum(MeshVerticesCoords[:, :])          #radius of the circular plane from where rays originate
        distance = rmax * 100.0                                                    #distance at which the circular plane is located (it should be out from the satellite body)
    else
        rmax = 0.0
        distance = 0.0
    end

    if VdirFlag == 1
        #direction vector
        dir = SV3(Vrel_v / norm(Vrel_v))   #direction defined in the opposite direction to velocity
    elseif VdirFlag == 0
        dir = Float64.(@SVector [1, 0, 0])
    end


    return MeshVerticesCoords, dir, rmax, distance

end

function filter_backfaces(geometry::HomogeneousGeometry{GT}, viewpoint::Viewpoint) where {GT}
    new_faces = filter(face -> is_visible(face, viewpoint), geometry.faces)
    HomogeneousGeometry{GT}(new_faces, geometry.surface_props)
end

function is_visible(face::TriangleFace, viewpoint::Viewpoint)
    normal = face.normal
    dot(normal, viewpoint.direction) < 0
end

shrink_viewpoint(geom::AbstractGeometry, viewpoint::Viewpoint) = Viewpoint(geom.rmax, 100geom.rmax, viewpoint.direction)
face_area(geometry::HomogeneousGeometry, idx) = geometry.faces[idx].area