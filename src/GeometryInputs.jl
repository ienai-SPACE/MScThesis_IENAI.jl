using FileIO
import FilePathsBase
using FilePathsBase: /
# using GeometryBasics

@enum GeometryType Convex NonConvex

abstract type AbstractGeometry{GT} end

"""
    Geometry{GT,F<:Face} <: AbstractGeometry{GT}

- `faces::Vector{F}`
- `rmax::Float64` 
"""
struct Geometry{GT,F<:Face} <: AbstractGeometry{GT}
    faces::Vector{F}
    rmax::Float64
    function Geometry{GT}(faces::Vector{F}) where {GT,F}
        rmax = maximum([_max_coord_euclidean(face) for face in faces])
        new{GT,F}(faces, rmax)
    end
end

"""
    HomogeneousGeometry{GT,F<:FaceGeometry} <: AbstractGeometry{GT}

- `faces::Vector{F}`
- `surface_props::SurfaceProps{T}`
- `rmax::Float64`
"""
struct HomogeneousGeometry{GT,F<:FaceGeometry} <: AbstractGeometry{GT}
    faces::Vector{F}
    rmax::Float64
    function HomogeneousGeometry{GT}(faces::Vector{F}) where {GT,F}
        rmax = maximum([_max_coord_euclidean(face) for face in faces])
        new{GT,F}(faces, rmax)
    end
end

is_convex(geo::AbstractGeometry{Convex}) = true
is_convex(geo::AbstractGeometry{NonConvex}) = false
is_homogeneous(geo::HomogeneousGeometry) = true
is_homogeneous(geo::Geometry) = false

_max_coord(face::Face) = _max_coord(face.geometry)
_max_coord(face::TriangleFace) = maximum([maximum(v) for v ∈ face.vertices])
_max_coord_euclidean(face::TriangleFace) = maximum([norm(v) for v ∈ face.vertices])
_max_coord_euclidean(face::Face) = _max_coord_euclidean(face.geometry)

get_point_data(point) = point.main.data
get_point_data(point::Point) = point.data

"""
    n_faces(geo::AbstractGeometry)

Return the number of faces in `geo`
"""
n_faces(geo::AbstractGeometry) = length(geo.faces)

unit_dict = Dict("m" => 1, "mm" => 1e-3)

"""
    load_geometry(path, is_convex::Bool, units::String)

# Inputs
- `path
- `surface_props::SurfaceProps`
- `is_convex::Bool`
- `units::String`
# Outputs
- `HomogeneousGeometry{gtype}(faces, surface_props)`
"""
function load_geometry(path, is_convex::Bool, units::String)
    scale_factor = unit_dict[units]
    mesh = load(path)
    faces = TriangleFace{Float64}[]
    for triangle in mesh
        points = map(i -> get_point_data(triangle.points[i]) |> SVector{3,Float64} |> p -> p * scale_factor, 1:3)
        face_geometry = TriangleFace(points...)
        push!(faces, face_geometry)
    end
    gtype = is_convex ? Convex : NonConvex
    HomogeneousGeometry{gtype}(faces)
end

"""
    load_geometry(path, materials_path, is_convex::Bool, units::String)

# Inputs
- `path
- `surface_props::SurfaceProps`
- `is_convex::Bool`
- `units::String`

# Outputs
- `Geometry{gtype}(faces)`
"""
function load_geometry(path, materials_path, is_convex::Bool, units::String)
    mesh_facets, materials = load_material_properties(materials_path)
    materials = [SurfaceAtomProperties(; atomic_mass=material["atomic mass"]) |> SurfaceProps for material in materials]
    scale_factor = unit_dict[units]
    mesh = load(path)
    faces = Face{TriangleFace{Float64},Float64}[]
    for (triangle, mat_index) in zip(mesh, mesh_facets)
        points = map(i -> get_point_data(triangle.points[i]) |> SVector{3,Float64} |> p -> p * scale_factor, 1:3)
        face_geometry = TriangleFace(points...)
        idx = mat_index["material_index"] + 1
        face = Face(face_geometry, materials[idx])
        push!(faces, face)
    end
    gtype = is_convex ? Convex : NonConvex
    indices = collect(1:length(faces))
    Geometry{gtype}(faces)
end



"""
    Viewpoint(geo::AbstractGeometry, azimuth, elevation)

# Inputs
- `geo::AbstractGeometry`
- `azimuth`
- `elevation`
# Outputs
- `Viewpoint(rmax, distance, azimuth, elevation)`
"""
function Viewpoint(geo::AbstractGeometry, azimuth, elevation)
    rmax = get_rmax(geo)
    distance = 100 * rmax
    Viewpoint(rmax, distance, azimuth, elevation)
end

"""
    Viewpoint(geo::AbstractGeometry, , Vdir)

# Inputs
- `geo::AbstractGeometry`
- `Vdir:: Vector`
# Outputs
- `Viewpoint(rmax, distance, azimuth, elevation)`
"""
function Viewpoint(geo::AbstractGeometry, Vdir)
    rmax = get_rmax(geo)
    distance = 100 * rmax
    Viewpoint(rmax, distance, Vdir)
end

get_rmax(geo::AbstractGeometry) = geo.rmax


export load_geometry

include("inputMesh.jl")

"""
    GeomInputs(Vrel_v, VdirFlag, convexFlag)

Define whether the object is convex or non-convex
Define if the direction of "particle impingement" is set manually or by the velocity vector
Load the the meshed element

# Input:
- `Vrel_v`      : [m/s] relative velocity vector
- `VdirFlag`   : flag for direction criterion
- `convexFlag` : flag = 1 for convex shape and flag = 0 for non-convex
# Output:
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
    mesh = load(pkg_path / "test" / "inputs_models_data" / "T_Sat_fineMesh.obj") #T_SatMesh  sphereMesh cubeMesh coneMesh T_Sat_fineMesh


    MeshVerticesCoords = finputMesh(mesh)
    # MeshVerticesCoords = Float64.([1 1 0 0 1 1 1 0 1; 1 1 0 1 0 -1 0 1 -1; 0.5 0.5 0 1 1 1 0 0 1])

    if convexFlag == 0
        maxDistance = euclidean_norm(MeshVerticesCoords)
        # rmax = maximum(MeshVerticesCoords[:, :])          #radius of the circular plane from where rays originate
        rmax = maxDistance
        println("rmax = ", rmax)
        distance = rmax * 100.0                           #distance at which the circular plane is located (it should be out from the satellite body)
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

"""
    euclidean_norm(MeshVerticesCoords::Matrix{Float64})

Calculate the euclidean norm from the datum of the body reference frame to the furthest point in the object

# Input
- `MeshVerticesCoords::Matrix{Float64}`
# Output
- 'rmax = maximum(_maxVertex)'
"""
function euclidean_norm(MeshVerticesCoords::Matrix{Float64})
    Ntri = lastindex(MeshVerticesCoords[:, 1])
    maxVertex = [0.0, 0.0, 0.0]
    _maxVertex = Vector{Float64}(undef, Ntri)
    for jj ∈ 1:Ntri
        for ii ∈ 0:2
            maxVertex[ii+1] = sqrt(MeshVerticesCoords[jj, 1+3*ii]^2 + MeshVerticesCoords[jj, 2+3*ii]^2 + MeshVerticesCoords[jj, 3+3*ii]^2)
        end
        _maxVertex[jj] = maximum(maxVertex)
    end
    maximum(_maxVertex)
end

"""
    angle_V_n(dir, normal)

Calculate the angle between two vectors (the velocity direction vector and the normal vector of the surface)

# Input
- `dir`
- `normal`          :either scalar or vector may be used
# Output
- `gamma = acos.(dp_v) `
"""
function angle_V_n(dir, normal)
    Vdir = -dir
    dp_v = [dot(Vdir, normal[ii]) / (norm(normal[ii]) * norm(Vdir)) for ii in 1:Int(length(normal))]
    acos.(dp_v)
end

"""
    filter_backfaces(geometry::HomogeneousGeometry{GT}, viewpoint::Viewpoint) where {GT}

Filter out all the non-forward facing face_vertices

# Input
- `geometry::HomogeneousGeometry{GT}`
- `viewpoint::Viewpoint`
# Outputs
- `HomogeneousGeometry{GT}(new_faces, new_indices)`
"""
function filter_backfaces(geometry::HomogeneousGeometry{GT}, viewpoint::Viewpoint) where {GT}
    new_faces = filter(face -> is_visible(face, viewpoint), geometry.faces)
    HomogeneousGeometry{GT}(new_faces)
end


"""
    filter_backfaces(geometry::Geometry{GT}, viewpoint::Viewpoint) where {GT}

Filter out all the non-forward facing face_vertices

# Input
- `geometry::Geometry{GT}`
- `viewpoint::Viewpoint`
# Outputs
- `Geometry{GT}(new_faces, new_indices)`
"""
function filter_backfaces(geometry::Geometry{GT}, viewpoint::Viewpoint) where {GT}
    new_faces = filter(face -> is_visible(face, viewpoint), geometry.faces)
    Geometry{GT}(new_faces)
end


"""
    is_visible(face::TriangleFace, viewpoint::Viewpoint)

Forward-facing check

# Input
- `face::TriangleFace`
- `viewpoint::Viewpoint`
# Output
- ::Bool
"""
function is_visible(face::TriangleFace, viewpoint::Viewpoint)
    normal = face.normal
    #viewpoint direction is in the same sense as oncoming ray beam
    dot(normal, viewpoint.direction) < -1e-3
end

is_visible(face::Face, viewpoint) = is_visible(face.geometry, viewpoint)

"""
    shrink_viewpoint(geom::AbstractGeometry, viewpoint::Viewpoint)
    
It calls the function `Viewpoint(geom.rmax, 100geom.rmax, viewpoint.direction)` to adapt the viewpoint to the filtered geometry
"""
shrink_viewpoint(geom::AbstractGeometry, viewpoint::Viewpoint) = Viewpoint(geom.rmax, 100geom.rmax, viewpoint.direction)
face_area(geometry::HomogeneousGeometry, idx) = geometry.faces[idx].area
face_vertices(geometry::HomogeneousGeometry, idx) = geometry.faces[idx].vertices
face_normal(geometry::HomogeneousGeometry, idx) = geometry.faces[idx].normal
# face_idx(geometry::HomogeneousGeometry, idx) = geometry.indices[idx]
face_area(geometry::Geometry, idx) = geometry.faces[idx].geometry.area
face_vertices(geometry::Geometry, idx) = geometry.faces[idx].geometry.vertices
face_normal(geometry::Geometry, idx) = geometry.faces[idx].geometry.normal


"""
    Grid{T}
    
- `alpha::Vector{T}`
- `phi::Vector{T}`
"""
struct Grid{T}
    alpha::Vector{T}
    phi::Vector{T}

    function Grid(step::T) where {T<:Real}
        α = collect(-π:step:π)
        ϕ = collect(-π/2:step:π/2)  #APPLY SYMMETRY IF POSSIBLE
        return new{T}(α, ϕ)
    end
end

"""
    interpolator(_look_up_table::Matrix, aa_in, pp_in)

Interpolator to access the look-up table

# Input
- `_look_up_table::Matrix`
- `aa_in`               : alpha [deg] value to be accessed in the look-up table
- `pp_in`               : phi [deg] value to be accessed in the look-up table
# Output
- `P(pp_in, aa_in)`     : Interpolated value of the look-up table
"""
function interpolator(_look_up_table::Matrix, aa_in, pp_in)
    _phi = _look_up_table[:, 1] |> length |> ss -> (180 / (ss - 1)) |> (step -> collect(-90:step:90))
    _alpha = _look_up_table[1, :] |> length |> ss -> (360 / (ss - 1)) |> (step -> collect(-180:step:180))
    P = BicubicInterpolator(_phi, _alpha, _look_up_table, StrictBoundaries())
    return P(pp_in, aa_in)  #[deg] P(phi, alpha)
end


"""
    face_barycenters(filtered_geo::HomogeneousGeometry, idx)

# Input
- `filtered_geo::HomogeneousGeometry` 
- `idx`
# Outputs
- Barycenter of the input triangle
"""
function face_barycenters(filtered_geo::HomogeneousGeometry, idx)
    vtx = filtered_geo.faces[idx].vertices
    (1 / 3) * (vtx[1] + vtx[2] + vtx[3])
end

"""
    face_barycenters(filtered_geo::Geometry, idx)

# Input
- `filtered_geo::Geometry` 
- `idx`
# Outputs
- Barycenter of the input triangle
"""
function face_barycenters(filtered_geo::Geometry, idx)
    vtx = filtered_geo.faces[idx].geometry.vertices
    (1 / 3) * (vtx[1] + vtx[2] + vtx[3])
end

"""
    getBaryceners(filtered_geo::HomogeneousGeometry, hit_idx)

# Input
- `filtered_geo::AbstractGeometry` 
- `hit_idx`
# Outputs
- Vector of barycenters of the filtered geometry
"""
function getBarycenters(filtered_geo::AbstractGeometry, hit_idx)
    b = Vector{Vector}(undef, length(hit_idx))
    for ii in 1:lastindex(hit_idx)
        b[ii] = face_barycenters(filtered_geo, hit_idx[ii])
    end
    return b
end

