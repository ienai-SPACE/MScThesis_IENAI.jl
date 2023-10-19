struct SolarCellGeometry{F<:FaceGeometry,T}
    geometry::Vector{F}
    Aproj::T
end
struct _SolarCellGeometry{T}
    geometry::T
    Aproj::T
end

function getSolarCellsGeo(geometry::Geometry, viewpoint, sampler)

    #solar cell identification
    geo_solarCell = filter_solarCell(geometry)
    solarCell_geometry = [geo_solarCell.faces[ii].geometry for ii ∈ 1:lastindex(geo_solarCell.faces)]

    #solar cells: back-face culling + ray-tracing
    rti_vec_sc, filtered_geometry_sc, Aray_sc = raytrace(geo_solarCell, viewpoint, sampler) #culling +  ray tracing
    if rti_vec_sc == 0
        Aproj_sc = 0.0
    else
        valid_rti_sc = rti_vec_sc |> Filter(rti -> rti.mode == FrontFaceIntersection) |> tcollect
        Aproj_sc = Aray_sc * length(valid_rti_sc)
    end

    SolarCellGeometry(solarCell_geometry, Aproj_sc)
end

function getSolarCellsGeo(geometry::HomogeneousGeometry, viewpoint, sampler)
    return _SolarCellGeometry(0, 0)
end
# function getSolarCellsGeo(geometry::HomogeneousGeometry)
#     return 0
# end

function getSolarCellsGeo(geometry::AbstractGeometry)
    if is_homogeneous(geometry) != true
        #solar cell identification
        filtered_geo_solarCell = filter_solarCell(geometry)
        solarCell_geometry = [filtered_geo_solarCell.faces[ii].geometry for ii ∈ 1:lastindex(filtered_geo_solarCell.faces)]
        return SolarCellGeometry(solarCell_geometry, 0)
    else
        solarCell_geometry = 0
        return _SolarCellGeometry(solarCell_geometry, 0)
    end
end