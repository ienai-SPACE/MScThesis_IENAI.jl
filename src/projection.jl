"""
    projection(geo::AbstractGeometry, viewpoint::Viewpoint, Ntri)

Projection of areas according to a viewpoint direction

# Inputs
- `geo::AbstractGeometry`
- `viewpoint::Viewpoint`
- `Ntri`
# Outputs
- vector of projected areas
"""
function projection(geo::HomogeneousGeometry, viewpoint::Viewpoint, Ntri)

    m = Map(jj -> begin
        n = geo.faces[jj].normal
        Vdir = -viewpoint.direction
        if (norm(n) * norm(Vdir)) == 0.0
            c = 0.0
        else
            c = dot(n, Vdir) / (norm(n) * norm(Vdir))
        end
        Aproj_i = c * geo.faces[jj].area
    end)

    return 1:Ntri |> m |> collect
end


function projection(geo::Geometry, viewpoint::Viewpoint, Ntri)

    m = Map(jj -> begin
        n = geo.faces[jj].geometry.normal
        Vdir = -viewpoint.direction
        if (norm(n) * norm(Vdir)) == 0.0
            c = 0.0
        else
            c = dot(n, Vdir) / (norm(n) * norm(Vdir))
        end
        Aproj_i = c * geo.faces[jj].geometry.area
    end)

    return 1:Ntri |> m |> collect
end
