"""
    projection(geo::AbstractGeometry, viewpoint::Viewpoint, Ntri)

Projection of areas according to a viewpoint direction

#INPUTS
- `geo::AbstractGeometry`
- `viewpoint::Viewpoint`
- `Ntri`
#OUTPUTS
- vector of projected areas
"""

function projection(geo::AbstractGeometry, viewpoint::Viewpoint, Ntri)

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