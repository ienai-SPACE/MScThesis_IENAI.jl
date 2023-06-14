function projection(geo::AbstractGeometry, viewpoint::Viewpoint, Ntri)

    m = Map(jj -> begin
        n = geo.faces[jj].normal
        Vdir = -viewpoint.direction
        c = dot(n, Vdir) / (norm(n) * norm(Vdir))
        Aproj_i = c * geo.faces[jj].area
    end)

    return 1:Ntri |> m |> collect
end