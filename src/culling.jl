"""
    culling(MeshVertices, Vdir)

Filter out all non-forward facing faces_hit_idx_nonunique

# Input
- `MeshVertices::Matrix`   : coordinates of each vertex in the format [x1y1z1x2y2z2x3y3z3;...]
- `Vdir::Vector`           : vector with the direction to be analyzed
# Output
- `triangles::Transpose{Float64,Matrix{Float64}}`       :triangles that pass the culling filter
"""
function culling(MeshVertices, Vdir)
    #Number of triangles
    Ntri = size(MeshVertices, 1)

    #pre-allocation
    counter = 0.0
    triangles = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    for jj ∈ 1:Ntri #for loop to iterate over all triangles/quads
        vertices = MeshVertices[jj, 1:9]

        #definition of the triangle vertices
        V1 = [vertices[1], vertices[2], vertices[3]]
        V2 = [vertices[4], vertices[5], vertices[6]]
        V3 = [vertices[7], vertices[8], vertices[9]]


        #indices should be numbered in c.c.w direction
        edge1 = V2 - V1
        edge2 = V3 - V1

        n = cross(edge1, edge2)

        dotProd = dot(Vdir, n)


        if dotProd > 1e-5   #theoretically > 0, but 1e-5 to account for numerical errors

            counter += 1
            if counter == 1
                triangles = [vertices[1], vertices[2], vertices[3], vertices[4], vertices[5], vertices[6], vertices[7], vertices[8], vertices[9]]
            else
                triangles = hcat(triangles, [vertices[1], vertices[2], vertices[3], vertices[4], vertices[5], vertices[6], vertices[7], vertices[8], vertices[9]])
            end
        end
    end

    return transpose(triangles)

end