using LinearAlgebra


""" function: areasConvex 

GOAL:
    - Obtain the area of a triangle/quad if the angle between the normal and the direction of assessment is 0 < θ < π/2
INPUT:
    - vertices :: vector with the coordinates of each vertex in the format [x1y1z1x2y2z2x3y3z3...]
    - dir      :: vector with the direction to be analyzed
OUTPUT:
    - area
    """


function areasConvex(vertices, dir)

    if length(vertices) == 9         #triangular mesh

        #definition of the triangle vertices
        V1 = [vertices[1], vertices[2], vertices[3]]
        V2 = [vertices[4], vertices[5], vertices[6]]
        V3 = [vertices[7], vertices[8], vertices[9]]


        #indices should be numbered in c.c.w direction
        edge1 = V2 - V1
        edge2 = V3 - V1
        n = cross(edge1, edge2)

        #θ = acos(dot(dir, n) / (norm(dir)*norm(n)))

        dotProd = dot(dir, n)

        if dotProd > 0
            area = norm(n) / 2        #area of the triangle
        else
            area = 0.0
        end

    elseif length(vertices) == 12     #quads

        #definition of the quad vertices
        V1 = [vertices[1], vertices[2], vertices[3]]
        V2 = [vertices[4], vertices[5], vertices[6]]
        V3 = [vertices[7], vertices[8], vertices[9]]
        V4 = [vertices[10], vertices[11], vertices[12]]



        #indices should be numbered in c.c.w direction
        edge1 = V2 - V1
        edge2 = V4 - V1
        n = cross(edge1, edge2)

        #θ = acos(dot(dir, n) / (norm(dir) * norm(n)))
        dotProd = dot(dir, n)

        if dotProd > 0
            edge3 = V4 - V3
            edge4 = V3 - V2
            n2 = cross(edge3, edge4)

            area = (norm(n) + norm(n2)) / 2                  #area of the quad
        else
            area = 0.0
        end

    end



    #Calculate the angle between triangle's normal and velocity direction
    crossProd = n
    dot_prod = dot(crossProd, dir)
    γ_dir = acos(dot_prod / (norm(crossProd) * norm(dir)))    #[radians]



    return [area, γ_dir]
end


#TESTING
#------------------------------------------------------
#=
vertices = [0 0 0 0 1 1 0 -1 1] #x1y1z1x2y2z2x3y3z3
dir = [1 1 0]
#results: [1.0, 0.7853981633974484 = pi/4]

#vertices = [0 0 0 0 1 0 0 1 1 0 0 1] 
#dir = [1 1 0]
#results: [1.0, 0.7853981633974484 = pi/4]


area, delta = areasConvex(vertices, dir)
=#