using LinearAlgebra

"""
    unitaryDirections(velocity, normal)

Obtain the directions of the drag, lift, shear and pressure coefficients on a plate

#INPUTS
-`velocity`
-`normal`
#OUTPUTS
-`u_D`
-`u_L`
-`u_tau`
-`u_P`
"""

function unitaryDirections(velocity, normal)

    #Depending on the numbering direction of the triangle normals will go in one or another direction ------> TO BE CHECKED
    #------------------------------------------------------------

    u_D = -velocity |> normalize

    u_L_num = cross(cross(u_D, normal), u_D)

    u_L = -u_L_num / norm(u_L_num)

    u_tau_num = cross(cross(u_D, normal), normal)
    u_tau = -u_tau_num / norm(u_tau_num)

    u_P = -normal

    if norm(u_L_num) < 1e-5
        u_L = [0.0, 0.0, 0.0,]
    end
    if norm(u_tau_num) < 1e-5
        u_tau = [0.0, 0.0, 0.0,]
    end

    return u_D, u_L, u_P, u_tau
end

#TEST----------------------------------------------


# velocity = [1, 1, 0]
# normal = [1, -1, 0]

# u_D, u_L, u_tau, u_P = unitaryDirections(velocity, normal)

# normal = [1, -1, 0]
# u_d = [0.7071067811865475, 0.7071067811865475, 0]
# u_l_dum = cross(u_d, normal)
# u_l_num = cross(u_l_dum, u_d)
# u_L = -u_l_num / norm(u_l_num)