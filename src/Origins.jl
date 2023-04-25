"""   function: Origins

GOAL: 
     - Define the coordinates of the origins of the rays
     - The origins are placed on a plane, perpendicular to the velocity vector
     - The source of origins has a circular shape

INPUTS:
     - Vrel
     - distance: distance from body center to perpendicular plane from where the rays iniciate (the plane should be outside of the satellite body)  
     - rmax: radius of the ray source on the perpendicular plane 
OUTPUTS:
      - O         :: coordinates on the perpendicular plane from where the rays iniciate
      - Norig     :: number of ray origins

      PENDING:
      *unify the 'Areas' module
      *verify against the spherical tri-mesh
"""

abstract type RaySampler end

struct GridFilter <: RaySampler
    ray_density::Float64 # rays / m²
end

struct MonteCarloSampler <: RaySampler
    ray_density::Float64 # rays / m²
end

struct FibonacciSampler <: RaySampler
    ray_density::Float64 # rays / m²
end

function generate_two_normals(u::SV3)
    u = normalize(u)
    v_trial = [1.0, 1.0, 1.0] |> normalize
    v_pre = cross(u, v_trial)
    if norm(v_pre) < 1e-3
        v_trial = [1.0, 1.0, -1.0] |> normalize
        v_pre = cross(u, v_trial)
    end
    v = v_pre |> normalize
    w = cross(u, v)
    (v, w)
end

function generate_ray_origins(gf::GridFilter, dir, rmax, distance)
    step = 1 / sqrt(gf.ray_density)
    x_grid = -rmax:step:rmax
    y_grid = -rmax:step:rmax
    ny = length(y_grid)
    nx = length(x_grid)
    x_coords = x_grid' .* ones(ny)
    y_coords = ones(nx)' .* y_grid
    size_x = size(x_coords)[1] * size(x_coords)[2]
    size_y = size(y_coords)[1] * size(y_coords)[2]
    x_c = reshape(x_coords, (size_x,))
    y_c = reshape(y_coords, (size_y,))
    points = collect(zip(x_c, y_c)) |> filter(p -> norm(p) <= rmax)
    center = dir * distance
    points_new = zeros(SVector{3,Float64}, length(points))
    v, w = generate_two_normals(dir)
    for (i, p) ∈ enumerate(points)
        points_new[i] = center + p[1] * v + p[2] * w
    end
    points_new
end

function generate_ray_origins_old(dir, rmax::T, distance) where {T}


    Dangle = 20           #Delta angle for the polar definition of origins

    #------pre-allocation-------------------
    xt = zeros(T, Int(360 / Dangle))
    yt = zeros(rmax, Int(360 / Dangle))
    zt = zeros(rmax, Int(360 / Dangle))
    Xg = zeros(rmax * Int(360 / Dangle), 3)
    #---------------------------------------

    λ = atan(dir[2] / dir[1])
    θ = asin(dir[3] / norm(dir))


    #rotation matrix: from the local ref. frame at the center of the satellite to the center of the source ref. frame
    R = [-sin(λ) cos(λ) 0; -sin(θ)*cos(λ) -sin(θ)*sin(λ) cos(θ); cos(θ)*cos(λ) cos(θ)*sin(λ) sin(θ)]

    x0_v = dir * distance

    #create origin coordinates on the plane's local ref. frame
    for gg ∈ 0:Int((360 / Dangle - 1))
        for rr ∈ 1:rmax
            xt[rr, gg+1] = rr * cos(gg * Dangle * π / 180)
            yt[rr, gg+1] = rr * sin(gg * Dangle * π / 180)
            zt[rr, gg+1] = 0
        end
    end

    #transform origin coordinates into the local satellite ref. frame    
    counter = 0
    for gg ∈ 0:Int((360 / Dangle - 1))
        for rr ∈ 1:rmax
            counter += 1
            Xg[counter, :] = transpose(R) * [xt[rr, gg+1]; yt[rr, gg+1]; zt[rr, gg+1]] + x0_v
        end
    end

    Norig = rmax * Int(360 / Dangle)
    return Xg, Norig
end