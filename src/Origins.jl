
# using LinearAlgebra, StaticArrays, SatelliteToolbox
# SV3{T} = SVector{3,T}

using Random

abstract type RaySampler end

"""
    Sampler <: RaySampler

Definition of the method to be used to populate the disk of ray origins

#Sampler options
-`GridFilter`
-`FibonacciSampler`
-`MonteCarlSampler`

#field
-`ray_density::Float64` : [rays / m²]
"""

struct GridFilter <: RaySampler
    ray_density::Float64 # rays / m²
end

struct MonteCarloSampler <: RaySampler
    ray_density::Float64 # rays / m²
end

struct FibonacciSampler <: RaySampler
    ray_density::Float64 # rays / m²  

end

"""
    generate_two_normals(u::SV3)

Generate an orthogonal plane defined by unit vectors `v` and `w` w.r.t. the input vector `u`

#INPUT:
-`u::SVector{3,T}`
#OUTPUT:
-`v::Vector`
-`w::Vector`
"""

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

"""
    generate_ray_origins(gf::GridFilter, dir, rmax, distance)

Create a uniformly distributed source of rays

#INPUTS:
-`gf::GridFilter` : density of rays [rays/m^2]
-`dir`
-`rmax`
-`distance`
#OUTPUT:
-`points_new::SVector{3,Float64}`
-`Norig::Int` : number of ray origins 
"""

function generate_ray_origins(gf::GridFilter, dir, rmax, distance)

    dir = -dir #same sense as velocity vector

    step = 1 / sqrt(gf.ray_density)
    x_grid = -rmax:step:rmax      #maybe rmax/2 is enough
    y_grid = -rmax:step:rmax
    ny = length(y_grid)
    nx = length(x_grid)
    x_coords = x_grid' .* ones(ny)
    y_coords = ones(nx)' .* y_grid
    size_x = size(x_coords)[1] * size(x_coords)[2]
    size_y = size(y_coords)[1] * size(y_coords)[2]
    x_c = reshape(x_coords, (size_x,))
    y_c = reshape(y_coords, (size_y,))
    points = filter(p -> norm(p) <= rmax, collect(zip(x_c, y_c)))
    #points = collect(zip(x_c, y_c)) |> filter(p -> norm(p) <= rmax)
    center = dir * distance
    points_new = zeros(SVector{3,Float64}, length(points))
    v, w = generate_two_normals(dir)
    for (i, p) ∈ enumerate(points)
        points_new[i] = center + p[1] * v + p[2] * w #`center` to displace the plane in the opposite direction to the satellite
    end
    points_new, Int(length(points))
end


"""
    generate_ray_origins(gf::FibonacciSampler, dir, rmax, distance)

Create a uniformly distributed source of rays using Fibonacci algorithm
Source: http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/

#INPUTS:
-`gf::FibonacciSampler` : density of rays [rays/m^2]
-`dir`
-`rmax`
-`distance`
#OUTPUT:
-`points_new::SVector{3,Float64}`
-`Norig::Int` : number of ray origins 
"""

function generate_ray_origins(gf::FibonacciSampler, dir, rmax, distance)

    dir = -dir #same sense as velocity vector

    goldenRatio = (1 + sqrt(5)) / 2
    n = ceil(gf.ray_density * (π * rmax^2)) #number of distributed points inside the disk of radius rmax

    x_v = LinRange(0, n, Int(n + 1)) / goldenRatio
    y_v = LinRange(0, n, Int(n + 1)) / n
    θ = 2 * π * x_v
    r = sqrt.(y_v) * rmax
    x_coords = cos.(θ) .* r
    y_coords = sin.(θ) .* r

    center = dir * distance
    points_new = zeros(SVector{3,Float64}, length(x_coords))
    v, w = generate_two_normals(dir)
    for ii ∈ 1:length(collect(x_coords))
        points_new[ii] = center + x_coords[ii] * v + y_coords[ii] * w #center to displace the plane in the opposite direction to the `dir`
    end

    points_new, Int(length(x_coords))
end


"""
    generate_ray_origins(gf::MonteCarloSampler, dir, rmax, distance)

Create a randomly distributed source of rays using Monte Carlo algorithm, based on test particle Monte Carlo method
Source: Xuhong Jina et al,"Monte Carlo simulation for aerodynamic coefficients of satellites in LowEarth Orbit"

#INPUTS:
-`gf::FibonacciSampler` : density of rays [rays/m^2]
-`dir`
-`rmax`
-`distance`
#OUTPUT:
-`points_new::SVector{3,Float64}`
-`Norig::Int` : number of ray origins 
"""

function generate_ray_origins(gf::MonteCarloSampler, dir, rmax, distance)

    dir = -dir #same sense as velocity vector

    n = ceil(gf.ray_density * (4 * rmax^2)) #number of randomly distributed points inside the disk of radius rmax
    # rng = Xoshiro(1234) # create a new Xoshiro random number generator object
    # R1 = rand(rng, Int(n))   # generate a random number using the rng object
    # R2 = rand(rng, Int(n))
    # R3 = rand(rng, Int(n))
    # R4 = rand(rng, Int(n))


    rng = Xoshiro() # create a new Xoshiro random number generator object
    R1 = rand(rng, Int(n))   # generate a random number using the rng object
    R2 = rand(rng, Int(n))

    Rc = rmax

    x_c = Rc * sqrt.(R1) .* cos.(2 * π * R2)
    y_c = Rc * sqrt.(R1) .* sin.(2 * π * R2)


    points = filter(p -> norm(p) <= rmax, collect(zip(x_c, y_c)))

    center = dir * distance
    points_new = zeros(SVector{3,Float64}, length(points))
    v, w = generate_two_normals(dir)
    for (i, p) ∈ enumerate(points)
        points_new[i] = center + p[1] * v + p[2] * w #center to displace the plane in the opposite direction to the `dir`

    end
    points_new, Int(length(points))
end


# export generate_two_normals

#TEST
#-------------------------------

# dir = SV3([1 / sqrt(2) 1 / sqrt(2) 0])
# rmax = 2
# distance = 10
# samplerG = GridFilter(10)
# samplerF = FibonacciSampler(10)
# samplerMC = MonteCarloSampler(1)
# PN, Norig = generate_ray_origins(samplerMC, dir, rmax, distance);







#=
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
=#

