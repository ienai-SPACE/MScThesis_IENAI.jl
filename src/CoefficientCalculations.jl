include("GSICalculations.jl")
include("VectorizeCoefficients.jl")


"""
    ElementInteractionProps{T}

- `δ::T`    : angle between oncoming direction vector and normal to the surface
- `SRF::T`  : molecular mass of the surface atom
- `Tw::T`   : temperature at the wall
"""

struct ElementInteractionProps{T}
    δ::T
    SRF::T
    Tw::T
end

ElementInteractionProps(surfprops::SurfaceProps, srfMat, angle) = ElementInteractionProps(angle, srfMat, surfprops.Tw)


abstract type InteractionGeometry end

"""
    InteractionGeometryHomo{T}

- `area::T`       : [m^2]
- `angle::T`      : [rad]
"""

struct InteractionGeometryHomo{T} <: InteractionGeometry
    area::T
    angle::T
end

"""
    InteractionGeometryHetero{T}

- `area::T`       : [m^2]
- `angle::T`      : [rad]
- `index::T`      : global index relating material property
"""
struct InteractionGeometryHetero{T} <: InteractionGeometry
    area::T
    angle::T
    index::T
end


"""
    AerodynamicCoefficients{T}

- `Cd::Float64`
- `Cl::Float64`
- `Cp::Float64`
- `Ctau::Float64`
"""
struct AerodynamicCoefficients{T}
    Cd::T
    Cl::T
    Cp::T
    Ctau::T
end

"""
    CoefficientsVectorized{T}

- `CD_norm::T`
- `CD_dir::Vector{T}`
- `CL_norm::T`
- `CL_dir::Vector{T}`
- `CP_norm::T`
- `CP_dir::Vector{T}`
- `CTau_norm::T`
- `CTau_dir::Vector{T}`

"""

struct CoefficientsVectorized{T}
    CD_norm::T
    CD_dir::Vector{T}
    CL_norm::T
    CL_dir::Vector{T}
    CP_norm::T
    CP_dir::Vector{T}
    CTau_norm::T
    CTau_dir::Vector{T}

    function CoefficientsVectorized(coeffs_vec::Vector{AbstractVector{Float64}})
        CD_norm = norm(coeffs_vec[1])
        CL_norm = norm(coeffs_vec[2])
        CP_norm = norm(coeffs_vec[3])
        CTau_norm = norm(coeffs_vec[4])

        if CD_norm == 0
            CD_dir = [0.0, 0.0, 0.0]
        else
            CD_dir = coeffs_vec[1] / norm(coeffs_vec[1])
        end
        if CL_norm == 0
            CL_dir = [0.0, 0.0, 0.0]
        else
            CL_dir = coeffs_vec[2] / norm(coeffs_vec[2])
        end
        if CP_norm == 0
            CP_dir = [0.0, 0.0, 0.0]
        else
            CP_dir = coeffs_vec[3] / norm(coeffs_vec[3])
        end
        if CTau_norm == 0
            CTau_dir = [0.0, 0.0, 0.0]
        else
            CTau_dir = coeffs_vec[4] / norm(coeffs_vec[4])
        end

        return CD_norm, CD_dir, CL_norm, CL_dir, CP_norm, CP_dir, CTau_norm, CTau_dir
    end
end

Base.:*(a::Real, c::AerodynamicCoefficients) = AerodynamicCoefficients(a * c.Cd, a * c.Cl, a * c.Cp, a * c.Ctau)
Base.:*(c::AerodynamicCoefficients, a::Real) = a * c
Base.:+(a::AerodynamicCoefficients, b::AerodynamicCoefficients) = AerodynamicCoefficients(a.Cd + b.Cd, a.Cl + b.Cl, a.Cp + b.Cp, a.Ctau + b.Ctau)
Base.:/(a::AerodynamicCoefficients, d::Real) = AerodynamicCoefficients(a.Cd / d, a.Cl / d, a.Cp / d, a.Ctau / d)

"""
    compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::InteractionGeometryHomo, Vrel_v)

Compute the toal drag, lift, pressure, and shear coefficients evaluated for a single area element

#INPUT:
- `surfprops::SurfaceProps`
- `gasprops::GasStreamProperties`
- `intgeo::InteractionGeometryHomo`
- `Vrel_v:: Vector{Float64}`     [m/s]
#OUTPUT:
- `AerodynamicCoefficients(Cd, Cl, Cp, Ctau)`
"""

function compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::InteractionGeometryHomo, Vrel_v, normals)
    element_interaction = ElementInteractionProps(surfprops, surfprops.m_srf, intgeo.angle)
    Cd, Cl, Cp, Ctau = DRIA_GSI(element_interaction, gasprops, Vrel_v, normals)
    AerodynamicCoefficients(Cd, Cl, Cp, Ctau)
end

"""
    compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::InteractionGeometryHetero, Vrel_v)

Compute the toal drag, lift, pressure, and shear coefficients evaluated for a single area element

#INPUT:
- `surfprops::SurfaceProps`
- `gasprops::GasStreamProperties`
- `intgeo::InteractionGeometryHetero`
- `Vrel_v:: Vector{Float64}`     [m/s]
#OUTPUT:
- `AerodynamicCoefficients(Cd, Cl, Cp, Ctau)`
"""

function compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::InteractionGeometryHetero, Vrel_v, normals)
    element_interaction = ElementInteractionProps(surfprops, intgeo.index, intgeo.angle)
    Cd, Cl, Cp, Ctau = DRIA_GSI(element_interaction, gasprops, Vrel_v, normals)
    AerodynamicCoefficients(Cd, Cl, Cp, Ctau)
end

""" 
    compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::Vector{<:InteractionGeometry}, Vrel_v)

Compute the toal drag, lift, pressure, and shear coefficients evaluated for all areas of the element

#INPUT:
- `surfprops::SurfaceProps`
- `gasprops::GasStreamProperties`
- `intgeo::Vector{<:InteractionGeometry}`
- `Vrel_v::Vector{Float64}`
#OUTPUT:
- `AerodynamicCoefficients{Float64}` : total Cd, Cl, Cp, Ctau
- `Atot::Float64`                    : [m^2]
- `Aproj:.Float64`                   : [m^2]
"""

function compute_coefficients(surfprops::SurfaceProps, gasprops::GasStreamProperties, intgeo::Vector{<:InteractionGeometry}, Vrel_v, normals)
    aero_coeffs = map(intgeo) do intgeo
        compute_coefficients(surfprops, gasprops, intgeo, Vrel_v, normals)
    end
    areas = map(intgeo) do x
        x.area
    end
    angles = map(intgeo) do x
        x.angle
    end
    # num = sum(aero_coeffs .* areas)
    num = aero_coeffs .* areas

    # print(normals, //)

    coeffs_vec = vectorizeCoeffs(num, normals, Vrel_v) #CD,CL,Cp, Ctau

    # print(coeffs_vec, //)
    # print(areas, //)
    # print(cos.(angles), //)
    # print(rad2deg.(angles), //)
    # return scaled_coefficients = sum(C*A_i)/Aref, A_tot
    CoefficientsVectorized(coeffs_vec ./ sum(areas .* abs.(cos.(angles)))), sum(areas), sum(areas .* abs.(cos.(angles)))
    # coeffs_vec ./ sum(areas .* cos.(angles)), sum(areas), sum(areas .* cos.(angles))
end

export InteractionGeometryHomo, InteractionGeometryHetero, compute_coefficients