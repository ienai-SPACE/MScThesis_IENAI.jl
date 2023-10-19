module SatelliteGeometryCalculations

using LinearAlgebra, StaticArrays, SatelliteToolbox, TickTock, Accessors, Transducers, BasicInterpolators
using GeometryBasics: Point

SV3{T} = SVector{3,T}

include("viewpoint.jl")
include("geometry.jl")
include("EnvironmentalConstants.jl")
include("PermanentProperties.jl")
include("SurfaceProps.jl")
include("GasStreamProps.jl")   # include("EnvironmentalCalcs.jl")
include("CoefficientCalculations.jl") # ----- include("GSICalculations.jl")
include("GeometryInputs.jl")
include("Areas.jl") # 
# include("MTalgorithm.jl")
# include("Origins.jl")
# include("Convex.jl")
# include("NonConvex.jl")

include("Orbit&DateInputs.jl")

include("sweep_v2.jl")
include("projection.jl")

include("aerodynamics.jl")
include("solar.jl")

# include("VectorizeCoefficients.jl")


end

