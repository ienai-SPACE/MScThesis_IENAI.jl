module SatelliteGeometryCalculations

using LinearAlgebra, StaticArrays, SatelliteToolbox, GeometryBasics, TickTock

SV3{T} = SVector{3,T}

include("geometry.jl")
include("EnvironmentalConstants.jl")
include("PermanentProperties.jl")
include("SurfaceProps.jl")
include("GasStreamProps.jl")   # include("EnvironmentalCalcs.jl")
include("CoefficientCalculations.jl") # ----- include("GSICalculations.jl")

include("Areas.jl") # 
# include("MTalgorithm.jl")
# include("Origins.jl")
# include("Convex.jl")
# include("NonConvex.jl")
include("GeometryInputs.jl")
include("Orbit&DateInputs.jl")
include("IlluminationConvex.jl")
include("IlluminationNonConvex.jl")
include("DragOrientationConvex.jl")

# include("VectorizeCoefficients.jl")

# include("SweepConvex.jl")
# include("SweepNonConvex.jl")

end

