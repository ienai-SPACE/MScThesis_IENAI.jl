module SatelliteGeometryCalculations

using LinearAlgebra, StaticArrays, SatelliteToolbox, TickTock, Accessors, Transducers
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
# include("IlluminationConvex.jl")        #NOT NEEDED: substituted by "areasSphericalInput.jl"
# include("IlluminationNonConvex.jl")     #NOT NEEDED: substituted by "areasSphericalInput.jl"
# include("DragOrientationConvex.jl")     #NOT NEEDED: substituted by "areasSphericalInput.jl"

include("areasSphericalInput.jl")
include("sweep.jl")

# include("VectorizeCoefficients.jl")

# include("SweepConvex.jl")
# include("SweepNonConvex.jl")

end

