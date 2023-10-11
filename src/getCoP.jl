"""
    getCoP(PCSA, intercept_info, barycenters)

Calculation of the center of pressure wrt the local ref.frame at which the vertices of the geometry are defined

# Inputs
- `PCSA` : projected cross-sectional area = Aproj
- `intercept_info::Vector{InteractionGeometryHomo{Float64}}` or `Vector{InteractionGeometry{Float64}}`
- `barycenters::Vector{Vectors}`

# Outputs
- `CoP::Vector`
"""
# make sure the areas and barycenters correspond to the same elements!!!!!
function getCoP(PCSA, intercept_info, barycenters)
    # for ii in 1:lastindex(barycenters)
    #     CoP[ii] = (1 / PCSA) * intercept_info[ii].area * barycenters[ii]
    # end
    println("length(barycenters) = ", length(barycenters))
    println("length(intercept_info) = ", length(intercept_info))
    CoP = (1 / PCSA) * sum([intercept_info[ii].area * cos(intercept_info[ii].angle) for ii in 1:lastindex(barycenters)] .* barycenters, dims=1)
    SVector(CoP[1][1], CoP[1][2], CoP[1][3])
end