push!(LOAD_PATH, "../src/")
using Documenter
using Revise
using SatelliteGeometryCalculations

ENV["GKSwstype"] = "100"
DocMeta.setdocmeta!(SatelliteGeometryCalculations, :DocTestSetup, :(using SatelliteGeometryCalculations); recursive=true)
Documenter.makedocs(root=joinpath(@__DIR__),
        source="src",
        build="build",
        clean=true,
        doctest=true,
        modules=Module[SatelliteGeometryCalculations],
        repo="",
        highlightsig=true,
        sitename="SatelliteGeometryCalculations",
        expandfirst=[],
        pages=[
                "Home" => "intro.md",
                "Reference" => [
                        "reference/types.md",
                        "reference/functions.md"
                ],
        ],
        format=Documenter.HTML(
                ansicolor=true,
                prettyurls=get(ENV, "CI", nothing) == "true",
                sidebar_sitename=false
        )
)
