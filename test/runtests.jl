using LumpedHydroModels, Test

indir = "$(@__DIR__)/.." |> abspath

include("test-linear_reservoir.jl")
include("test-param.jl")
include("test-run_XAJ.jl")
# include("test-xaj_routing.jl")
