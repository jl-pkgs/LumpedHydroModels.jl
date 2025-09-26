using LumpedHydroModels, Test

indir = "$(@__DIR__)/.." |> abspath

include("test-param.jl")
include("test-run_XAJ.jl")
