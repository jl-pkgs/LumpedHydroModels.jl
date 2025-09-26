module LumpedHydroModels

using ModelParams, UnPack
using Parameters

using SpecialFunctions: gamma
using DataFrames: DataFrame

include("Routing/linear_reservoir.jl")
include("Routing/uh_gamma.jl")

include("XAJ_model.jl")
include("XAJ_module.jl")
include("UpdateW.jl")

end # module LumpedHydroModels
