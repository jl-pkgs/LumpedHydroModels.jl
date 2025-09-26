module LumpedHydroModels

export Params
using Parameters
using SpecialFunctions: gamma
using DataFrames: DataFrame

include("Parameters.jl")
include("DataType.jl")
include("Routing/linear_reservoir.jl")
include("Routing/uh_gamma.jl")

include("XAJ_model.jl")
include("XAJ_module.jl")
include("XAJ_UpdateW.jl")

export XAJ

end # module LumpedHydroModels
