module LumpedHydroModels

export Params
using Parameters
using SpecialFunctions: gamma
import DataFrames: DataFrame

# 2.5x faster power method
"Faster method for exponentiation"
@fastmath pow(x::Real, y::Real) = exp(y * log(x))

include("Parameters.jl")
include("DataType.jl")
include("Routing/linear_reservoir.jl")
include("Routing/uh_gamma.jl")

include("ultilize.jl")
include("XAJ_model.jl")
include("XAJ_module.jl")
include("XAJ_UpdateW.jl")

export XAJ

end # module LumpedHydroModels
