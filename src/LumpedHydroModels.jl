module LumpedHydroModels

export Params
using Parameters
using SpecialFunctions: gamma
import FieldMetadata: @bounds, bounds, @units, units
import DataFrames: DataFrame

# 2.5x faster power method
"Faster method for exponentiation"
@fastmath pow(x::Real, y::Real) = exp(y * log(x))

include("Parameters.jl")
include("DataType.jl")
include("macro.jl")
include("Routing/Routing.jl")

include("ultilize.jl")
include("XAJ_model.jl")
include("XAJ_UpdateW.jl")

export XAJ
export GammaUH, LinearReservoirRouting, MultiRouting

end # module LumpedHydroModels
