export linear_Cs2K, linear_K2Cs
export linear_reservoir_uh
export uh_gamma, uh_linear
export conv_uh, guess_uh


using SpecialFunctions: gamma

include("linear_reservoir.jl")
include("unit_hydrograph.jl")
