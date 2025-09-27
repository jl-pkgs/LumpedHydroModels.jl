export linear_Cs2K, linear_K2Cs
export linear_reservoir_uh
export uh_gamma, uh_linear
export conv_uh, guess_uh
export AbstractSingleRouting, MultiRouting, route_flow

export LinearReservoir, LinearReservoirSM, LinearReservoirGW

using SpecialFunctions: gamma

include("linear_reservoir.jl")
include("unit_hydrograph.jl")


abstract type AbstractSingleRouting{FT} <: AbstractRouting{FT} end

abstract type LinearReservoir{FT} <: AbstractSingleRouting{FT} end


@bounds @units @with_kw mutable struct MultiRouting{FT} <: AbstractRouting{FT}
  RS::AbstractSingleRouting{FT} = RoutingVoid{FT}()
  RI::AbstractSingleRouting{FT} = RoutingVoid{FT}()
  RG::AbstractSingleRouting{FT} = RoutingVoid{FT}()
end

@bounds @units @with_kw mutable struct GammaUH{FT} <: AbstractSingleRouting{FT}
  α::FT = 2.5 | (1.0, 5.0) | "-"      # 伽马分布形状参数
  θ::FT = 1.0 | (0.1, 3.0) | "hr"     # 伽马分布尺度参数 (时间)
  τ_max::Int = 10 | (5, 50) | "t"     # 单位线长度 (时间步)
  UH::Vector{FT} = uh_gamma.(1:τ_max; α, θ) |> normalize_uh
end

function init_routing(routing::GammaUH{FT}) where {FT<:Real}
  (; α, θ, τ_max) = routing
  routing.UH .= uh_gamma.(1:τ_max; α, θ) |> normalize_uh
end

# CS::FT = 0.7 | (0.01, 0.9) | "-"      # 汇流参数 
# CI::FT = 0.7 | (0.5, 0.99) | "-"      # 壤中流，线性水库汇流参数
# CG::FT = 0.98 | (0.95, 0.998) | "-"   # 地下径流，线性水库汇流参数

@bounds @units @with_kw mutable struct LinearReservoirSM{FT} <: LinearReservoir{FT}
  CS::FT = 0.7 | (0.5, 0.99) | "-"    # 线性水库消退系数
  nseg::Int = 1 | (1, 10) | "-"       # 水库串联数
  # K::FT = 0.5 | (0.1, 5.0) | "hr"   # 蓄水常数 (可选，如果要用 K 定义)
  # dt::FT = 1.0 | (0.1, 5.0) | "hr"  # 时间步长 (可选)
end

@bounds @units @with_kw mutable struct LinearReservoirGW{FT} <: LinearReservoir{FT}
  CS::FT = 0.98 | (0.95, 0.998) | "-"  # 线性水库消退系数
  nseg::Int = 1 | (1, 10) | "-"        # 水库串联数
end


"""
对单个径流成分进行汇流计算
"""
function route_flow(routing::RoutingVoid{FT}, I::AbstractVector{FT}; ignored...) where {FT<:Real}
  return I # RoutingVoid 不改变输入
end


function route_flow(routing::GammaUH{FT}, I::AbstractVector{FT}; ignored...) where {FT<:Real}
  return conv_uh(I, routing.UH; start=1)
end


function route_flow(routing::LinearReservoir{FT}, I::AbstractVector{FT}; dt::FT=FT(1)) where {FT<:Real}
  # CS_calc = linear_K2Cs(routing.K; dt) # 如果启用 K 参数
  (; CS, nseg) = routing
  return linear_reservoir(I; CS, nseg, dt)
end


function route_flow(routing::MultiRouting{FT}, rs::V, ri::V, rg::V;
  dt::FT=FT(1)) where {FT<:Real,V<:AbstractVector{FT}}

  Qrs = route_flow(routing.RS, rs; dt)
  Qri = route_flow(routing.RI, ri; dt)
  Qrg = route_flow(routing.RG, rg; dt)
  return Qrs .+ Qri .+ Qrg
end
