export XAJ, StateXAJ, OutputsXAJ, run_XAJ

# Table 2-11
@bounds @units @with_kw mutable struct XAJ{FT} <: AbstractHydroModel{FT}
  K::FT = 0.95 | (0.2, 1.5) | "-"       # 蒸发折算系数
  C::FT = 0.14 | (0.05, 0.20) | "-"     # 深层蒸发系数

  IM::FT = 0.1 | (0.01, 0.3) | "-"      # 不透水面积比例

  WUM::FT = 15.0 | (5.0, 20.0) | "mm"   # 上层土壤平均蓄水容量
  WLM::FT = 85.0 | (10.0, 90.0) | "mm"  # 下层
  WDM::FT = 20.0 | (10.0, 120.0) | "mm" # 深层

  SM::FT = 5.0 | (10.0, 60.0) | "mm"    # 自由水蓄水容量

  B::FT = 0.3 | (0.1, 0.6) | "-"        # 蓄水容量曲线的指数参数
  EX::FT = 0.7 | (0.5, 2.0) | "-"       # 自由水蓄水容量曲线指数

  KI::FT = 0.3 | (0.01, 0.7) | "-"      # 壤中流出流系数
  KG::FT = 0.3 | (0.01, 0.7) | "-"      # 地下水出流系数

  CS::FT = 0.7 | (0.01, 0.9) | "-"      # 汇流参数
  CI::FT = 0.7 | (0.5, 0.99) | "-"      # 壤中流，线性水库汇流参数
  CG::FT = 0.98 | (0.95, 0.998) | "-"   # 地下径流，线性水库汇流参数

  routing::AbstractRouting{FT} = RoutingVoid{FT}()
end
# XAJ() = XAJ{Float64}()

@with_kw mutable struct StateXAJ{FT} <: AbstractHydroState{FT}
  ET::FT = 0.0
  EU::FT = 0.0
  EL::FT = 0.0
  ED::FT = 0.0
  PE::FT = 0.0    # 净雨, P - ET

  WU::FT = 0.0
  WL::FT = 2.2
  WD::FT = 20.0
  W::FT = WU + WL + WD

  FR::FT = 0.0
  R::FT = 0.0     # 透水界面径流
  R_IM::FT = 0.0  # 不透水界面径流
  RS::FT = 0.0
  RI::FT = 0.0
  RG::FT = 0.0
  S::FT = 0.0     # 自由水蓄量
end
@DefOutputs StateXAJ Rsim


function run_XAJ(P::Vector{T}, PET::Vector{T}; param::XAJ{T}, state=nothing) where {T<:Real}
  ntime = length(P)
  output = OutputsXAJ{T}(; ntime)
  (; RS, RI, RG, FR) = output

  isnothing(state) && (state = StateXAJ{FT}())

  for t = 1:ntime
    _P = P[t]
    _PET = PET[t]
    FR1 = t == 1 ? 0.0 : FR[t-1]

    Evapotranspiration!(state, _P, _PET; param) # EU, EL, ED, ET, PE
    Runoff(state; param)                        # R, R_IM, FR
    UpdateW(state, _P; param)                   # WU, WL, WD, W

    Runoff_divide3S(state, FR1; param)          # RS, RI, RG, S
    output[t] = state
  end
  output.Rsim = RS + RI + RG
  output |> DataFrame
end
