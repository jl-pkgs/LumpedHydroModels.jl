abstract type AbstractHydroModel{FT} <: AbstractModel{FT} end
abstract type AbstractHydroState{FT} <: AbstractModel{FT} end

@bounds @units @with_kw mutable struct XAJ{FT} <: AbstractHydroModel{FT}
  K::FT = 0.0 | (0.1, 1.5) | "-"      # 蒸发折算系数
  C::FT = 0.0 | (0.9, 1.1) | "-"      # 深层蒸发系数

  IM::FT = 0.0 | (0.01, 0.1) | "-"    # 不透水面积比例

  UM::FT = 0.0 | (5.0, 20.0) | "mm"   # 上层土壤平均蓄水容量
  LM::FT = 0.0 | (60.0, 90.0) | "mm"  # 下层
  DM::FT = 0.0 | (60.0, 120.0) | "mm" # 深层

  # WM::FT = 0.0     # 水蓄量
  SM::FT = 0.0 | (10.0, 25.0) | "mm"  # 自由水蓄水容量

  B::FT = 0.0 | (0.1, 0.4) | "-"      # 蓄水容量曲线的指数参数
  EX::FT = 0.0 | (1.0, 1.5) | "-"     # 自由水蓄水容量曲线指数

  KI::FT = 0.0 | (0.01, 0.7) | "-"    # 壤中流出流系数
  KG::FT = 0.0 | (0.01, 0.7) | "-"    # 地下水出流系数

  CS::FT = 0.0 | (0.01, 0.9) | "-"    # 汇流参数
  CI::FT = 0.0 | (0.5, 0.99) | "-"    # 壤中流，线性水库汇流参数
  CG::FT = 0.0 | (0.95, 0.998) | "-"  # 地下径流，线性水库汇流参数
end


@with_kw mutable struct StateXAJ{FT} <: AbstractHydroState{FT}
  ET::FT = 0.0
  EU::FT = 0.0
  EL::FT = 0.0
  ED::FT = 0.0

  PE::FT = 0.0    # 净雨, P - ET

  W::FT = 0.0
  WU::FT = 0.0
  WL::FT = 0.0
  WD::FT = 0.0

  R::FT = 0.0     # 透水界面径流
  R_IM::FT = 0.0  # 不透水界面径流
  RS::FT = 0.0
  RI::FT = 0.0
  RG::FT = 0.0
  S::FT = 0.0     # 自由水蓄量
end


function XAJ_main(P::Vector{T}, PET::Vector{T}, params::XAJ{T}) where {T<:Real}
  (; K) = params

  ntime = length(P)
  S = zeros(T, ntime)
  FR = zeros(T, ntime)

  RS_IM = zeros(T, ntime)
  RS = zeros(T, ntime)
  RI = zeros(T, ntime)
  RG = zeros(T, ntime)
  ET = zeros(T, ntime)

  state = StateXAJ{FT}()

  for t = 2:ntime
    S1 = S[t-1]
    FR1 = FR[t-1]
    FR = FR[t]

    _P = P[t]
    _PET = PET[t] * K

    EU, EL, ED = calc_evap(_P, _PET, WU, WL; WLM, C)
    ET[t] = EU + EL + ED
    PE = _P - _PET

    R, RS_IM[t] = calc_runoff(PE, W; B, IM, WM)
    WU, WL, WD, W = Storage_Update_Zhong(P, R, state; params)
    # RS[t], RI[t], RG[t], S[t] = flow_divide(PE, R, S1, FR1, FR; params)
  end

  # RS = Rs + RS_IM
  ## 何时做汇流？
  Rs, Ri, Rg # area ==> Q
end
