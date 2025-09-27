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

  # routing::AbstractRouting{FT} = RoutingVoid{FT}()
  routing::MultiRouting{FT} = MultiRouting{FT}(
    RS=GammaUH{FT}(),                # 伽马单位线
    RI=LinearReservoirSM{FT}(), # 线性水库, 使用 CI 初始值
    RG=LinearReservoirGW{FT}()  # 线性水库, 使用 CG 初始值
  )
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

  isnothing(state) && (state = StateXAJ{T}())

  for t = 1:ntime
    _P = P[t]
    _PET = PET[t]
    FR1 = t == 1 ? 0.0 : FR[t-1]

    Evapotranspiration!(state, _P, _PET; param) # EU, EL, ED, ET, PE
    Runoff!(state; param)                        # R, R_IM, FR
    UpdateW!(state, _P; param)                   # WU, WL, WD, W

    Runoff_divide3S!(state, FR1; param)          # RS, RI, RG, S
    output[t] = state
  end

  Rsim = route_flow(param.routing, RS, RI, RG) # 假设 dt=1 个时间步
  output.Rsim = Rsim
  output |> DataFrame
end


"""
三层蒸发模型计算

参数:
- `P`   : 流域平均降水量 (mm)
- `PET` : 潜在蒸散发量 (mm)
- `WLM` : 下层土壤平均蓄水容量 (mm)
- `C`   : 深层蒸发系数 (无量纲)  
- `WU`  : 上层土壤初始含水量 (mm)
- `WL`  : 下层土壤初始含水量 (mm)
"""
function Evapotranspiration!(state::StateXAJ{T}, P::T, PET::T; param::XAJ{T}) where {T<:Real}
  (; WL, WU) = state
  (; K, C, WLM) = param
  PET = PET * K

  EU = EL = ED = T(0)
  # Eqs. 2-28，三层蒸发模式
  if WU + P >= PET
    EU = PET
  else
    EU = WU + P
    if WL >= C * WLM                      # 第二阶段，线性下降
      EL = C * (PET - EU)
    elseif C * (PET - EU) <= WL < C * WLM # 第三阶段，气体扩散, L2 有足够的水
      EL = C * (PET - EU)
    else                                  # 第三阶段，气体扩散, L2 没有足够的水
      EL = WL
      ED = C * (PET - EU) - EL # 气体扩散的剩余蒸发能力
    end
  end
  ET = EU + EL + ED
  PE = P - ET
  # 再进行一个限制，蒸发不能超过现在的水量
  @pack! state = EU, EL, ED, ET, PE
end


"""
降雨径流计算

- `B`  : 蓄水容量曲线的指数参数 (无量纲)
- `IM` : 不透水面积系数 (0-1)
- `WM` : 流域平均蓄水容量 = um+lm+dm (mm)
- `W`  : 流域初始土壤含水量 = wu0+wl0+wd0 (mm)
"""
function Runoff!(state::StateXAJ{T}; param::XAJ{T}) where {T<:Real}
  (; B, IM, WUM, WLM, WDM) = param
  (; W, PE) = state
  WM = WUM + WLM + WDM

  # 保护性代码，防止 (1 - W / WM) 为负
  # 使用 eps(T) 作为机器精度，或一个很小的数如 1e-10
  # 确保 W/WM 不超过 1.0
  ratio_WM = W / WM
  if ratio_WM > 1.0
    @warn "W > WM in Runoff! (W=$(W), WM=$(WM)), setting ratio_WM = 1.0 - eps(T). This might indicate an issue in UpdateW! or initial conditions."
    ratio_WM = 1.0 - eps(T)
  elseif ratio_WM < 0.0
    @warn "W < 0 in Runoff! (W=$(W)), setting ratio_WM = 0.0. This might indicate an issue in UpdateW! or initial conditions."
    ratio_WM = 0.0
  end

  WMM = WM * (1 + B) # 流域平均含水量W与最大储水量的关系, Eq. 2-54
  a = WMM * (1 - (1 - ratio_WM)^(1 / (1 + B))) # Eq. 2-58

  FR = 1 - (1 - a / WMM)^B # 2-55，产流面积α0
  R = T(0)
  if PE > 0
    if PE + a < WMM
      R = PE + W - WM + WM * (1 - (PE + a) / WMM)^(B + 1) # Eq. 2-60
    else
      R = PE + W - WM # Eq. 2-61
    end
    # FR = R / PE            # TODO, 2-59, 检查哪一个更可靠
  end

  R = max(R * (1 - IM), 0) # 自然产量
  R_IM = max(PE * IM, 0)   # 不透水界面全部产流

  @pack! state = R, R_IM, FR
end


function Runoff_divide3S!(state::StateXAJ{T}, FR1::T; param::XAJ{T}) where {T<:Real}
  (; SM, EX, KI, KG) = param
  (; PE, R, R_IM, FR) = state # 在计算R时，FR已经进行了更新
  S1 = state.S # 上一时刻

  # FR = PE / R
  RS = RI = RG = T(0.0)
  if FR == 0
    S = T(0.0)
  else
    SMM = SM * (1 + EX)
    S = S1 * FR1 / FR # 产流面积上的平均自由水蓄量

    AU = SMM * (1 - (1 - S / SM)^(1 / (EX + 1)))
    if PE > 0
      if PE + AU < SMM
        RS = FR * (PE + SM * (1 - (PE + AU) / SMM)^(EX + 1) - SM + S)
        S = S + (R - RS) / FR
      else
        RS = FR * (PE + S - SM) # 多出来的水都转化为RS
        S = SM
      end
    end
    RS = RS + R_IM # 不透水界面加到RS
    RI = KI * S * FR
    RG = KG * S * FR
    S = S * (1 - KI - KG)
  end

  @pack! state = RS, RI, RG, S
end
