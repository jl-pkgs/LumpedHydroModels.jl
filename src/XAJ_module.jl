# https://github.com/OuyangWenyu/hydromodel/blob/11d9b7a0a6faa01e8a9d01dc7f296073e6a7f702/hydromodel/models/xaj.py

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
function Evapotranspiration!(P::T, PET::T, state::StateXAJ{T}; param::XAJ{T}) where {T<:Real}
  (; WL, WU) = state
  (; C, WLM) = param

  EU, EL, ED = T(0)
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
function Runoff(PE::T, WM, state::StateXAJ{T}; param::XAJ{T}) where {T<:Real}
  (; B, IM) = param
  (; W) = state

  WMM = WM * (1 .+ B) # 流域平均含水量W与最大储水量的关系, Eq. 2-54
  a = WMM * (1 - (1 - W / WM)^(1 / (1 + B))) # Eq. 2-58

  FR = 1 - (1 - a / WMM)^b # 2-55，产流面积α0
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


function Runoff_divide3S(state::StateXAJ{T}, FR1::T; param::XAJ{T}) where {T<:Real}
  (; SM, EX, KI, KG) = param
  (; PE, R, FR) = state # 在计算R时，FR已经进行了更新
  S1 = state.S # 上一时刻

  # FR = PE / R
  if FR == 0
    RS, RI, RG, S = 0.0
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

    RI = KI * S * FR
    RG = KG * S * FR
    S = S * (1 - KI - KG)
  end

  @pack! state = RS, RI, RG, S
end
