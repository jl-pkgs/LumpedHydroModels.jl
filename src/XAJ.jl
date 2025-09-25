"""
https://github.com/OuyangWenyu/hydromodel/blob/11d9b7a0a6faa01e8a9d01dc7f296073e6a7f702/hydromodel/models/xaj.py

XinAnJiang水文模型Julia实现
作者: kongdd (从Python翻译)
"""

using SpecialFunctions: gamma
using DSP: conv

abstract type AbstractHydroState{FT} <: AbstractModel{FT} end

@with_kw mutable struct StateXAJ{FT} <: AbstractHydroState{FT}
  Eu::FT = 0.0
  El::FT = 0.0
  Ed::FT = 0.0
  WU::FT = 0.0
  WL::FT = 0.0
  WD::FT = 0.0
  W::FT = 0.0
  RS::FT = 0.0
  RI::FT = 0.0
  RG::FT = 0.0
  S::FT = 0.0
end



"""
降雨径流计算

- `B`  : 蓄水容量曲线的指数参数 (无量纲)
- `IM` : 不透水面积系数 (0-1)
- `WM` : 流域平均蓄水容量 = um+lm+dm (mm)
- `W`  : 流域初始土壤含水量 = wu0+wl0+wd0 (mm)
"""
function calc_runoff(PE::T, W::T; B::T, IM::T, WM::T) where {T<:Real,V<:Vector{T}}
  WMM = WM * (1 .+ B) # 流域平均含水量W与最大储水量的关系, Eq. 2-54
  a = WMM * (1 - (1 - W / WM)^(1 / (1 + B))) # Eq. 2-58

  R = T(0)
  if PE > 0
    if PE + a < WMM
      R = PE + W - WM + WM * (1 - (PE + a) / WMM)^(B + 1) # Eq. 2-60
    else
      R = PE + W - WM # Eq. 2-61
    end
  end
  R = max(R * (1 - IM), 0) # 自然产量
  R_im = max(PE * IM, 0) # 不透水界面全部产流
  return R, R_im
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
function calc_evap(P::T, Ep::T, WU::T, WL::T; WLM::T, C::T) where {T<:Real}
  Eu, El, Ed = T(0)
  # Eqs. 2-28，三层蒸发模式
  if WU + P >= Ep
    Eu = Ep
  else
    Eu = WU + P
    if WL >= C * WLM                     # 第二阶段，线性下降
      El = C * (Ep - Eu)
    elseif C * (Ep - Eu) <= WL < C * WLM # 第三阶段，气体扩散
      El = C * (Ep - Eu)      # L2 有足够的水蒸发
    else
      El = WL
      Ed = C * (Ep - Eu) - El # 气体扩散的剩余蒸发能力
    end
  end
  # 再进行一个限制，蒸发不能超过现在的水量
  return Eu, El, Ed
end


# 土壤蓄水量更新, TODO 存在另一种写法
function Storage_Update_Zhong(WU::V, WL::V, WD::V,
  P::V, R::V, Eu::V, El::V, Ed::V;
  WUM::T, WLM::T, WDM::T) where {T<:Real,V<:Vector{T}}

  # 蒸发过程已经考虑了W的限制，因此这里主要处理的是水多的情景
  if WU + P - Eu - R <= WUM
    WU += P - Eu - R
    WL -= El
    WD -= Ed
    exceed_UL = 0.0
  else
    WU = WUM
    exceed_UL = WU + P - Eu - R - WUM # 多余出来的水

    if WL - El + exceed_UL <= WLM
      WL = WL - El + exceed_UL
      WD -= Ed
      exceed_LD = 0.0
    else
      WL = WLM
      exceed_LD = WL + exceed_UL - EL - WLM

      if WD - Ed + exceed_LD <= WDM
        WD = WD - Ed + exceed_LD
      else
        WD = WDM # TODO: 多出来的水去哪了？放到地下？
      end
    end
  end
  W = WU + WL + WD
  return WU, WL, WD, W
end



# 径流产生
"""
- k: 蒸散发折算系数，潜在蒸散发与参考作物蒸发的比值 (无量纲)
- b: 蓄水容量曲线指数参数 (无量纲)
- im: 不透水面积系数 (0-1)
- um: 上层土壤平均蓄水容量 (mm)
- lm: 下层土壤平均蓄水容量 (mm)
- dm: 深层土壤平均蓄水容量 (mm)  
- c: 深层蒸发系数 (无量纲)
"""
function generation(P::V, PET::V;
  K::T, B::T, IM::T, WUM::T, WLM::T, WDM::T, C::T,
  wu0::Union{V,Nothing}=nothing,
  wl0::Union{V,Nothing}=nothing,
  wd0::Union{V,Nothing}=nothing) where {T<:Real,V<:Vector{T},M<:Matrix{T}}

  WMM = WUM + WLM + WDM

  Eu, El, Ed = calc_evap(wu0, wl0, P, PET .* K; WLM, C)
  ET = Eu + El + Ed

  PE = max.(prcp - ET, 0)
  r, rim = calc_runoff(w0, PE; wm, B, IM)
  wu, wl, wd = Storage_Update(wu0, wl0, wd0, eu, el, ed, prcp .- e, r; WUM, WLM, WDM)

  return (r, rim, e, pe), (wu, wl, wd)
end


function sources(PE::T, R::T, S1::T, FR1::T, FR::T; params) where {T<:Real}
  (; SM, EX, KI, KG) = params
  # FR = PE / R
  if FR == 0
    RS, RI, RG, S = 0.0
    return RS, RI, RG, S, FR
  end

  SMM = SM * (1 + EX)
  S = S1 * FR1 / FR

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
  return RS, RI, RG, S, FR
end

# 主函数
function xaj(p_and_e::A, params::XAJ; return_state::Bool=false,
  warmup_length::Integer=365) where {T<:Real,A<:Array{T,3},M<:Matrix{T}}

  ## 需要修正ki, kg的参数
  (; KI, KG) = params
  if (KI + KG > 1)
    PRECISION = 1e-5
    scale = (1 - PRECISION) / (KI + KG)
    KI, KG = KI * scale, KG * scale
  end

  # 预热
  if warmup_length > 0
    _, _, wu0, wl0, wd0, s0, fr0, qi0, qg0 = xaj(p_and_e[1:warmup_length, :, :], params;
      return_state=true, warmup_length=0)
  else
    wu0, wl0, wd0 = T(0.5) .* um, T(0.5) .* lm, T(0.5) .* dm
    s0, fr0 = T(0.5) .* sm, fill(T(0.1), size(ex))
    qi0, qg0 = fill(T(0.1), size(ci)), fill(T(0.1), size(cg))
  end

  # 主循环
  inputs = p_and_e[warmup_length+1:end, :, :]
  n_steps, n_basins = size(inputs)[1:2]

  runoff_ims = zeros(T, n_steps, n_basins)
  rss = zeros(T, n_steps, n_basins)
  ris = zeros(T, n_steps, n_basins)
  rgs = zeros(T, n_steps, n_basins)
  es = zeros(T, n_steps, n_basins)

  wu, wl, wd = wu0, wl0, wd0
  s, fr = s0, fr0
  qi, qg = qi0, qg0

  for i in 1:n_steps
    (r, rim, e, pe), (wu, wl, wd) = generation(inputs[i:i, :, :];
      k, b, im, um, lm, dm, c,
      wu0=wu, wl0=wl, wd0=wd)
    (rs, ri, rg), (s, fr) = sources(pe, r; s0=s, fr0=fr, sm, ex, ki, kg)

    runoff_ims[i, :] = rim
    fact = 1 - im
    rss[i, :] = rs .* fact
    ris[i, :] = ri .* fact
    rgs[i, :] = rg .* fact
    es[i, :] = e
  end

  # 汇流计算
  qs = zeros(T, n_steps, n_basins)
  qt = zeros(T, n_steps, n_basins)

  qi, qg = qi0, qg0
  for i in 1:n_steps
    qi = linear_reservoir(ris[i, :], ci; last_y=qi)
    qg = linear_reservoir(rgs[i, :], cg; last_y=qg)
    qt[i, :] = rss[i, :] .+ runoff_ims[i, :] .+ qi .+ qg
  end

  for j in 1:n_basins
    lag = Int(round(l[j]))
    qs[1:lag, j] = qt[1:lag, j]

    for i in (lag+1):n_steps
      qs[i, j] = cs[j] * qs[i-1, j] + (1 - cs[j]) * qt[i-lag, j]
    end
  end

  q_sim = reshape(qs, n_steps, n_basins, 1)
  es_out = reshape(es, n_steps, n_basins, 1)

  return return_state ? (q_sim, es_out, wu, wl, wd, s, fr, qi, qg) : (q_sim, es_out)
end
