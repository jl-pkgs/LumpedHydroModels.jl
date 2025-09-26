# 土壤蓄水量更新, TODO 存在另一种写法
function UpdateW_Zhong!(P::T, state::StateXAJ; param::XAJ) where {T<:Real}
  (; R, WU, WL, WD, EU, EL, ED) = state
  (; WUM, WLM, WDM) = param

  # 蒸发过程已经考虑了W的限制，因此这里主要处理的是水多的情景
  if WU + P - EU - R <= WUM
    WU += P - EU - R
    WL -= EL
    WD -= ED
    exceed_UL = 0.0
  else
    WU = WUM
    exceed_UL = WU + P - EU - R - WUM # 多余出来的水

    if WL - EL + exceed_UL <= WLM
      WL = WL - EL + exceed_UL
      WD -= ED
      exceed_LD = 0.0
    else
      WL = WLM
      exceed_LD = WL + exceed_UL - EL - WLM

      if WD - ED + exceed_LD <= WDM
        WD = WD - ED + exceed_LD
      else
        WD = WDM # TODO: 多出来的水去哪了？放到地下？
      end
    end
  end
  W = WU + WL + WD
  @pack! state = WU, WL, WD, W
  # return WU, WL, WD, W
end


# 书上的做法
function UpdateW(P::T, state::StateXAJ; param::XAJ) where {T<:Real}
  (; R, WU, WL, WD, ET) = state
  (; WUM, WLM, WDM) = param

  Δ = P - ET - R
  WU += Δ

  if WU > WUM # 供水充足
    WL += WU - WUM
    WU = WUM

    if WL > WLM
      WD += WL - WLM
      WL = WLM
    end
  end

  if WU < 0 # 供水不足
    WL += WL
    WU = 0.0

    if WL < 0
      WD += WL # TODO: 若WD < 0，应该如何处理？
      WL = 0.0
    end
  end

  W = WU + WD + WL
  @pack! state = WU, WL, WD, W
end
