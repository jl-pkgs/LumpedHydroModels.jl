#import "@local/modern-cug-report:0.1.2": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#set par(leading: 1.24em, spacing: 1.14em)

= 1 *新安江模型主要原理*

#include("Figure1_蓄水曲线.typ")

#box-blue[
  #set par(leading: 0.6em, spacing: 0.5em)
  涨水过程中，所有格点的土壤含水层、蓄水量均匀上涨。对于某一点$P(alpha, W)$，饱和部分的比例$alpha$，
  未饱和部分比例为$1 - alpha$。
]

#let wm = $"WM'"$
#let WM = $"WM"$
#let WMM = $"WMM"$

#let Alpha = $bold(alpha)$
#let delta(x) = $Delta #x$

$ Alpha = 1 - (1 - wm / WMM)^b $

$
  WM = integral_0^WMM (1 - Alpha) thin d"WM'"
  = -WMM/(b+1) (1 - wm/WMM)^(b+1) |_0^WMM
  = WMM / (1 + b)
$ <eq_WM>

其中，$wm$为流域所有网格中最大土壤含水量；$Alpha$为包气带蓄水容量$≤wm$的比例。
$Alpha_0$为蓄满的比例，$a$为最大土壤含水量。
// #box-blue[
//   右侧为未饱和的部分。
//   左侧是已经饱和的部分。
// ]

$ Alpha_0 = 1 - (1 - a / WMM)^b $

$
  W = integral_0^a (1 - Alpha ) thin d"WM'" = WM [1 - (1 - a / WMM)^(b + 1)]
$ <eq_w>

$ WM (1 - a / WMM)^(b + 1) = WM - W $

写成$a$的形式：
$
  (1 - a / WMM)^(b + 1) = 1 - W / WM ==> 1 - a/WMM = (1 - W / WM)^(1/(b+1)) ==> \
  a = WMM [1 - (1 - W / WM)^(1/(b+1))]
$

== 1.1 产流过程

$
  R = integral_a^(a + "PE") Alpha thin d"WM'"
  = integral_a^(a + "PE") [1 - (1 - wm / WMM)^b] thin d"WM'" \
  = "PE" - integral_a^(a + "PE") (1 - wm / WMM)^b thin d"WM'" \
  = "PE" + WMM/(b+1) (1 - wm/WMM)^(b+1) |_a^(a+"PE") \
  = "PE" + WMM/(b + 1) [(1 - (a + "PE")/WMM)^(b+1) - (1 - a/WMM)^(b+1)] \
  = "PE" + WM [(1 - (a + "PE")/WMM)^(b+1) - (1 - a/WMM)^(b+1)] \
  = "PE" + WM (1 - (a + "PE")/WMM)^(b+1) + W - WM
$

// #pagebreak()
== 1.2 径流分割

#let mm = $"mm"$
#let EX = $"EX"$
#let AU = $"AU"$
#let FR = $"FR"$
#let PE = $"PE"$
#let RS = $"RS"$
#let RI = $"RI"$
#let RG = $"RG"$

对自由自由水蓄量进行径流分割。#Blue[注意这里建模时$Alpha$指的是产流面积上的比例，后续在计算$RS$时要乘以产流面积$FR$]。因为只有产流产流区域才有自由水蓄量。

$ Alpha = 1 - (1 - S / S_mm)^EX $

$ S_m = S_mm / (1 + EX) $

$ S = S_m [1 - (1 - AU / S_mm)^(EX + 1)] $

$ AU = S_mm [1 - (1 - S / S_m)^(1/(EX + 1))] $ <eq_au_raw>

#v(-0.8em)

#box-red[
  新安江模型将土壤层的蓄水分为两种：
  - *张力水*（或称“土壤含水量”），土壤中被毛管力吸附的水分，低于田间持水量，不直接产流
  - *自由水*，超过田间持水量的重力水，可自由下渗或形成径流，可产流（地表/地下），#Blue[仅存在于已蓄满（产流）区域]

  *未产流区域有没有“自由水蓄量”。*
]

因此，式#[@eq_au_raw]需要考虑产流面积的不同，进行修正。

$ AU = S_mm [1 - (1 - (S times FR_1 \/ FR) / S_m)^(1/(EX + 1))] $

$
  (1 - AU / S_mm )^(EX + 1) = 1 - (S times FR_1 \/ FR) / S_m
$

$FR_1, FR$分别是上一时段和本时段的产流面积。

$
  "RS"/FR = integral_AU^{AU+PE} Alpha thin d S
  = integral_AU^{AU+PE} [1 - (1 - S / S_mm)^EX] d S \
  = PE + S_mm / (EX + 1) (1 - S / S_mm)^(EX + 1) |_AU^{AU+PE} \
  = PE + S_m (1 - S / S_mm)^(EX + 1) |_AU^{AU+PE} \
  = PE + S_m (1 - (PE + AU) / S_mm)^(EX + 1) - S_m (1 - AU / S_mm)^(EX + 1) \
  = PE + S_m (1 - (PE + AU) / S_mm)^(EX + 1) - S_m + S times FR_1 / FR
$

径流，地表没有流走的部分，加到自由水蓄量$S$中。

$ S = S_1 FR_1 / FR + (R - RS) / FR $

$
  RI = K_I times S times FR \
  RG = K_G times S times FR \
$

$ S = S(1 - K_I - K_G) $

= 2 *汇流方法*

== 2.1 线性水库

$ dv(S, t) = I - Q = (I_1 + I_2) / 2 - (Q_1 + Q_2) / 2 $

$
  Q = S / K, Delta S = S_2 - S_1 = K(Q_2 - Q_1)
$

因此，

$
  K(Q_2 - Q_1) / (Delta t) = (I_1 + I_2) / 2 - (Q_1 + Q_2) / 2 \
  2K(Q_2 - Q_1) = Delta t [(I_1 + I_2) - (Q_1 + Q_2) ] \
  (2K + delta(t)) Q_2 = (2K - delta(t)) Q_1 + delta(t) (I_1 + I_2) \ 
  Q_2 = (2K - delta(t))/ (2K + delta(t)) Q_1 + delta(t) / (2K + delta(t)) (I_1 + I_2) \
$

$
  Q_2 = C_s Q_1 + (1 - C_S) (I_1 + I_2) / 2, C_s = (2K - delta(t))/ (2K + delta(t))
$

== 2.2 滞后线性水库

1. 水流像“快递”一样，从产流点到出口需要固定时间 τ ，途中不堆积、不变形。

2. 出口处有一个“理想水库”，其出流与蓄水量成正比：$Q = S / K$

$ dv(S, t) = I(t) - Q(t) $

$ K dv(Q(t), t) + Q(t) = I(t) $

产流 R(t) 先经过滞时 τ → 变成 R(t−τ)，R(t−τ) 再进入线性水库 → 输出 Q(t)

$ K dv(Q(t), t) + Q(t) = R(t - tau) $

$ K (Q_i - Q_{i-1}) / delta(t) + Q_i = R_{i-L} $

$
  Q_i = underbrace(frac(K, K + Delta t), c_s) Q_{i-1} +
  underbrace(frac(Delta t, K + Delta t), 1 - c_s) R_{i - L}
$

$ Q_i = c_s Q_{i-1} + (1 - c_s) R_{i - L} $

== 2.3 单位线

2参数gamma分布，作为单位线：

#let uh = $bold("uh")$
$
  bold(uh)(tau) = 1 / (Gamma(a) theta^a) tau^(a - 1) e^(- tau/theta)
$

$a$是形状参数，默认为2.5；$theta$是时间参数，1h。

众数：$(a - 1) theta$，方差：$a theta^2$

$
  q(t) = integral_0^{tau_max} uh(tau) thin R(t - tau) thin d tau
$

#let q = $gamma$
#let tmax = "tmax"

*uh从1开始*
$
  Q_1 & = r_1 uh_1 \
  Q_2 & = r_2 uh_1 + r_1 uh_2 \
  Q_3 & = r_3 uh_1 + r_2 uh_2 + r_1 uh_3 \
  Q_4 & = r_4 uh_1 + r_3 uh_2 +r_2 uh_3 + r_1 uh_4 \
  Q_t & = r_t uh_1 + r_{t-1} uh_2 + ... + r_{t - tmax + 1} uh_{tmax}
$

*uh从0开始*
$
  Q_1 & = r_1 uh_0 \
  Q_2 & = r_2 uh_0 + r_1 uh_1 \
  Q_3 & = r_3 uh_0 + r_2 uh_1 + r_1 uh_2 \
  Q_4 & = r_4 uh_0 + r_3 uh_2 +r_2 uh_2 + r_1 uh_3 \
  Q_t & = r_t uh_0 + r_{t-1} uh_1 + ... + r_{t - tmax + 1} uh_{tmax-1} + r_{t - tmax} uh_{tmax}
$
