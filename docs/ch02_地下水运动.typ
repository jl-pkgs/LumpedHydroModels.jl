#import "@preview/mannot:0.3.0": *
#import "@local/modern-cug-report:0.1.2": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "")

#set par(leading: 1.24em, spacing: 1.14em)

#let delta(x) = $Delta #x$
#let boxed = markrect.with(outset: (x: 0.3em, y: 0.3em))

= 1 *基本概念*

在承压含水层中，$S = S_s times M$，$S_s$为比储水系数，M为含水层厚度，$S$比储水系数 。

在非承压含水层中，$S = S_y$。

== 1.1 承压含水层


$ delta(W) = S_s b thin delta(h) $

$ arrow(q) = -k b thin nabla h $

$
  pdv(W, t) = - nabla dot arrow(q) = k b thin nabla^2 h, ==>
$

$ boxed(pdv(h, t) = k/S_s pdv(h, z, 2)) $

其中$h$为水位，向下为负。

== 1.2 非承压含水层

#let Seff = $S_"eff"$

$ delta(W) = S_y thin delta(h) $

$ arrow(q) = -k h thin nabla h $

$ pdv(W, t) = - nabla dot arrow(q) $

$ S_y pdv(h, t) = nabla (k h nabla h) $

若假设$h approx D$（平均厚度），则$pdv(h, t) = (k D) / S_y pdv(h, z, 2)$。定义“等效比储水系数”$Seff = S_y / D, ==>$

$ boxed(pdv(h, t) = k / Seff pdv(h, z, 2)) $ <eq_unconfined>

=== 1.2.1 公式求解

由式#[@eq_unconfined]作为出发点，求解偏微分方程。

令：

$ h(x, t) = X(x) T(t) $

带入方程，

$ X T' = k D / S_y X'' T ==> T'/T = (k D) / S_y X'' / X = - lambda $

$lambda$的由来：上式想成立，必须满足左右都等于常数。

$ T' + lambda T = 0 ==> T(t) = e^{- lambda t} $

$ X'' + (S_y lambda )/ (k D) X = 0 $

令$mu^2 = (S_y lambda) / (k D)$，则：

$
  lambda = (mu^2 k D) / S_y
$

$
  X'' + mu^2 X = 0 ==> X(x) = A cos(mu x) + B sin (mu x)
$

#box-blue[
  *边界条件：*
  1. $h(0, t) = 0$ 以河道为水平面进行建模，$==> X(0) = 0$

    $X(0) = A cos(0) + B sin(0) = A ==> boxed(X(x) = B sin (mu x))$

  2. $pdv(h, x)(B, t) =0 ==> X'(B) = 0$

    $X'(x) = B mu cos(mu x) ==> X'(B) = B mu cos(mu x) = 0$，也即是$cos(mu B) = 0$，$==> mu B = pi/2, (3 pi)/2, (5 pi)/2$。

    $
      boxed(mu_n = ((2n - 1) pi) / (2B))
    $
]

$ lambda = (mu^2 k D) / S_y = (k D) / S_y [((2n - 1) pi) / (2B)]^2 $

$
  h_n(x, t) &= X(x) T(t) = B sin(mu x) e^{- lambda t} \
  &= B sin (((2n - 1) pi) / (2B) x ) e^{- (k D) / S_y [((2n - 1) pi) / (2B)]^2 t}
$


因为上式是线性齐次偏微分方程，满足叠加原理，所有特解的线性组合仍是解。

$
  h(x, t) &= sum_(n = 1)^(oo) C_n thin h_n (x, t) \
  &= sum_(n = 1)^(oo) C_n thin B sin (((2n - 1) pi) / (2B) x ) e^{- (k D) / S_y [((2n - 1) pi) / (2B)]^2 t}
$

=== 1.2.2 解的简化

另一个初始条件：$h(x, 0) = D$，将t =0，带入上式：

$ D = sum_(n = 1)^(oo) C_n thin B sin (((2n - 1) pi) / (2B) x ) $

$
  phi.alt_n(x) = sin( ((2n - 1) pi x) / (2B) ), quad n = 1, 2, 3, dots
$

两边同乘以 $phi.alt_m(x)$ ，在 [0, B] 上积分：

$
  integral_0^B D phi.alt_m(x) d x = sum_(n = 1)^(oo) C_n integral_0^B phi.alt_n(x) phi.alt_m(x) d x
$

$
  integral_0^B D phi.alt_m(x) d x = C_m integral_0^B [phi.alt_m(x)]^2 d x \
  => C_m = frac( integral_0^B D phi.alt_m(x) d x, integral_0^B [phi.alt_m(x)]^2 d x )
$

$
  C_n = frac(2, B) integral_0^B D sin( ((2n - 1) pi x)/(2B) ) d x = frac(4D, (2n - 1) pi )
$

#let ecomp = $ -(k_0 D) / n_e ( ((2n-1) pi) / (2B) )^2 t $

$
  H(x, t) = sum_(n=1)^(oo) frac(4D, (2n-1)pi) sin((2n-1)pi x / (2B))
  e^{-(k_0 D / n_e) ((2n-1)pi / (2B))^2 t}
$

$
  (∂H)/(∂x) &= sum_(n=1)^(oo) frac(4D, (2n-1)pi) frac((2n-1)pi, 2B)
  cos((2n-1)pi x / (2B)) e^{ecomp} \
  &= sum_(n=1)^(oo) frac(2D, B) cos((2n-1)pi x / (2B)) e^{ecomp}
$

$x = 0$处：

$
  (∂H)/(∂x) |_(x=0) &= sum_(n=1)^(oo) frac(2D, B) e^{ecomp} \
  &= frac(2D, B) sum_(n=1)^(oo) e^{-(k_0 D pi^2 (2n-1)^2 t)/(4 n_e B^2)}
$

$
  q(t) = k_0 D frac(2D, B) sum_(n=1)^(oo) e^{-(k_0 D pi^2 (2n-1)^2 t)/(4 n_e B^2)}
$

在实际应用中，指数衰减项中$n=1$衰减最慢，主导长期行为：

$
  q(t) approx k_0 D frac(2D, B) e^{-(k_0 D pi^2 t)/(4 n_e B^2)} = frac(2 k_0 D^2, B) e^{-(pi^2 k_0 D t)/(4 n_e B^2)}
$

引入经验权重因子，最终解：

$ boxed(q(t) = 2 k_0 p frac(D^2, B) e^{- (pi^2 k_0 p D t) / (4 n_e B^2) }) $

#pagebreak()

= 2 *注水试验*

