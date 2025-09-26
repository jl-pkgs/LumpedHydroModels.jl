# 线性水库汇流：应用于地下水汇流
# linear_reservoir
function linear_reservoir_low(I::AbstractVector; CS=NaN, K=1.0, dt=1, nseg::Int=3)
  isnan(CS) && (CS = (2K - dt) / (2K + dt))

  n = length(I)
  Q = zeros(n, nseg + 1)
  Q[:, 1] .= I
  Q[1, :] .= I[1]

  # 设计一种更加节省内存的结构
  @inbounds for i = 2:n
    for j = 2:nseg+1
      I1 = Q[i-1, j-1]
      I2 = Q[i, j-1]
      Q1 = Q[i-1, j]

      Q[i, j] = CS * Q1 + (1 - CS) * (I1 + I2) / 2
    end
  end
  return Q
end


function linear_reservoir(I::AbstractVector; CS=NaN, K=1.0, dt=1, nseg::Int=3)
  isnan(CS) && (CS = (2K - dt) / (2K + dt))

  n = length(I)
  Qtemp = zeros(2, nseg + 1)  # 只保存两行：一行表示t-1时刻，一行表示t时刻
  Qtemp[1, :] .= I[1]         # 初始化：将第一个时刻的所有河段流量设为I[1]

  Qout = zeros(n)             # 储存每个时间最后一段的出流结果
  Qout[1] = I[1]

  @inbounds for i = 2:n
    Qtemp[2, 1] = I[i]        # update Q_{t}
    for j = 2:nseg+1
      I1 = Qtemp[1, j-1]
      I2 = Qtemp[2, j-1]
      Q1 = Qtemp[1, j]
      Qtemp[2, j] = CS * Q1 + (1 - CS) * (I1 + I2) / 2
    end

    Qout[i] = Qtemp[2, end]     # 只储存最后一段的Q
    for j in 1:nseg+1
      Qtemp[1, j] = Qtemp[2, j]
    end
    # Qtemp[1, :] .= Qtemp[2, :]  # 更新t-1的数据
  end
  return Qout
end



function linear_reservoir_uh(I::AbstractVector; K=1, dt=1, nseg::Int=3, n_uh=30)
  uh = coef_Pm_linear.(0:n_uh; K, dt, n=nseg)
  conv_uh(I, uh)
end


function guess_uh(n=10; atol=0.01, maxn=100, fun=coef_Pm_linear, param::NamedTuple)
  count = 0
  while true
    uh = fun.(0:n; param...) # K, dt, n=nseg
    if 1 - sum(uh) < atol
      return n, sum(uh), uh
    end
    if count >= maxn
      @warn("Not converge! maxn = $maxn reached!")
      return n, sum(uh), uh
    end
    count += 1
    n += 1
  end
end


"""
线性水库单位线

# References
1. 沈冰，水文学原理，2008，P171, Eq. 9-105
> Eq. 9-109, 公式存在bug，`(n-1)`应是`(n-1)!`
"""
coef_Pm_linear(t; dt, K, n) = dt / (K * factorial(n - 1)) * pow(t / K, n - 1) * exp(-t / K)

export linear_reservoir, linear_reservoir_uh,
  linear_reservoir_low,
  coef_Pm_linear,
  guess_uh
