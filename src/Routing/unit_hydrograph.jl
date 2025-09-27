"""
## 第1个时刻
- `start = 1`: 无法传播, Qt = I[t-1] * uh[1] + ..., 
  Q[t] += I[t-τ] * uh[τ]
  t - τ >= 1       ==>   τ <= t - 1

- `start = 0`: 瞬时传播，Qt = I[t] * uh[1] + ..., 
  Q[t] += I[t-τ+1] * uh[τ]
  t - τ + 1 >= 1   ==>   τ <= t
"""
function conv_uh(I::AbstractVector{T}, uh::AbstractVector{T}; start::Int=1) where {T<:Real}
  ntime = length(I)
  Q = zeros(ntime)
  τ_max = length(uh)

  τ0 = 1 - start # (start = 1: τ0 = 0), (start = 1: τ0 = 1)
  for t = 1:ntime
    for τ = 1:τ_max # min(τ_max, t - start)
      k = t - τ + τ0
      k <= 0 && (k = 1) # 假设I是前期为I[1]
      Q[t] += I[k] * uh[τ]
    end
  end
  return Q
end


# 用于山坡汇流
function uh_gamma(τ::Real; α::T=2.5, θ::T=1.0)::T where {T<:Real}
  1 / (gamma(α) * θ^α) * τ^(α - 1) * exp(-τ / θ)
end

function normalize_uh(uh::Vector{T}) where {T<:Real}
  S::T = sum(uh)
  uh ./ S
end

# function uh_gamma(; α::T=2.5, θ::T=1.0) where {T<:Real}
#   e = map(τ_max -> begin
#       uh = uh_gamma.(1:τ_max; α, θ)
#       1 - sum(uh)
#     end, 6:24)
#   DataFrame(; τ_max=6:24, e)
# end


"""
线性水库单位线

# References
1. 沈冰，水文学原理，2008，P171, Eq. 9-105
> Eq. 9-109, 公式存在bug，`(n-1)`应是`(n-1)!`
"""
uh_linear(t::Real; K=1, dt=1, nseg::Int=1) = dt / (K * factorial(nseg - 1)) * pow(t / K, nseg - 1) * exp(-t / K)


function linear_reservoir_uh(I::AbstractVector; K=1, dt=1, nseg::Int=3, τ_max=30)
  uh = uh_linear.(1:τ_max; K, dt, nseg)
  conv_uh(I, uh)
end


function guess_uh(τ_max=10; atol=0.01, maxn=100, fun=uh_linear, param::NamedTuple)
  count = 0
  while true
    uh = fun.(1:τ_max; param...) # K, dt, n=nseg
    if 1 - sum(uh) < atol
      return τ_max, sum(uh), uh
    end

    if count >= maxn
      @warn("Not converge! maxn = $maxn reached!")
      return τ_max, sum(uh), uh
    end
    count += 1
    τ_max += 1
  end
end
