using SpecialFunctions: gamma

function _uh_gamma(τ::Real; α::T=2.5, θ::T=1.0) where {T<:Real}
  1 / (gamma(α) * θ^α) * τ^(α - 1) * exp(-τ / θ)
end

function uh_gamma(τ_max::Int; α::T=2.5, θ::T=1.0) where {T<:Real}
  map(τ -> _uh_gamma(τ; α, θ), 0:τ_max)
end

function uh_gamma(; α::T=2.5, θ::T=1.0) where {T<:Real}
  e = map(τ_max -> begin
      uh = map(τ -> _uh_gamma(τ; α, θ), 0:τ_max)
      1 - sum(uh)
    end, 6:24)
  DataFrame(τ_max=6:24, e)
end


"""
A -> B (I -> O)
B -> C 
"""
function uh_conv(I::AbstractVector{T}, uh::AbstractVector{T}) where {T<:Real}
  ntime = length(I)
  O = zeros(T, ntime)
  τ_max = length(uh) - 1 # index begin from 0

  for t in 1:ntime
    q = T(0)
    for τ = maximum(t + 1, 0):τ_max
      # τ >= maximum(t + 1, 0)
      q += I[t-τ] * uh[τ]
    end
    O[t] = q
  end
  return O
end
