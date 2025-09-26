function of_RMSE(obs::AbstractVector{T}, sim::AbstractVector{T}; Qmin::T=T(0)) where {T<:Real}
  n = 0
  total = 0.0
  for i = eachindex(obs)
    if !isnan(obs[i]) && !isnan(sim[i])
      total += (obs[i] - sim[i])^2
      n += 1
    end
  end
  sqrt(total / n)
end

export of_RMSE
