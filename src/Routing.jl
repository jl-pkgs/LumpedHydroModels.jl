# 线性水库
function linear_reservoir(x::V, weight::V; last_y::Union{V,Nothing}=nothing) where {T<:Real,V<:Vector{T}}
  isnothing(last_y) && (last_y = fill(T(0.001), size(weight)))
  return weight .* last_y .+ (1 .- weight) .* x
end


# 单位线卷积
function uh_conv(x::A, uh::A) where {T<:Real,A<:Array{T,3}}
  time_len, batch_size, _ = size(x)
  outputs = zeros(T, size(x))

  for i in 1:batch_size
    uh_vec = uh[:, i, 1]
    input_vec = x[:, i, 1]
    conv_result = conv(input_vec, uh_vec)
    outputs[:, i, 1] = conv_result[1:time_len]
  end
  return outputs
end


# Gamma单位线
function uh_gamma(a::A, theta::A, len_uh::Integer=15) where {T<:Real,A<:Array{T,3}}
  m = size(a)
  len_uh > m[1] && throw(ArgumentError("单位线长度过长"))

  aa = max.(0, a[1:len_uh, :, :]) .+ T(0.1)
  theta_adj = max.(0, theta[1:len_uh, :, :]) .+ T(0.5)

  t = reshape(repeat(collect(T(0.5):1:(len_uh-T(0.5))), outer=(1, m[2])), (len_uh, m[2], 1))

  denominator = gamma.(aa) .* (theta_adj .^ aa)
  w = (1 ./ denominator) .* (t .^ (aa .- 1)) .* exp.(-t ./ theta_adj)
  w = w ./ sum(w, dims=1)
  return w
end
