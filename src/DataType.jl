import DataFrames: DataFrame
export ModelNetwork
export to_mat, DataFrame;

abstract type AbstractHydroModel{FT} <: AbstractModel{FT} end
abstract type AbstractHydroState{FT} <: AbstractModel{FT} end
abstract type AbstractHydroOutputs{FT} <: AbstractModel{FT} end


@with_kw mutable struct ModelNetwork{FT} <: AbstractHydroModel{FT}
  models::Vector{AbstractHydroModel{FT}}
  # routing::Vector      # 汇流参数
  # info_node::DataFrame # 河道节点
end


function Base.getindex(x::AbstractHydroState, names::Vector{Symbol})
  [getfield(x, name) for name in names]
end

function Base.getindex(x::AbstractHydroState, index::Vector{Int})
  names = fieldnames(typeof(x))[index]
  [getfield(x, name) for name in names]
end

function Base.setindex!(res::AbstractHydroOutputs, r::AbstractHydroState, i::Int64)
  fields = fieldnames(typeof(r))
  @inbounds for f in fields
    getfield(res, f)[i] = getfield(r, f)
  end
  return res
end


function to_mat(res::AbstractHydroOutputs{T}) where {T<:Real}
  TYPE = typeof(res)
  names = fieldnames(TYPE)[2:end] |> collect
  data = map(i -> getfield(res, i), names)
  data = cat(data..., dims=2)
  data, names
end

function DataFrame(res::AbstractHydroOutputs{T}) where {T<:Real}
  data, names = to_mat(res)
  DataFrame(data, names)
end
