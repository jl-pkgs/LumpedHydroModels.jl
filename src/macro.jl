macro DefOutputs(StateStruct, extra_fields...)
  state_name = StateStruct  # 获取名称（Symbol）
  output_name = Symbol(replace(string(state_name), r"^State" => "Outputs"))

  state_type = getfield(__module__, state_name)
  fieldnames_list = collect(fieldnames(state_type))

  field_expressions = []

  # 1. ntime 字段（必须第一个，因为其他字段依赖它）
  push!(field_expressions, :(ntime::Int = 100))

  # 2. 从 StateStruct 复制所有字段 → 转为 Vector{FT}
  for fname in fieldnames_list
    push!(field_expressions, :($fname::Vector{FT} = zeros(FT, ntime)))
  end

  # 3. 添加额外字段（如 Rsim）
  for extra in extra_fields
    if isa(extra, Symbol)
      # 默认初始化为 zeros(FT, ntime)
      push!(field_expressions, :($extra::Vector{FT} = zeros(FT, ntime)))
    elseif isa(extra, Expr) && extra.head === :(=)
      # 支持自定义默认值，例如 Rsim = something
      push!(field_expressions, extra)
    else
      error("额外字段格式错误：$(extra)")
    end
  end

  # 构建完整的 struct 定义，并显式继承 AbstractHydroOutputs{FT}
  return esc(quote
    @with_kw mutable struct $output_name{FT} <: AbstractHydroOutputs{FT}
      $(field_expressions...)
    end
  end)
end
