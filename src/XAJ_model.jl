@bounds @units @with_kw mutable struct XAJ{FT} <: AbstractHydroModel{FT}
  K::FT = 0.8 | (0.2, 1.5) | "-"        # 蒸发折算系数
  C::FT = 0.1 | (0.05, 0.20) | "-"      # 深层蒸发系数

  IM::FT = 0.1 | (0.01, 0.3) | "-"      # 不透水面积比例

  WUM::FT = 10.0 | (5.0, 20.0) | "mm"   # 上层土壤平均蓄水容量
  WLM::FT = 30.0 | (10.0, 90.0) | "mm"  # 下层
  WDM::FT = 6.0 | (10.0, 120.0) | "mm"  # 深层

  SM::FT = 5.0 | (10.0, 60.0) | "mm"    # 自由水蓄水容量

  B::FT = 0.4 | (0.1, 0.6) | "-"        # 蓄水容量曲线的指数参数
  EX::FT = 0.7 | (0.5, 2.0) | "-"       # 自由水蓄水容量曲线指数

  KI::FT = 0.3 | (0.01, 0.7) | "-"      # 壤中流出流系数
  KG::FT = 0.3 | (0.01, 0.7) | "-"      # 地下水出流系数

  CS::FT = 0.7 | (0.01, 0.9) | "-"      # 汇流参数
  CI::FT = 0.7 | (0.5, 0.99) | "-"      # 壤中流，线性水库汇流参数
  CG::FT = 0.98 | (0.95, 0.998) | "-"   # 地下径流，线性水库汇流参数
end


@with_kw mutable struct StateXAJ{FT} <: AbstractHydroState{FT}
  ET::FT = 0.0
  EU::FT = 0.0
  EL::FT = 0.0
  ED::FT = 0.0

  PE::FT = 0.0    # 净雨, P - ET

  W::FT = 0.0
  WU::FT = 0.0
  WL::FT = 0.0
  WD::FT = 0.0

  R::FT = 0.0     # 透水界面径流
  R_IM::FT = 0.0  # 不透水界面径流
  RS::FT = 0.0
  RI::FT = 0.0
  RG::FT = 0.0
  S::FT = 0.0     # 自由水蓄量
end

@with_kw mutable struct OutputsXAJ{FT} <: AbstractHydroOutputs{FT}
  ntime::Int = 100
  ET::FT = zeros(ntime)
  EU::FT = zeros(ntime)
  EL::FT = zeros(ntime)
  ED::FT = zeros(ntime)

  PE::FT = zeros(ntime)    # 净雨, P - ET

  W::FT = zeros(ntime)
  WU::FT = zeros(ntime)
  WL::FT = zeros(ntime)
  WD::FT = zeros(ntime)

  R::FT = zeros(ntime)     # 透水界面径流
  R_IM::FT = zeros(ntime)  # 不透水界面径流
  RS::FT = zeros(ntime)
  RI::FT = zeros(ntime)
  RG::FT = zeros(ntime)
  Rsim::FT = zeros(ntime)
  S::FT = zeros(ntime)     # 自由水蓄量
end



function run_XAJ(P::Vector{T}, PET::Vector{T}; param::XAJ{T}) where {T<:Real}
  ntime = length(P)
  output = OutputsXAJ{T}(; ntime)
  (; RS, RI, RG) = output

  state = StateXAJ{FT}()

  for t = 2:ntime
    _P = P[t]
    _PET = PET[t]

    Evapotranspiration!(state, _P, _PET; param)
    Runoff(state; param)
    UpdateW(state, _P; param)

    Runoff_divide3S(state, FR[t-1]; param)
    output[t] = state
  end
  output.Rsim = RS + RI + RG
  output
end
