# 测试：蒸发、updateW，产流
using RTableTools, LumpedHydroModels, DataFrames

d = fread("$indir/data/dat_Table2-11.csv")
replace_missing!(d, 0.0)

@testset "run_XAJ" begin
  P = d.P
  PET = d.Eo

  param = XAJ(; K=0.95, C=0.14, B=0.3,
    WUM=15.0, WLM=85.0, WDM=20.0,
    IM=0.0)
  state = StateXAJ(; WU=0.0, WL=2.2, WD=20.0)

  r = run_XAJ(P, PET; param, state)
  r = round.(r, digits = 1)
  @test sum(r.R) ≈ 39.3
end

# 汇总：各输出列的总量（单行宽格式）
# values = map(sum, eachcol(r))
# values = round.(values, digits=1)
# s = NamedTuple{Tuple(Symbol.(names(r)))}(values)
# DataFrame([s])

# r[:, [:EU, :EL, :ED]]
# r[:, [:WU, :WL, :WD]]
# r
