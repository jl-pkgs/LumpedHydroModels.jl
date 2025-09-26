@testset "ModelParams update!" begin
  FT = Float64
  model = ModelNetwork{FT}(; models=[XAJ(), XAJ()])

  params = Params(model)
  paths = [
    [:models, 1, :K],
    [:models, 2, :WDM],
  ]
  parvalues = [0.5, 90.0]

  @time update!(model, paths, parvalues; params)  
  @test model.models[1].K == 0.5
  @test model.models[2].WDM == 90.0
end

## 同时增加功能，哪些可以设置为共同参数，哪些设置为变异参数
# vars_common = 
# vars_vary = [:KI, :KG]
# par_common = DataFrame(; name=[:K, :C, :IM], value=[0.9, 0.2, 0.3])
# par_vary = DataFrame(; path=paths, value=parvalues)
