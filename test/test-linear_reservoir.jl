using Test

@testset "linear_reservoir" begin
  nseg = 3
  dt = 1
  K = 3.0
  param = (; K, dt, n=nseg)

  I = [2300, 2340, 2400, 2480, 2520, 2600, 2700, 2810, 2900, 3010, 3190, 3350, 3600, 4500, 6000, 7000,
    7520, 8100, 8800, 9300, 9500, 9700, 9700, 9650, 9550, 9430, 9250, 9100, 9070, 9000]

  uh = guess_uh(10; param, atol=0.01)[end]
  # Q_test = conv_uh(I, uh)[:, end]
  Q_test = linear_reservoir_uh(I; K, dt, nseg, n_uh=30)
  Q_true = linear_reservoir_low(I; K, dt)[:, end]
  @test of_RMSE(Q_true, Q_test) < 20

  n, accu, uh = guess_uh(10; atol=0.01, param)
  @test abs(accu - 0.9906) < 0.01


  ## test performance
  nseg = 10
  I = rand(1000)
  @time Q_kong = linear_reservoir_low(I; K, dt)[:, end]
  @time Q_shi = linear_reservoir(I; K, dt)
  @test all(abs.(Q_kong - Q_shi) .< 1e-8)
end

# using Plots
# plot(Q_test, label="Q_test")
# plot!(Q_true, label="Q_true")
# uh = guess_uh(10; param, atol=0.01)
