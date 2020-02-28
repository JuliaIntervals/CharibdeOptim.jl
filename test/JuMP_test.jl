using Test, JuMP
using CharibdeOptim
using IntervalArithmetic

#= @testset "Using JuMP syntax for Constrained Optimisation" begin
      model = Model(with_optimizer(CharibdeOptim.Optimizer))
      @variable(model, -4 <= x <= 4)
     @NLconstraint(model, x+y<=4)
      @NLobjective(model, Max, -(x-4)^2-(y-4)^2)
      optimize!(model)

      @test JuMP.termination_status(model) == MOI.OPTIMAL

      @test JuMP.objective_value(model) ⊆ -8.01 .. -7.99

      @test JuMP.value(y) ⊆ (1.99 .. 2.01)
end
=#

@testset "Using JuMP syntax by using only one worker " begin
      model = Model(with_optimizer(CharibdeOptim.Optimizer, workers = 1))
      @variable(model, 1<=x<=2)
      @variable(model, 1<=y<=2)
      @NLobjective(model, Min, x^2+y^2)
      optimize!(model)

      @test JuMP.termination_status(model) == MOI.OPTIMAL
      @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
      @test JuMP.objective_value(model) ≈ 2.0
      @test JuMP.value(x) ≈ 1.0
      @test JuMP.value(y) ≈ 1.0
end

@testset "Optimising difficult problem using JuMP" begin
      model = Model(with_optimizer(CharibdeOptim.Optimizer))

      @variable(model, 2 <= x1 <= 3)
      @variable(model, 3 <= x2 <= 4)
      @variable(model, 9 <= x3 <= 14)
      @variable(model, 2 <= x4 <= 3)
      @variable(model, 3 <= x5 <= 4)
      @variable(model, 9 <= x6 <= 14)
      @variable(model, 2 <= x7 <= 3)
      @variable(model, 3 <= x8 <= 4)
      @variable(model, 9 <= x9 <= 14)

      @NLobjective(model, Min, x1^2 + x2^2 + x3^4 - x4^7 - 200x5 - x6^5 - x7^9 + x8^5 - 8x9^3)

      optimize!(model)

      @test JuMP.termination_status(model) == MOI.OPTIMAL
      @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
      @test JuMP.objective_value(model) ≈ -575629.0
      @test JuMP.value(x1) ⊆ 1.99 .. 2.01
      @test JuMP.value(x2) ⊆ 2.99 .. 3.01
      @test JuMP.value(x3) ⊆ 8.99 .. 9.01
      @test JuMP.value(x4) ⊆ 2.99 .. 3.01
      @test JuMP.value(x5) ⊆ 3.99 .. 4.01
      @test JuMP.value(x6) ⊆ 13.99 .. 14.01
      @test JuMP.value(x7) ⊆ 2.99 .. 3.01
      @test JuMP.value(x8) ⊆ 2.99 .. 3.01
      @test JuMP.value(x9) ⊆ 13.99 .. 14.01
end
