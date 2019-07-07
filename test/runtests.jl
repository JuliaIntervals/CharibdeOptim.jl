using Distributed
addprocs(2)
@everywhere using CharibdeOptim
@everywhere using IntervalArithmetic
using Test, JuMP

@testset "CharibdeOptim tests " begin

      @testset "Optimising by Interval Bound and Contract Algorithm" begin
            (global_min, minimisers)= ibc_minimise((x,y)->x^2 + y^2, IntervalBox(2..3, 3..4))
            @test global_min ⊆ 13 .. 13.01
            @test minimisers[1] ⊆ (2.0 .. 2.001) × (3.0 .. 3.001)
      end

      @testset "Optimising by Charibde: A hybrid approach" begin
            (global_min, minimisers) = charibde_min((x,y)->x^2+y+1, IntervalBox(1..2, 2..3))
            @test global_min ⊆ 4.0 .. 4.01
            @test minimisers[1] ⊆ (1..1.001) × (2..2.001)

            (global_min, minimisers, info)= charibde_min((x,y)->x^2 + y^2, IntervalBox(2..3, 3..4))
            @test global_min ⊆ 13 .. 13.01
            @test minimisers[1] ⊆ (2.0 .. 2.001) × (3..3.001)
      end

      @testset "JuMP syntax" begin
            model = Model(with_optimizer(CharibdeOptim.Optimizer))
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

end
