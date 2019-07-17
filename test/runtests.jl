using CharibdeOptim
using IntervalArithmetic
using Test, JuMP


@testset "Optimising by Interval Branch and Contract Algorithm" begin
      (global_min, minimisers)= ibc_minimise(X->((x,y)= X;x^2 + y^2), IntervalBox(2..3, 3..4))
      @test global_min ⊆ 13 .. 13.01
      @test minimisers[1] ⊆ (2.0 .. 2.001) × (3.0 .. 3.001)
end


@testset "Optimising by Charibde (A hybrid approach) using only one worker" begin

      (global_min, minimisers, info) = charibde_min(X->((x,y)=X;x^2+y+1), IntervalBox(1..2, 2..3), workers = 1)
      @test global_min ⊆ 4.0 .. 4.01
      @test minimisers[1] ⊆ (1..1.001) × (2..2.001)

      (global_min, minimisers, info)= charibde_min(X->((x,y)=X;x^2 + y^2), IntervalBox(2..3, 3..4), workers = 1)
      @test global_min ⊆ 13 .. 13.01
      @test minimisers[1] ⊆ (2.0 .. 2.001) × (3..3.001)
end


@testset "Using JuMP syntax by having only one worker " begin
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


@testset "Using JuMP syntax by having two workers" begin
      model = Model(with_optimizer(CharibdeOptim.Optimizer))  # By default workers is set to 2
      @variable(model, 1<=x<=2)
      @variable(model, 2<=y<=3)
      @NLobjective(model, Min, x^3+2y+5)
      optimize!(model)

      @test JuMP.termination_status(model) == MOI.OPTIMAL
      @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
      @test JuMP.objective_value(model) ≈ 10.0
      @test JuMP.value(x) ≈ 1.0
      @test JuMP.value(y) ≈ 2.0
end


 # No need to add worker because a worker is already added while running testset "Using JuMP syntax by having two workers".
 # Otherwise we have to add a worker by using 'Distributed.addprocs(1)' and load the package on each worker
 # by '@everywhere using CharibdeOptim'.

@testset "Optimising by Charibde (A hybrid approach) using 2 workers" begin
      (global_min, minimisers, info) = charibde_min(X->((x,y)=X;x^3 + 2y + 5), IntervalBox(2..4, 2..3))
      @test global_min ⊆ 17.0 .. 17.01
      @test minimisers[1] ⊆ (2..2.001) × (2..2.001)

      (global_min, minimisers, info)= charibde_min(X->((x,y)=X;x^2 + y^2), IntervalBox(2..3, 3..4))
      @test global_min ⊆ 13 .. 13.01
      @test minimisers[1] ⊆ (2.0 .. 2.001) × (3..3.001)

end

using ModelingToolkit
using IntervalConstraintProgramming

@testset "Using Interval bound and contract algorithm for Constrained Optimisation" begin
      vars = ModelingToolkit.@variables x y
      C1 = Constraint(vars, x+y, -Inf..4)
      C2 = Constraint(vars, x+3y, -Inf..9)

      (maxima, maximisers) = ibc_maximise(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4),[C1, C2])
      @test maxima ⊆ -6.75 .. -6.74
      @test maximisers[1] ⊆ (2.92 .. 2.93) × (1.63 .. 1.64)
end
