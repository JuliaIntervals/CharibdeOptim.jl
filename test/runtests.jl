using Test, JuMP

using Distributed
addprocs(2)
@everywhere using CharibdeOptim
@everywhere using IntervalArithmetic

@testset "Using Charibde for Constrained Optimsation" begin
      @everywhere using IntervalArithmetic
      @everywhere using ModelingToolkit

      @everywhere vars = ModelingToolkit.@variables x y
      @everywhere C1 = constraint(vars, x+y, -Inf..4)
      @everywhere C2 = constraint(vars, x+3y, -Inf..9)
      @everywhere constraints = [C1, C2]
      (maxima, maximisers, info) = charibde_max(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4), constraints)

      @test_skip maxima ⊆ -8.1 .. -7.9
      @test_skip maximisers[1] ⊆ (1.9 .. 2.1) × (1.9 .. 2.1)
end

@testset "Using Interval bound and contract algorithm for Constrained Optimisation" begin
      vars = ModelingToolkit.@variables x y
      C1 = constraint(vars, x+y, -Inf..4)
      C2 = constraint(vars, x+3y, -Inf..9)

      (maxima, maximisers, info) = ibc_maximise(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4),[C1, C2])
      @test_skip maxima ⊆ -8.1 .. -7.9
      @test_skip maximisers[1] ⊆ (1.9 .. 2.1) × (1.9 .. 2.1)
end


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


@testset "Using JuMP syntax by using only one worker " begin           #for using two workers just dont pass 'workers' arguments as its value is set to 2
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


 # No need to add worker because a worker is already added while running testset "Using Charibde for Constrained Optimsation".
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
