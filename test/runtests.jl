using Test, JuMP

using Distributed
addprocs(2)
@everywhere using CharibdeOptim
@everywhere using IntervalArithmetic

@testset "Using Charibde for Constrained Optimsation" begin
      @everywhere using ModelingToolkit

      @everywhere vars = ModelingToolkit.@variables x y
      @everywhere C1 = constraint(vars, x+y, -Inf..4)
      @everywhere C2 = constraint(vars, x+3y, -Inf..9)
      @everywhere constraints = [C1, C2]

      prob = ConstrainedOptimisationProblem(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4), constraints)
      maxima, maximisers, info = maximise_by_charibde(prob)
      @test maxima ⊆ -8.01 .. -7.99
      @test maximisers[1] ⊆ (1.99 .. 2.01) × (1.99 .. 2.01)
end


@testset "Using Interval bound and contract algorithm for Constrained Optimisation" begin
      vars = ModelingToolkit.@variables x y
      C1 = constraint(vars, x+y, -Inf..4)
      C2 = constraint(vars, x+3y, -Inf..9)

      prob = ConstrainedOptimisationProblem(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4), [C1, C2])
      (maxima, maximisers, info) = maximise_by_ibc(prob)
      @test maxima ⊆ -8.01 .. -7.99
      @test maximisers[1] ⊆ (1.99 .. 2.01) × (1.99 .. 2.01)
end


@testset "Optimising by Interval Branch and Contract Algorithm" begin

      prob = OptimisationProblem(X->((x,y)= X;x^2 + y^2), IntervalBox(2..3, 3..4))

      (global_min, minimisers)= minimise_by_ibc(prob)
      @test global_min ⊆ 13 .. 13.01
      @test minimisers[1] ⊆ (2.0 .. 2.001) × (3.0 .. 3.001)
end

@testset "Optimising by Differencial Evolution" begin
      prob = OptimisationProblem(X->((x,y)= X;x^2 + y^2), IntervalBox(2..3, 3..4))
      (global_min, minimiser)= minimise_by_diffevol(prob, iterations = 40)
      @test global_min ⊆ 13 .. 13.01
      @test minimiser ⊆ (2.0 .. 2.001) × (3.0 .. 3.001)
end


@testset "Optimising by Charibde (A hybrid approach) using only one worker" begin

      prob = OptimisationProblem(X->((x,y)=X;x^2+y+1), IntervalBox(1..2, 2..3))
      (global_min, minimisers, info) = minimise_by_charibde(prob, workers = 1)
      @test global_min ⊆ 4.0 .. 4.01
      @test minimisers[1] ⊆ (1..1.001) × (2..2.001)

      prob = OptimisationProblem(X->((x,y)=X;x^2 + y^2), IntervalBox(2..3, 3..4))
      (global_min, minimisers, info)= minimise_by_charibde(prob, workers = 1)
      @test global_min ⊆ 13 .. 13.01
      @test minimisers[1] ⊆ (2.0 .. 2.001) × (3..3.001)
end


@testset "Optimising by Charibde (A hybrid approach) using 2 workers" begin

      prob = OptimisationProblem(X->((x,y)=X;x^3 + 2y + 5), IntervalBox(2..4, 2..3))
      (global_min, minimisers, info) = minimise_by_charibde(prob)
      @test global_min ⊆ 17.0 .. 17.01
      @test minimisers[1] ⊆ (2..2.001) × (2..2.001)

      prob = OptimisationProblem(X->((x,y)=X;x^2 + y^2), IntervalBox(2..3, 3..4))
      (global_min, minimisers, info)= minimise_by_charibde(prob)
      @test global_min ⊆ 13 .. 13.01
      @test minimisers[1] ⊆ (2.0 .. 2.001) × (3..3.001)

end

@testset "Optimising difficult problems using Charibde" begin
      f = X->((x1,x2,x3,x4,x5,x6,x7,x8,x9)=X;x1^2 + x2^2 + x3^4 - x4^7 - 200x5 - x6^5 - x7^9 + x8^5 - 8x9^3)
      X = IntervalBox(2..3, 3..4, 9..14, 2..3, 3..4, 9..14, 2..3, 3..4, 9..14)

      prob = OptimisationProblem(f, X)
      (global_min, minimisers, info) = minimise_by_charibde(prob)

      @test global_min ⊆ -575630 .. -575628
      @test minimisers[1] ⊆ (1.99 .. 2.01) × (2.99 .. 3.01) × (8.99 .. 9.01) × (2.99 .. 3.01) × (3.99 .. 4.01) × (13.99 .. 14.01) × (2.99 .. 3.01) × (2.99 .. 3.01) × (13.99 .. 14.01)

end

@testset "Optimising using JuMP syntax" begin
      include("JuMP_test.jl")
end
