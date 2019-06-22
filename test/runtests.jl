using Distributed
addprocs(2)
@everywhere using CharibdeOptim
@everywhere using IntervalArithmetic
using Test

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

end
