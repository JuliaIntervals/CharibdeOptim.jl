using Charibde
using IntervalArithmetic
using Test

@testset "Charibde tests " begin

      @testset "Optimising by Differential Evolution Algorithm" begin
            (global_min, X_best) = DiffEvolution_min((x,y)->x^2+y+1, IntervalBox(1..2, 2..3))
            @test global_min ⊆ 4.0 .. 4.001
            @test length(X_best) == 2
            @test X_best[1] ⊆ 1.0 .. 1.001
            @test X_best[2] ⊆ 2.0 .. 2.001
      end

      @testset "Optimising by Interval Bound and Contract Algorithm" begin
            (global_min, minimisers)= IBC_min( X -> ( (x,y) = X; x^2 + y^2 ) , IntervalBox(2..3, 3..4))
            @test global_min ⊆ 13 .. 13.01
            @test length(minimisers) == 2
            @test minimisers[1] ⊆ (2.0 .. 2.001) × (3.0 .. 3.001)
      end

end
