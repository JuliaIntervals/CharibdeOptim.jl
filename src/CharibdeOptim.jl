module CharibdeOptim

export  minimise_by_charibde, maximise_by_charibde
export minimise_by_ibc, minimise_by_diffevol, maximise_by_ibc, maximise_by_diffevol
export constraint, OptimisationProblem, ConstrainedOptimisationProblem

using IntervalArithmetic
using Distributed
using IntervalOptimisation: SortedVector, filter_elements!
using IntervalConstraintProgramming
using ModelingToolkit
using MathOptInterface
using StaticArrays
using MacroTools

import Base: invokelatest, push!
import ForwardDiff: gradient

struct OptimisationProblem{N, T}
    f::Function
    X::IntervalBox{N, T}
end


struct Constraint{T}
   bound::Interval{T}
   C::Contractor
end

struct ConstrainedOptimisationProblem{N, T}
    f::Function
    X::IntervalBox{N, T}
    constraints::Vector{Constraint{T}}
end
"""Usage:
```
For Unconstrained Optimsation:
  f = X->((x,y)=X;x^3 + 2y + 5)
  A = IntervalBox(2..4, 2..3)
  prob = OptimisationProblem(f, A)

  (global_min, minimisers, info) = minimise_by_ibc(prob)
  (global_max, maximisers, info) = maximise_by_ibc(prob)
  (global_min, minimisers) =  minimise_by_diffevol(prob)
  (global_max, maximisers) = maximiser_by_diffevol(prob)

For Constrained Optimisation:
  f = X->((x,y)=X;-(x-4)^2-(y-4)^2)
  A = IntervalBox(-4..4, -4..4)

  vars = ModelingToolkit.@variables x y
  C1 = Constraint(vars, x+y, -Inf..4)
  C2 = Constraint(vars, x+3y, -Inf..9)

  prob = ConstrainedOptimisation(f, A, [C1, C2])

  (global_min, minimisers, info) = minimise_by_ibc(prob)
  (global_max, maximisers, info) = maximise_by_diffevol(prob)
  (global_min, minimisers) =  minimise_by_diffevol(prob)
  (global_max, maximisers) = maximiser_by_diffevol(prob)

minimise_by_ibc/maximise_by_ibc, minimise_by_diffevol/maximise_by_diffevol find the global minimum/maximum value of the function in given search space by using Interval Bound & Contract(IBC)  and Differential Evolution algorithm
```
"""
function minimise_by_ibc(prob::OptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, T}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), structure = SortedVector, debug = false, tol=1e-6, ibc_ind = true) where{N, T}

    vars = [Variable(Symbol("x",i))() for i in 1:length(prob.X)]
    g(x...) = prob.f(x)
    C = BasicContractor(vars, g)
    invokelatest(ibc_minimise, prob.f, prob.X, C, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, structure = structure, debug = debug, tol=tol, ibc_ind = ibc_ind)
end

function maximise_by_ibc(prob::OptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, T}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), structure = SortedVector, debug = false, tol=1e-6, ibc_ind = true) where{N, T}
    vars = [Variable(Symbol("x",i))() for i in 1:length(prob.X)]
    g(x...) = prob.f(x)
    C = BasicContractor(vars, g)
    invokelatest(ibc_maximise, prob.f, prob.X, C, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, structure = structure, debug = debug, tol=tol, ibc_ind = ibc_ind)
end

function minimise_by_ibc(prob::ConstrainedOptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, T}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), structure = SortedVector, debug = false, tol=1e-6, ibc_ind = true) where{N, T}
    vars = [Variable(Symbol("x",i))() for i in 1:length(prob.X)]
    g(x...) = prob.f(x)
    C = BasicContractor(vars, g)
    invokelatest(ibc_minimise, prob.f, prob.X, prob.constraints, C, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, structure = structure, debug = debug, tol=tol, ibc_ind = ibc_ind)
end

function maximise_by_ibc(prob::ConstrainedOptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, T}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), structure = SortedVector, debug = false, tol=1e-6, ibc_ind = true) where{N, T}
    vars = [Variable(Symbol("x",i))() for i in 1:length(prob.X)]
    g(x...) = prob.f(x)
    C = BasicContractor(vars, g)
    invokelatest(ibc_maximise, prob.f, prob.X, prob.constraints, C, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, structure = structure, debug = debug, tol=tol, ibc_ind = ibc_ind)
end


minimise_by_diffevol(prob::OptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), np = 10*N, debug = false, de_ind = true, iterations = 20) where{N, T} = diffevol_minimise(prob.f, prob.X, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, np = np, debug = debug, de_ind = de_ind, iterations = iterations)

maximise_by_diffevol(prob::OptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), np = 10*N, debug = false, de_ind = true, iterations = 20) where{N, T} = diffevol_maximise(prob.f, prob.X, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, np = np, debug = debug, de_ind = de_ind, iterations = iterations)

minimise_by_diffevol(prob::ConstrainedOptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), np = 10*N, debug = false, de_ind = true, iterations = 20) where{N, T} = diffevol_minimise(prob.f, prob.X, prob.constraints, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, np = np, debug = debug, de_ind = de_ind, iterations = iterations)

maximise_by_diffevol(prob::ConstrainedOptimisationProblem{N, T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), np = 10*N, debug = false, de_ind = true, iterations = 20) where{N, T} = diffevol_maximise(prob.f, prob.X, prob.constraints, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, np = np, debug = debug, de_ind = de_ind, iterations = iterations)


function constraint(vars, constraint_expr::Operation, bound::Interval{T}; epsilon = 1e-4) where{T}
   C = Contractor(vars, constraint_expr)
   if diam(bound) == 0.0
       bound = Interval(bound.lo - epsilon, bound.hi + epsilon)
   end
   Constraint{T}(bound, C)
end

mutable struct Information
    de_to_ibc::Int
    ibc_to_de::Int
    iterations::Int
end


include("IBC.jl")
include("DifferentialEvolution.jl")
include("ConstrainedDifferentialEvolution.jl")
include("ConstrainedIBC.jl")
include("BoundEnsure.jl")
include("GenerateRandom.jl")
include("SecondOrder.jl")


function create_channels(X::IntervalBox{N, T}, workers::Int64) where{N, T}
    if workers > 1
        return (RemoteChannel(()->Channel{Tuple{IntervalBox{N, T}, T}}(1)), RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(1)))
    else
        return (Channel{Tuple{IntervalBox{N, T}, T}}(1), Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(1))
    end
end

"""Usage:
```
f = X->((x,y)=X;x^3 + 2y + 5)
A = IntervalBox(2..4, 2..3)
prob = OptimisationProblem(f, A)
(global_min, minimisers, info) = minimise_by_charibde(prob)
(global_max, maximisers, info) = maximise_by_charibde(prob)

charibde_min/charibde_max find the global minimum/maximum value of the function in given search space by using two algorithms Interval Bound & Contract(IBC) and Differential Evolution inparallel
TODO: Fix the parallel implementation of charibde_min/max for Constrained Optimisation
```
"""
function minimise_by_charibde(prob::OptimisationProblem{N, T}; workers = 2, tol = 1e-6, np = N*10, debug = false) where{N,T}

    worker_ids = Distributed.workers()
    if workers > 1
        if worker_ids[1] == 1
            error("Not enough workers available: Add one worker and load the package on it")
        end

        (chnl1, chnl2) = create_channels(prob.X, workers)                     #IBC recieve element from chnl1 and DiffEvolution from chnl2

        r1 = remotecall(minimise_by_diffevol, worker_ids[1], prob, chnl1, chnl2, np = np, de_ind = false)
        r2 = remotecall(minimise_by_ibc, worker_ids[2], prob, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol = tol, debug = debug, ibc_ind = false)
        return fetch(r2)
    else
        (chnl1, chnl2) = create_channels(prob.X, workers)

        r1 = @async minimise_by_diffevol(prob, chnl1, chnl2, np = np, de_ind = false)
        r2 = @async minimise_by_ibc(prob, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol = tol, debug = debug, ibc_ind = false)
        return fetch(r2)
    end
end


function minimise_by_charibde(prob::ConstrainedOptimisationProblem{N, T}; workers = 2, tol = 1e-6, np = N*10, debug = false) where{N,T}

    worker_ids = Distributed.workers()
    if workers > 1
        if nprocs() < 3
            error("Not enough workers available: Session should have atleast 2 workers more")
        end

        (chnl1, chnl2) = create_channels(prob.X, workers)

        r1 = remotecall(minimise_by_diffevol, worker_ids[1], prob, ibc_chnl = chnl1, diffevol_chnl = chnl2, np = np, de_ind = false)
        r2 = remotecall(minimise_by_ibc, worker_ids[2], prob, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol = tol, debug = debug, ibc_ind = false)
        return fetch(r2)
    else
        (chnl1, chnl2) = create_channels(X, workers)

        r1 = @async minimise_by_diffevol(prob, ibc_chnl = chnl1, diffevol_chnl = chnl2, np = np, de_ind = false)
        r2 = @async minimise_by_ibc(prob, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol= tol, debug = debug, ibc_ind = false)
        return fetch(r2)
    end
end


function maximise_by_charibde(prob::OptimisationProblem{N, T}; workers = 2, tol = 1e-6, np = N*10, debug = false) where{N,T}
    new_prob = OptimisationProblem(x -> -prob.f(x), prob.X)
    bound, minimizers, info = minimise_by_charibde(new_prob, workers = workers, tol = tol, np = np, debug = debug)
    return -bound, minimizers, info
end

function maximise_by_charibde(prob::ConstrainedOptimisationProblem{N, T}; workers = 2, tol = 1e-6, np = N*10, debug = false) where{N,T}
    new_prob = ConstrainedOptimisationProblem(x-> -prob.f(x), prob.X, prob.constraints)
    bound, minimizers, info = minimise_by_charibde(new_prob, workers = workers, tol = tol, np = np, debug = debug)
    return -bound, minimizers, info
end

include("MOI_wrapper.jl")

end
