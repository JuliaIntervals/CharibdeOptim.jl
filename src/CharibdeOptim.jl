module CharibdeOptim

export Constraint, charibde_min, charibde_max, ibc_maximise
export ibc_minimise, diffevol_minimise

using IntervalArithmetic
using Distributed
using IntervalOptimisation: SortedVector, filter_elements!
using IntervalConstraintProgramming
using ModelingToolkit
using MathOptInterface
using StaticArrays

import Base: invokelatest, push!

struct Constraint{T}
   bound::Interval{T}
   C::BasicContractor
end

function Constraint(vars, constraint_expr::Operation, bound::Interval{T}; epsilon = 1e-4) where{T}
   C = BasicContractor(vars, constraint_expr)
   if diam(bound) == 0.0
       bound = Interval(bound.lo - epsilon, bound.hi + epsilon )
   end
   return Constraint{T}(bound, C)
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
include("MaxDist.jl")
include("BoundEnsure.jl")
include("GenerateRandom.jl")


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
(global_min, minimisers, info) = charibde_min(f, A)
(global_max, maximisers, info) = charibde_max(f, A)

charibde_min/charibde_max find the global minimum/maximum value of the function in given search space by using two algorithms Interval Bound & Contract(IBC) and Differential Evolution inparallel
TODO: Fix the parallel implementation of charibde_min/max for Constrained Optimisation
```
"""
function charibde_min(f::Function, X::IntervalBox{N,T}; workers = 2, tol = 1e-3, np = N*10, debug = false) where{N,T}

    worker_ids = Distributed.workers()
    if workers > 1
        if worker_ids[1] == 1
            error("Not enough workers available: Add one worker and load the package on it")
        end

        (chnl1, chnl2) = create_channels(X, workers)                     #IBC recieve element from chnl1 and DiffEvolution from chnl2

        remotecall(diffevol_minimise, worker_ids[1], f, X, chnl1, chnl2, np = np)
        return ibc_minimise(f, X, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol = tol, debug = debug)
    else
        (chnl1, chnl2) = create_channels(X, workers)

        @async diffevol_minimise(f, X, chnl1, chnl2, np = np)
        return ibc_minimise(f, X, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol = tol, debug = debug)
    end
end


function charibde_min(f::Function, X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}; workers = 2, tol = 1e-3, np = N*10, debug = false) where{N,T}

    worker_ids = Distributed.workers()
    if workers > 1
        if worker_ids[1] == 1
            error("Not enough workers available: Add one worker and load the package on it")
        end

        (chnl1, chnl2) = create_channels(X, workers)

        remotecall(diffevol_minimise, worker_ids[1], f, X, constraints, chnl1, chnl2, np = np)
        return ibc_minimise(f, X, constraints, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol = tol, debug = debug)
    else
        (chnl1, chnl2) = create_channels(X, workers)

        @async diffevol_minimise(f, X, constraints, chnl1, chnl2, np = np)
        return ibc_minimise(f, X, constraints, ibc_chnl = chnl1, diffevol_chnl = chnl2, tol= tol, debug = debug)
    end
end


function charibde_max(f::Function, X::IntervalBox{N,T}; workers = 2, tol = 1e-3, np = N*10, debug = false) where{N,T}
    bound, minimizers, info = charibde_min(x -> -f(x), X, workers = workers, tol = tol, np = np, debug = debug)
    return -bound, minimizers, info
end

function charibde_max(f::Function, X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}; workers = 2, tol = 1e-3, np = N*10, debug = false) where{N,T}
    bound, minimizers, info = charibde_min(x -> -f(x), X, constraints, workers = workers, tol = tol, np = np, debug = debug)
    return -bound, minimizers, info
end

include("MOI_wrapper.jl")

end
