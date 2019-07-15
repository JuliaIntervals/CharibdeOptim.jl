module CharibdeOptim

export ConstraintCond, charibde_min, charibde_max
export ibc_minimise, diffevol_minimise

using IntervalArithmetic
using Distributed
using IntervalOptimisation: SortedVector, filter_elements!
using IntervalConstraintProgramming
using ModelingToolkit
using MathOptInterface
using StaticArrays

import Base: invokelatest, push!

struct ConstraintCond{T}
   bound::Interval{T}
   C::BasicContractor
end

function ConstraintCond(vars::Vector{Operation}, constraint_expr::Operation, bound::Interval{T}) where{T}
   C = BasicContrator(vars, constraint_expr)
   return ConstraintCond{T}(bound, C)
end

mutable struct Information
    de_to_ibc::Int
    ibc_to_de::Int
    iterations::Int
end


include("IBC.jl")
include("DiffrentialEvolution.jl")
include("ConstraintDifferentialEvolution.jl")
include("BoundEnsure.jl")
include("GenerateRandom.jl")
include("MaxDist.jl")



function charibde_min(f::Function, X::IntervalBox{N,T}; workers = 2, debug = false) where{N,T}

    worker_ids = Distributed.workers()
    if workers > 1
        if worker_ids[1] == 1
            error("Not enough workers available: Add one worker and load the package on it")
        end

        chnl1 = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(1))                       #IBC recieve element from this channel
        chnl2 = RemoteChannel(()->Channel{Tuple{SArray{Tuple{N},Float64,1,N},Float64}}(1))            #DiffEvolution recieve element from this channel

        remotecall(diffevol_minimise, worker_ids[1], f, X, chnl1, chnl2)
        return ibc_minimise(f, X, debug = debug, ibc_chnl = chnl1, diffevol_chnl = chnl2)
    else
        chnl1 = Channel{Tuple{IntervalBox{N,T}, Float64}}(1)
        chnl2 = Channel{Tuple{SArray{Tuple{N},Float64,1,N},Float64}}(1)

        @async diffevol_minimise(f, X, chnl1, chnl2)
        return ibc_minimise(f, X, debug = debug, ibc_chnl = chnl1, diffevol_chnl = chnl2)
    end
end

function charibde_min(f::Function, X::IntervalBox{N,T}, C::Vector{ConstraintCond{T}}; workers = 2, debug = false) where{N,T}

    worker_ids = Distributed.workers()
    if workers > 1
        if worker_ids[1] == 1
            error("Not enough workers available: Add one worker and load the package on it")
        end

        chnl1 = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(1))                       #IBC recieve element from this channel
        chnl2 = RemoteChannel(()->Channel{Tuple{SArray{Tuple{N},Float64,1,N},Float64}}(1))            #DiffEvolution recieve element from this channel

        remotecall(diffevol_minimise, worker_ids[1], f, X, C, chnl1, chnl2)
        return ibc_minimise(f, X, C, debug = debug, ibc_chnl = chnl1, diffevol_chnl = chnl2)
    else
        chnl1 = Channel{Tuple{IntervalBox{N,T}, Float64}}(1)
        chnl2 = Channel{Tuple{SArray{Tuple{N},Float64,1,N},Float64}}(1)

        @async diffevol_minimise(f, X, C, chnl1, chnl2)
        return ibc_minimise(f, X, C, debug = debug, ibc_chnl = chnl1, diffevol_chnl = chnl2)
    end
end


function charibde_max(f::Function, X::IntervalBox{N,T}; workers = 2, debug = false) where{N,T}
    return charibde_min(x -> -f(x), X, workers, debug)
end

function charibde_max(f::Function, X::IntervalBox{N,T}, C::Vector{ConstraintCond{T}}; workers = 2, debug = false) where{N,T}
    return charibde_min(x -> -f(x), X, C, workers, debug)
end

include("MOI_wrapper.jl")

end
