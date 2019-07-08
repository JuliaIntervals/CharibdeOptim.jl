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

include("IBC.jl")
include("DiffrentialEvolution.jl")
include("ConstraintDifferentialEvolution.jl")
include("BoundEnsure.jl")
include("GenerateRandom.jl")
include("MaxDist.jl")


function charibde_min(f::Function, X::IntervalBox{N,T}) where{N,T}

    len = length(X)
    chnl1 = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(1))                       #IBC recieve element from this channel
    chnl2 = RemoteChannel(()->Channel{Tuple{MArray{Tuple{N},Float64,1,N},Float64}}(1))        #DiffEvolution recieve element from this channel

    r1 = remotecall(ibc_minimise, 2, f, X, ibc_chnl = chnl1, diffevol_chnl = chnl2)
    r2 = remotecall(diffevol_minimise, 3, f, X, chnl1, chnl2)
    return fetch(r1)
end

function charibde_max(f::Function, X:: T) where{T}
    return charibde_min(x -> -f(x), X)
end

include("MOI_wrapper.jl")

end
