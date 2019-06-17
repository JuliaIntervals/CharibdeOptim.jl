module CharibdeOptim

export ConstraintCond, charibde_min, charibde_max
export ibc_minimise, diffevol_minimise

using IntervalArithmetic
using Distributed
using IntervalOptimisation: HeapedVector, SortedVector, filter_elements!
using IntervalConstraintProgramming
using ModelingToolkit

import Base: invokelatest, push!

include("IBC.jl")
include("DiffrentialEvolution.jl")
include("ConstraintDifferentialEvolution.jl")
include("BoundEnsure.jl")
include("GenerateRandom.jl")


function charibde_min(f::Function, X::T) where{T}

    chnl1 = RemoteChannel(()->Channel{Tuple{T, Float64}}(10))                       #IBC recieve element from this channel
    chnl2 = RemoteChannel(()->Channel{Tuple{Vector{Float64}, Float64}}(10))        #DiffEvolution recieve element from this channel

    r1 = remotecall(diffevol_minimise, 2, f, X, ibc_chnl = chnl1, diffevol_chnl = chnl2)
    r2 = remotecall(ibc_minimise, 3, f, X, ibc_chnl = chnl1, diffevol_chnl = chnl2)

    return (fetch(r2))
end

function charibde_max(f::Function, X:: T) where{T}
    return charibde_min(x -> -f(x), X)
end

end
