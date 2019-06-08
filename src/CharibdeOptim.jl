module CharibdeOptim

export ConstraintCond, charibde_min
export ibc_minimise, ibc_maximise, diffevol_minimise, diffevol_maximise

using IntervalArithmetic
using Distributed
using IntervalOptimisation: HeapedVector, SortedVector, filter_elements!
using IntervalConstraintProgramming
using ModelingToolkit

import Base.push!

include("IBC.jl")
include("DiffrentialEvolution.jl")
include("ConstraintDifferentialEvolution.jl")
include("BoundEnsure.jl")
include("GenerateRandom.jl")

ibc_chnl = RemoteChannel(()->Channel(10))        #IBC recieve element from this channel from DiffEvolution
diffevol_chnl = RemoteChannel(()->Channel(10))    #DiffEvolution recieve element from this channel


function charibde_min(f::Function, X:: T) where{T}

    ibc_chnl = RemoteChannel(()->Channel(10))        #IBC recieve element from this channel from DiffEvolution
    diffevol_chnl = RemoteChannel(()->Channel(10))    #DiffEvolution recieve element from this channel

    r1 = remotecall(diffevol_minimise, 2, f, X, ibc_chnl, diffevol_chnl)
    r2 = remotecall(ibc_minimise, 3, f, X, ibc_chnl, diffevol_chnl)

    R2 = fetch(r2)
    R1 = fetch(r1)
    return (R2, R1)
end

end
