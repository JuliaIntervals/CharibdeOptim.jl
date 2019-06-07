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


function charibde_min(f::Function, X:: T) where{T}
    r1 = remotecall(diffevol_minimise, 2, f, X)
    r2 = remotecall(ibc_minimise, 3, f, X)
    return (fetch(r1), fetch(r2))
end

end
