include("problems.jl")
using BenchmarkTools

using Distributed
addprocs(2)
@everywhere using CharibdeOptim
@everywhere using IntervalArithmetic


const SUITE = BenchmarkGroup()


for solver in (charibde_min, ibc_minimise)
    S = SUITE[string(solver)]
    s = S["Normal problems"]
    for prblm in (prblm_1, prblm_2, prblm_3)
        s[string(prblm)] = @benchmarkable($(solver), $(prblm.func), $(prblm.domain))
