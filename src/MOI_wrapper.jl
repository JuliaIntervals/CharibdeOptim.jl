import MathOptInterface

const MOI = MathOptInterface


mutable struct VariableInfo
    lower_bound::Float64
    upper_bound::Float64
end

VariableInfo() = VariableInfo(-Inf, Inf)

mutable struct Optimizer <: MOI.AbstractOptimizer
    result
    variable_info::Vector{VariableInfo}
    sense::MOI.OptimizationSense
    objective::Union{MOI.SingleVariable,MOI.ScalarAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64},Nothing}
    variable_bound::Vector{Tuple{MOI.SingleVariable, MOI.Interval}}
end

Optimizer() = Optimizer(nothing, [], MOI.FEASIBILITY_SENSE, nothing, [])

MOI.get(::Optimizer, ::MOI.SolverName) = "CharibdeOptim"

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)

function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i in 1:length(model.variable_info)]
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    model.sense = sense
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction, func::Union{MOI.SingleVariable, MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction})
    model.objective = func
    return
end


function MOI.add_constraint(model::Optimizer, var::MOI.SingleVariable, bound::MOI.Interval)
    vi = var.variable
    model.variable_info[vi.value].lower_bound = bound.lower
    model.variable_info[vi.value].upper_bound = bound.upper
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(vi.value)
end

function MOI.add_variable(model::Optimizer)
    push!(model.variable_info, VariableInfo())
    return MOI.VariableIndex(length(model.variable_info))
end
function MOI.add_variables(model::Optimizer, n::Int)
    return [MOI.add_variable(model) for i in 1:n]
end

function eval_function(var::MOI.SingleVariable, x)
    return x[var.variable.value]
end


function MOI.optimize!(model::Optimizer)

    function obj_func(X...)
       f2 = x -> MOI.Utilities.evalvariables(vi -> x[vi.value], model.objective)
       return f2(X)
    end

    X = Interval[]

    for var in model.variable_info
        lower = var.lower_bound
        upper = var.upper_bound
        push!(X, lower..upper)
    end
    if model.sense == MOI.MIN_SENSE
        model.result = ibc_minimise(obj_func, IntervalBox(X...))
    elseif model.sense == MOI.MAX_SENSE
        model.result = ibc_maximise(obj_func, IntervalBox(X...))
    else
        error("Min or Max Sense is not set")
    end
end
