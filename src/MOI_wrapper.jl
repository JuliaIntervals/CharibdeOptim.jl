import MathOptInterface

const MOI = MathOptInterface


mutable struct VariableInfo
    lower_bound::Float64
    upper_bound::Float64
end

VariableInfo() = VariableInfo(-Inf, Inf)

mutable struct Optimizer <: MOI.AbstractOptimizer
    result
    nlp_data::MOI.NLPBlockData
    variable_info::Vector{VariableInfo}
    sense::MOI.OptimizationSense
    objective::Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}, Nothing}
end

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

Optimizer() = Optimizer(nothing, empty_nlp_data(), [], MOI.FEASIBILITY_SENSE, nothing)

MOI.get(::Optimizer, ::MOI.SolverName) = "CharibdeOptim"

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true


MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)

function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i in 1:length(model.variable_info)]
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    model.sense = sense
    return
end

function MOI.set(model::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    model.nlp_data = nlp_data
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction,func::Union{MOI.SingleVariable, MOI.ScalarAffineFunction,
                             MOI.ScalarQuadraticFunction})
    model.objective = func
    return
end

function MOI.is_empty(model::Optimizer)
    return isempty(model.variable_info) && model.sense == MOI.FEASIBILITY_SENSE
end

function MOI.empty!(model::Optimizer)
    model.result = nothing
    empty!(model.variable_info)
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
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

function eval_function(aff::MOI.ScalarAffineFunction, x)
    function_value = aff.constant
    for term in aff.terms
        # Note the implicit assumtion that VariableIndex values match up with
        # x indices. This is valid because in this wrapper ListOfVariableIndices
        # is always [1, ..., NumberOfVariables].
        function_value += term.coefficient*x[term.variable_index.value]
    end
    return function_value
end

function eval_function(quad::MOI.ScalarQuadraticFunction, x)
    function_value = quad.constant
    for term in quad.affine_terms
        function_value += term.coefficient*x[term.variable_index.value]
    end
    for term in quad.quadratic_terms
        row_idx = term.variable_index_1
        col_idx = term.variable_index_2
        coefficient = term.coefficient
        if row_idx == col_idx
            function_value += 0.5*coefficient*x[row_idx.value]*x[col_idx.value]
        else
            function_value += coefficient*x[row_idx.value]*x[col_idx.value]
        end
    end
    return function_value
end

function eval_objective(model::Optimizer, x)
    # The order of the conditions is important. NLP objectives override regular
    # objectives.
    if model.nlp_data.has_objective
        return MOI.eval_objective(model.nlp_data.evaluator, x)
    elseif model.objective !== nothing
        return eval_function(model.objective, x)
    else
        # No objective function set. This could happen with FEASIBILITY_SENSE.
        return 0.0
    end
end


function MOI.optimize!(model::Optimizer)

    obj_func(X...) = eval_objective(model, X)
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
    return
end
