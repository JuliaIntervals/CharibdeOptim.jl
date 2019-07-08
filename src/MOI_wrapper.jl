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

MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

Optimizer() = Optimizer(nothing, empty_nlp_data(), [], MOI.FEASIBILITY_SENSE, nothing)

MOI.get(::Optimizer, ::MOI.SolverName) = "CharibdeOptim"

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.LessThan{Float64}}) = true

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

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F <: Union{MOI.SingleVariable, MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction}}
    model.objective = func
    return
end

function MOI.is_empty(model::Optimizer)
    return isempty(model.variable_info) && model.nlp_data.evaluator isa EmptyNLPEvaluator && model.sense == MOI.FEASIBILITY_SENSE
end

function MOI.empty!(model::Optimizer)
    model.result = nothing
    empty!(model.variable_info)
    model.nlp_data = empty_nlp_data()
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
end

function MOI.add_constraint(model::Optimizer, var::MOI.SingleVariable, bound::MOI.Interval)
    vi = var.variable
    model.variable_info[vi.value].lower_bound = bound.lower
    model.variable_info[vi.value].upper_bound = bound.upper
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval}(vi.value)
end

function MOI.add_constraint(model::Optimizer, var::MOI.SingleVariable, hi::MOI.LessThan{Float64})
    vi = var.variable
    model.variable_info[vi.value].upper_bound = hi.upper
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(vi.value)
end

function MOI.add_constraint(model::Optimizer, var::MOI.SingleVariable, low::MOI.GreaterThan{Float64})
    vi = var.variable
    model.variable_info[vi.value].lower_bound = low.lower
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}(vi.value)
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


function substitute_variables(expr::Expr)
    if expr.head == :ref && length(expr.args) == 2 && expr.args[1] == :x
        return :(x[$(expr.args[2].value)])
    else
        for (index, arg) in enumerate(expr.args)
            expr.args[index] = substitute_variables(arg)
        end
    end
    return expr
end
substitute_variables(arg) = arg

function eval_objective(model::Optimizer, eval_expr::Union{Function, Nothing}, x)
    if eval_expr !== nothing
        return invokelatest(eval_expr, x)
    elseif model.objective !== nothing
        return eval_function(model.objective, x)
    else
        return 0.0
    end
end

function MOI.optimize!(model::Optimizer; debug =true, chnl1 = nothing, chnl2 = nothing)

    eval_expr = nothing
    if model.nlp_data.has_objective
        MOI.initialize(model.nlp_data.evaluator, [:ExprGraph])
        expr = MOI.objective_expr(model.nlp_data.evaluator)
        expr = substitute_variables(expr)
        eval_expr = eval(:(x -> $(expr)))
    end

    obj_func(x...) = eval_objective(model, eval_expr, x)

    X = [Interval(var.lower_bound, var.upper_bound) for var in model.variable_info]
    search_space = IntervalBox(X...)

    if myid() == 2
        if model.sense == MOI.MIN_SENSE
            diffevol_minimise(obj_func, search_space, chnl1, chnl2)
        elseif model.sense == MOI.MAX_SENSE
            diffevol_maximise(obj_func, search_space, chnl1, chnl2)
        else
            error("Min or Max Sense is not set")
        end
    elseif myid() == 1
        chnl1 = RemoteChannel(()->Channel{Tuple{IntervalBox, Float64}}(1))
        chnl2 = RemoteChannel(()->Channel{Tuple{Vector{Float64}, Float64}}(1))

        remotecall(MOI.optimize!, 2, model, debug = debug, chnl1 = chnl1, chnl2 = chnl2)

        if model.sense == MOI.MIN_SENSE
            model.result = ibc_minimise(obj_func, search_space, debug = debug, ibc_chnl = chnl1, diffevol_chnl = chnl2)
        elseif model.sense == MOI.MAX_SENSE
            model.result = ibc_maximise(obj_func, search_space, debug = debug, ibc_chnl = chnl1, diffevol_chnl = chnl2)
        else
            error("Min or Max Sense is not set")
        end
    end
end


function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.result == nothing
        return MOI.OPTIMIZE_NOT_CALLED
    else
        return MOI.OPTIMAL
    end
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    if model.result === nothing
        error("ObjectiveValue not available.")
    end
    return model.result[1].lo
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    return (model.result !== nothing) ? 1 : 0
end

function MOI.get(model::Optimizer, ::MOI.PrimalStatus)
    if model.result === nothing
        return MOI.NO_SOLUTION
    else
        return MOI.FEASIBLE_POINT
    end
end

function MOI.get(model::Optimizer, ::MOI.DualStatus)
        return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if model.result === nothing
        error("VariablePrimal not available.")
    end
    return model.result[2][1][vi.value].lo
end
