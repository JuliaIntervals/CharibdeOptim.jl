function hc4(X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}, tol=1e-5) where{N, T}
    n = length(constraints)
    while true
        X_temp = X
        for i in 1:n
            X = invokelatest(constraints[i].C,  constraints[i].bound, X)
        end
        if isempty(X) || sum(dist.(X, X_temp)) < tol
            break
        end
    end

    new_constraints = Constraint{T}[]

    for i in 1:n
        if !(invokelatest(constraints[i].C, X) ⊆ constraints[i].bound)
            push!(new_constraints, constraints[i])
        end
    end
    return new_constraints, X
end

function contraction(f::Function, C, global_min::Float64, X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}, tol=1e-5) where {N, T}
    lb = -Inf
    while true
        X_temp = X
        X = invokelatest(C, -Inf..global_min, X)
        lb = inf(f(X))
        constraints, X = hc4(X, constraints)

        if isempty(X) || sum(dist.(X, X_temp)) < tol
            break
        end
    end

    return lb, X, constraints
end

function generate_random_feasible_point(X::IntervalBox{N, T}, constraints::Vector{Constraints{T}}) where{N, T}
    point = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:length(X)]

    for j in 1:length(constraints)
        if !(invokelatest(constraints[j].C, point) ⊆ constraints[j].bound)
            point = generate_random_feasible_point(X, constraints)
            break
        end
    end

    return point
end



function ibc_minimise(f::Function , X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)), diffevol_chnl = Nothing, structure = SortedVector, debug = false, tol=1e-6) where{N, T}

    vars = [Variable(Symbol("x",i))() for i in 1:length(X)]
    g(x...) = f(x)
    C = BasicContractor(vars, g)

    working = structure([(X, inf(f(X)))], x->x[2]) # list of boxes with corresponding lower bound, arranged according to selected structure :

    minimizers = IntervalBox{N,T}[]
    global_min = ∞  # upper bound

    info = Information(0, 0, 0)
    num_bisections = 0

    while !isempty(working)

        if isready(ibc_chnl)
            from_diff = take!(ibc_chnl)     # Receiving best individual from ibc_minimise
            if debug
                println("Box recevied from DifferentialEvolution: ", from_diff)
            end
            global_min = min(from_diff[2], global_min)
            info.de_to_ibc = info.de_to_ibc + 1
        end
        (X, X_min) = popfirst!(working)

        if debug
            println("New search-space : ", X)
        end


        X_min, X, constraints = contraction(f, C, global_min, X, constraints)

        if debug
            println("Contracted search_space: ", X)
        end

        if X_min > global_min    # X_min == inf(f(X))
            continue
        end

        # find candidate for upper bound of global minimum by just evaluating a point in the interval:
        m = sup(f(Interval.(generate_random_feasible_point(X, constraints))))   # evaluate at feasible point

        if m < global_min
            global_min = m
            x_best = SVector(mid(X))
            if diffevol_chnl != Nothing
                if debug
                    println("Box send to DifferentialEvolution: ", x_best )
                end
                put!(diffevol_chnl, (x_best, global_min))  # sending best individual to diffevol
                info.ibc_to_de = info.ibc_to_de + 1
            end
        end


        filter_elements!(working , (X, global_min) )   # Remove all boxes whose lower bound is greater than the current one:

        if diam(X) < tol
            push!(minimizers, X)
        else
            X1, X2 = bisect(X)
            push!( working, (X1, inf(f(X1))) )
            push!( working, (X2, inf(f(X2))) )
            num_bisections += 1
        end
        info.iterations= info.iterations + 1
    end

    if debug
        println("IBC search is ended")
    end

    if diffevol_chnl != Nothing
        if isready(diffevol_chnl)
            take!(diffevol_chnl)
            put!(diffevol_chnl,(SVector(mid(X)), Inf) )
        else
            put!(diffevol_chnl,(SVector(mid(X)), Inf) )
        end
        if debug
            println("DifferentialEvolution is also terminated")
        end
    end

    if isready(ibc_chnl)
        take!(ibc_chnl)
    end

    #lower_bound = minimum(inf.(f.(minimizers)))

    return Interval(global_min, global_min), minimizers, info

end


function ibc_maximise(f::Function, X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}; debug = false, tol = 1e-6) where{N, T}
    bound, minimizer, info = ibc_minimise(x -> -f(x), X, constraints, debug = debug, tol = tol)
    return -bound, minimizer, info
end
