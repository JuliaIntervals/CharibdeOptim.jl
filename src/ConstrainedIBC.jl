function hc4(X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}, tol=1e-5) where{N, T}
    n = length(constraints)
    while true
        X_temp = X
        for i in 1:n
            X = constraints[i].C(constraints[i].bound, X)
        end
        if isempty(X) || sum(dist.(X, X_temp)) < tol
            break
        end
    end

    #new_constraints = Constraint{T}[]

    #for i in 1:n
    #    if !(invokelatest(constraints[i].C, X) ⊆ constraints[i].bound)
    #        push!(new_constraints, constraints[i])
    #    end
    #end
    #return new_constraints, X
    return constraints, X
end

function contraction(f::Function, C, global_min::Float64, X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}, tol=1e-5) where {N, T}
    lb = -Inf
    while true
        X_temp = X
        X = C(-Inf..global_min, X)
        lb = inf(f(X))
        constraints, X = hc4(X, constraints)

        if isempty(X) || sum(dist.(X, X_temp)) < tol
            break
        end
    end

    return lb, X, constraints
end

function generate_random_feasible_point(X::IntervalBox{N, T}, constraints::Vector{Constraint{T}}) where{N, T}
    for i in 1:30
        point = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:length(X)]      # discover a random point in interval box X
        for j in 1:length(constraints)
            if !(constraints[j].C(point) ⊆ constraints[j].bound)
                break
            end
            if j == length(constraints)
                return (true, point)    # return a feasible point
            end
        end
        if i == 30
            return (false, point)    # returns a infeasible point
        end
    end

end

function convex_hull(vec::Vector{IntervalBox{N,T}}) where{N, T}

    x_big = Interval{T}[]
    num_variables = N
    num_boxes = length(vec)

    for i in 1:num_variables
        lower_bound = Inf
        upper_bound = -Inf
        for j in 1:num_boxes
            if vec[j][i].lo < lower_bound
                lower_bound = vec[j][i].lo
            end
            if upper_bound < vec[j][i].hi
                upper_bound = vec[j][i].hi
            end
        end
        append!(x_big, Interval(lower_bound, upper_bound))
    end
    return IntervalBox(x_big...)
end


function check_feasiblity(point::SVector{N,T}, constraints::Vector{Constraint{T}}) where{N, T}
    for constraint in constraints
        if !(constraint.C(point) ∈ constraint.bound)
            return false
        end
    end
    return true
end

function ibc_minimise(f::Function , X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}, C::BasicContractor; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, T}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), structure = SortedVector, debug = false, tol=1e-6, ibc_ind = true) where{N, T}


    working = structure([(X, inf(f(X)))], x->x[2]) # list of boxes with corresponding lower bound, arranged according to selected structure :
    minimizers = IntervalBox{N,T}[]
    global_min = ∞  # upper bound

    info = Information(0, 0, 0)
    num_bisections = 0

    while !isempty(working)

        info.iterations= info.iterations + 1

        if !ibc_ind
            if isready(ibc_chnl)
                from_diff = take!(ibc_chnl)     # Receiving best individual from ibc_minimise
                if debug
                    println("Box recevied from DifferentialEvolution: ", from_diff)
                end
                global_min = min(from_diff[2], global_min)
                info.de_to_ibc = info.de_to_ibc + 1
            end
        end

        (X, X_min) = popfirst!(working)

        gra = gradient(f, X.v)
        (lbb, fcb, cb) = bauman_form(X, f, gra)      # ----- second order form

        if global_min > lbb
            X_min = max(X_min, lbb)
            if fcb < global_min && check_feasiblity(cb, constraints)
                global_min = fcb
                x_best = cb
            end
        else
            continue
        end

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

        status, output = generate_random_feasible_point(X, constraints)   # finding a feasible point in the interval box if there present any

        if status                           # true if we got any feasible point
            feas_point = output
            m = sup(f(Interval.(feas_point)))  # find candidate for upper bound of global minimum by evaluating f on a feasible point
            if m < global_min
                global_min = m
                x_best = SVector(mid(X))
                if !ibc_ind
                    if debug
                        println("Box send to DifferentialEvolution: ", x_best )
                    end
                    if info.iterations % 200 == 0
                       x_big = convex_hull(working.data)
                        put!(diffevol_chnl, (x_best, global_min, x_big))  # sending best individual to diffevol
                    else
                        put!(diffevol_chnl, (x_best, global_min, nothing))
                    end
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
        else
            if diam(X) < tol
                push!(minimizers, X)
            else
                X1, X2 = bisect(X)
                push!( working, (X1, inf(f(X1))) )
                push!( working, (X2, inf(f(X2))) )
                num_bisections += 1
            end
        end

    end

    if debug
        println("IBC search is ended")
    end

    if !ibc_ind
        if isready(diffevol_chnl)
            take!(diffevol_chnl)
            put!(diffevol_chnl,(SVector(mid(X)), Inf, nothing) )
        else
            put!(diffevol_chnl,(SVector(mid(X)), Inf, nothing) )
        end
        if debug
            println("DifferentialEvolution is also terminated")
        end
        if isready(ibc_chnl)
            take!(ibc_chnl)
        end
    end

    lower_bound = minimum(inf.(f.(minimizers)))

    return Interval(lower_bound, global_min), minimizers, info

end


function ibc_maximise(f::Function , X::IntervalBox{N,T}, constraints::Vector{Constraint{T}}, C::BasicContractor; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, T}}(0)),
               diffevol_chnl = RemoteChannel(()->Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}(0)), structure = SortedVector, debug = false, tol=1e-6, ibc_ind = true) where{N, T}
    bound, minimizer, info = ibc_minimise(x -> -f(x), X, constraints, C, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, debug = debug, tol = tol, ibc_ind = ibc_ind)
    return -bound, minimizer, info
end
