function ibc_minimise(g::Function, X::T; ibc_chnl = RemoteChannel(()->Channel{Tuple{T, Float64}}(0)), diffevol_chnl = Nothing, tol=1e-3 ) where {T}

    vars = [Variable(Symbol("x",i))() for i in 1:length(X)]
    C = Contractor(vars, g)
    f = X -> g(X...)

    working = [X] # list of boxes with corresponding lower bound, arranged according to selected structure :
    x_best = IntervalBox(Interval.(mid.(X)))
    minimizers = T[]
    global_min = ∞  # upper bound

    iter = 0

    while !isempty(working)

        if isready(ibc_chnl)
            from_diff = take!(ibc_chnl)     # Receiving best individual from ibc_minimise
            if from_diff[2] < global_min
                (x_best, global_min) = from_diff
                sort!(working, by = x -> maxdist(x, x_best))
            end
        end

        X = pop!(working)

        X = invokelatest(C, -∞..global_min, X)                        # Contracting the box by constraint f(X) < globla_min
        X_min = inf(f(X))

        if X_min > global_min    # X_min == inf(f(X))
            continue
        end

        # find candidate for upper bound of global minimum by just evaluating a point in the interval:
        m = sup(f(Interval.(mid.(X))))   # evaluate at midpoint of current interval

        if m < global_min
            global_min = m
            x_best = IntervalBox(Interval.(mid.(X)))
            if diffevol_chnl != Nothing put!(diffevol_chnl, (Vector{Float64}(mid(X)), global_min)) end   # sending best individual to DiffEvaluation
        end

        working = filter(x -> inf(f(x))<global_min, working)

        if diam(X) < tol
            push!(minimizers, X)
        else
            X1, X2 = bisect(X)
            push!(working, X1, X2)
            sort!(working, by = x -> maxdist(x, x_best))
        end
        iter = iter + 1
    end

    println(iter)
    if diffevol_chnl != Nothing
        if isready(diffevol_chnl)
            take!(diffevol_chnl)
            put!(diffevol_chnl,(Vector{Float64}(mid(X)), Inf) )
        else
            put!(diffevol_chnl,(Vector{Float64}(mid(X)), Inf) )
        end
    end

    if isready(ibc_chnl)
        take!(ibc_chnl)
    end

    lower_bound = minimum(inf.(f.(minimizers)))

    return Interval(lower_bound,global_min), minimizers
    #return (global_min, x_best)
end
