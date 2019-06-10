function ibc_minimise(g::Function, X::T; ibc_chnl = RemoteChannel(()->Channel{Tuple{T, Float64}}(0)), diffevol_chnl = RemoteChannel(()->Channel{Tuple{Vector{Float64}, Float64}}(100)), structure = SortedVector, tol=1e-3 ) where {T}

    # list of boxes with corresponding lower bound, arranged according to selected structure :
    f = X -> g(X...)

    working = structure([(X, inf(f(X)))], x->x[2])
    minimizers = T[]
    global_min = ∞  # upper bound

    num_bisections = 0

    while !isempty(working)

        if isready(ibc_chnl)
            from_diff = take!(ibc_chnl)     # Receiving best individual from ibc_minimise
            global_min = min(from_diff[2], global_min)
        end

        (X, X_min) = popfirst!(working)

        if X_min > global_min    # X_min == inf(f(X))
            continue
        end

        # find candidate for upper bound of global minimum by just evaluating a point in the interval:
        m = sup(f(Interval.(mid.(X))))   # evaluate at midpoint of current interval

        if m < global_min
            global_min = m
            put!(diffevol_chnl, (Vector{Float64}(mid(X)), global_min)) # sending best individual to DiffEvaluation
        end

        # Remove all boxes whose lower bound is greater than the current one:
        filter_elements!(working , (X, global_min) )

        if diam(X) < tol
            push!(minimizers, X)
        else
            X1, X2 = bisect(X)
            push!( working, (X1, inf(f(X1))) )
            push!( working, (X2, inf(f(X2))) )
            num_bisections += 1
        end

    end

    lower_bound = minimum(inf.(f.(minimizers)))

    return Interval(lower_bound,global_min), minimizers
end
