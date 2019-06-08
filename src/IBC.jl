function ibc_minimise(f::Function, X::T, ibc_chnl, diffevol_chnl; structure = SortedVector, tol=1e-3 ) where {T}

    # list of boxes with corresponding lower bound, arranged according to selected structure :
    working = structure([(X, ∞)], x->x[2])
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
        m = f(mid(X)...)   # evaluate at midpoint of current interval

        if m < global_min
            global_min = m
            put!(diffevol_chnl, (mid(X), global_min)) # sending best individual to DiffEvaluation
        end

        # Remove all boxes whose lower bound is greater than the current one:
        filter_elements!(working , (X, global_min) )

        if diam(X) < tol
            push!(minimizers, X)
        else
            X1, X2 = bisect(X)
            push!( working, (X1, inf(f(X1...))) )
            push!( working, (X2, inf(f(X2...))) )
            num_bisections += 1
        end

    end

    lower_bound = inf(f(minimizers[1]...))
    for i in 2:length(minimizers)
        lower_bound = min(inf(f(minimizers[i]...)), lower_bound)
    end

    return Interval(lower_bound, global_min), minimizers
end


function ibc_maximise(f, X::T; structure = HeapedVector, tol=1e-3 ) where {T}
    bound, minimizers = minimise(x -> -f(x), X, structure, tol)
    return -bound, minimizers
end
