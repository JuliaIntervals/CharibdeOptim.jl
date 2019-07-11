mutable struct Information
    de_to_ibc::Int
    ibc_to_de::Int
    iterations::Int
end


function ibc_minimise(g::Function , X::IntervalBox{N,T}; debug = false,  ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)), diffevol_chnl = Nothing, structure = SortedVector, tol=1e-3 ) where{N, T}
  
    vars = [Variable(Symbol("x",i))() for i in 1:length(X)]
    C = Contractor(vars, g)
    f = x->g(x...)

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

        A = -∞..global_min
        X = invokelatest(C, A, X)                        # Contracting the box by constraint f(X) < globla_min
        X_min = inf(f(X))

        if debug
            println("Contracted search_space: ", X)
        end

        if X_min > global_min    # X_min == inf(f(X))
            continue
        end

        # find candidate for upper bound of global minimum by just evaluating a point in the interval:
        m = sup(f(Interval.(mid.(X))))   # evaluate at midpoint of current interval

        if m < global_min
            global_min = m
            if diffevol_chnl != Nothing
                put!(diffevol_chnl, (SVector(mid(X)), global_min))  # sending best individual to diffevol
                info.ibc_to_de = info.ibc_to_de + 1
                if debug
                    println("Box send to DifferentialEvolution: ", x_best)
                end
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
            put!(diffevol_chnl,(MVector(mid(X)), Inf) )
        else
            put!(diffevol_chnl,(MVector(mid(X)), Inf) )
        end
        if debug
            println("DifferentialEvolution is also terminated")
        end
    end

    if isready(ibc_chnl)
        take!(ibc_chnl)
    end

    lower_bound = minimum(inf.(f.(minimizers)))

    if (info.ibc_to_de, info.de_to_ibc) == (0, 0)
        return Interval(lower_bound,global_min), minimizers
    else
        return Interval(lower_bound,global_min), minimizers, info
    end

end
