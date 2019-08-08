"""Usage:
```
For Unconstrained Optimsation:
  f = X->((x,y)=X;x^3 + 2y + 5)
  A = IntervalBox(2..4, 2..3)
  (global_min, minimisers, info) = ibc_minimise(f, A)
  (global_max, maximisers, info) = ibc_maximise(f, A)

For Constrained Optimisation:
  f = X->((x,y)=X;-(x-4)^2-(y-4)^2)
  A = IntervalBox(-4..4, -4..4)

  vars = ModelingToolkit.@variables x y
  C1 = Constraint(vars, x+y, -Inf..4)
  C2 = Constraint(vars, x+3y, -Inf..9)

  (global_min, minimisers, info) = ibc_minimise(f, A, [C1, C2])
  (global_max, maximisers, info) = ibc_maximise(f, A, [C1, C2])

ibc_minimise/ibc_maximise find the global minimum/maximum value of the function in given search space by using Interval Bound & Contract(IBC) algorithm
```
"""
function ibc_minimise(f::Function , X::IntervalBox{N,T}; debug = false,  ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)), diffevol_chnl = Nothing, structure = SortedVector, tol=1e-6 ) where{N, T}

    vars = [Variable(Symbol("x",i))() for i in 1:length(X)]
    g(x...) = f(x)
    C = BasicContractor(vars, g)


    working = [(X, inf(f(X)))] # list of boxes with corresponding lower bound, arranged according to selected structure :

    minimizers = IntervalBox{N,T}[]
    global_min = ∞  # upper bound
    x_best = mid(X)

    info = Information(0, 0, 0)
    num_bisections = 0

    while !isempty(working)

        info.iterations= info.iterations + 1

        if isready(ibc_chnl)
            from_diff = take!(ibc_chnl)     # Receiving best individual from ibc_minimise
            if debug
                println("Box recevied from DifferentialEvolution: ", from_diff)
            end
            if from_diff[2] < global_min
                x_best = from_diff[1]
                global_min = from_diff[2]
                working = sort(working, by = x-> maxdist(x_best, x[1]))
            end

            info.de_to_ibc = info.de_to_ibc + 1
        end
        (X, X_min) = pop!(working)

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
            x_best = SVector(mid(X))
            working = sort(working, by = x-> maxdist(x_best, x[1]))
            if diffevol_chnl != Nothing
                if debug
                    println("Box send to DifferentialEvolution: ", x_best )
                end
                if info.iterations % 200 == 0
                   x_big = convex_hull(working)
                    put!(diffevol_chnl, (x_best, global_min, x_big))  # sending best individual to diffevol
                else
                    put!(diffevol_chnl, (x_best, global_min, nothing))
                end
                info.ibc_to_de = info.ibc_to_de + 1
            end
        end


        working = filter(x-> x[2] < global_min, working)   # Remove all boxes whose lower bound is greater than the current one:

        if diam(X) < tol
            push!(minimizers, X)
        else
            X1, X2 = bisect(X)
            push!( working, (X1, inf(f(X1))) )
            push!( working, (X2, inf(f(X2))) )
            num_bisections += 1
        end
    end

    if debug
        println("IBC search is ended")
    end

    if diffevol_chnl != Nothing
        if isready(diffevol_chnl)
            take!(diffevol_chnl)
            put!(diffevol_chnl,(SVector(mid(X)), Inf, nothing) )
        else
            put!(diffevol_chnl,(SVector(mid(X)), Inf, nothing) )
        end
        if debug
            println("DifferentialEvolution is also terminated")
        end
    end

    if isready(ibc_chnl)
        take!(ibc_chnl)
    end

    lower_bound = minimum(inf.(f.(minimizers)))

    return Interval(lower_bound,global_min), minimizers, info, x_best


end

function ibc_maximise(f::Function , X::IntervalBox{N,T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)), diffevol_chnl = Nothing, structure = SortedVector, debug = false, tol=1e-6) where{N, T}
    bound, minimizer, info = ibc_minimise(x -> -f(x), X, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, debug = debug, tol = tol)
    return -bound, minimizer, info
end
