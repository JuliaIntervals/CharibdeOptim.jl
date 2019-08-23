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


    working = structure([(X, inf(f(X)))], x->x[2]) # list of boxes with corresponding lower bound, arranged according to selected structure :

    minimizers = IntervalBox{N,T}[]
    global_min = ∞  # upper bound

    info = Information(0, 0, 0)
    num_bisections = 0

    while !isempty(working)

        info.iterations= info.iterations + 1

        if isready(ibc_chnl)
            from_diff = take!(ibc_chnl)     # Receiving best individual from ibc_minimise
            if debug
                println("Box recevied from DifferentialEvolution: ", from_diff)
            end
            global_min = min(from_diff[2], global_min)
            info.de_to_ibc = info.de_to_ibc + 1
        end

        (X, X_min) = popfirst!(working)
        g = ForwardDiff.gradient(f, X)
        (lbb, fcb, cb) = bauman_form(X, f, g)      # ----- second order form

        if global_min > lbb
            lb = max(lb, lbb)
            if fcb < global_min:
                global_min = fcb
                x_best = SVector(cb)
            end
        else
            continue
        end                                        # ------

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
            if diffevol_chnl != Nothing
                if debug
                    println("Box send to DifferentialEvolution: ", x_best )
                end
                if info.iterations % 200 == 0
                    put!(diffevol_chnl, (x_best, global_min, X))  # sending best individual to diffevol
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

    return Interval(lower_bound,global_min), minimizers, info


end

function ibc_maximise(f::Function , X::IntervalBox{N,T}; ibc_chnl = RemoteChannel(()->Channel{Tuple{IntervalBox{N,T}, Float64}}(0)), diffevol_chnl = Nothing, structure = SortedVector, debug = false, tol=1e-6) where{N, T}
    bound, minimizer, info = ibc_minimise(x -> -f(x), X, ibc_chnl = ibc_chnl, diffevol_chnl = diffevol_chnl, debug = debug, tol = tol)
    return -bound, minimizer, info
end
