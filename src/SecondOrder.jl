function bauman_form(X::T, f::function, g::Array{Interval{Float64},1}) where{T}
    cb = Float64[]

    for i in 1:length(G)
        if g[i] >= 0
            append!(cb, X[i].lo)
        elseif g[i] <= 0
            append!(cb, X[i].hi)
        else
            append!(cb,(g[i].hi * X[i].lo - g[i].lo * X[i].hi)/(g[i].hi - g[i].lo))
        end
    end

    (sum1, sum2)  = (0, 0)
    for i in 1:length(g)
        sum1 = sum1 + (g[i].hi)*(X[i].lo - cb[i])
        sum2 = sum2 + (g[i].lo)*(X[i].hi - cb[i])
    end

    lb = f(cb...) + min(sum1, sum2)

    return (lb, f(cb...), cb)



function second_order(X::T, f::function, lb, ub, g::Array{Interval{Float64},1}) where{T}

    (lbb, fcb, cb) = bauman_form(X, f, g)
    if ub > lbb:
        lb = max(lb, lbb)
        if fcb < ub:
            ub = fcb
        end
    end

    return (lb, ub)
