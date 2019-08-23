function bauman_form(X::T, f::Function, g::Array{Interval{Float64},1}) where{T}
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

    lbb = f(cb...) + min(sum1, sum2)

    return (lbb, f(cb...), cb)

end 
