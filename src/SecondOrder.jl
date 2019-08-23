function bauman_form(X::T, f::Function, g) where{T}
    cb = Float64[]

    for i in 1:length(g)
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

    func_val = f(cb)
    static_cb = SVector{length(cb)}(cb)

    lbb = func_val + min(sum1, sum2)

    return (lbb, func_val, static_cb)

end
