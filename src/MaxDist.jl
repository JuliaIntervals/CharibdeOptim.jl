function maxdist(X::T, x::T) where{T}
    d = 0.0
    for i in 1:length(x)
        if x[i].lo < X[i].lo
            d = d + (X[i].lo - x[i].lo)^2
        elseif  X[i].hi < x[i].lo
            d = d + (x[i].lo - X[i].hi)^2
        end
    end
    return d
end
