function maxdist(x::SVector{N, T}, X::IntervalBox{N, T}) where{N, T}
    d = 0.0
    for i in 1:N
        if x[i] < X[i].lo
            d = d + (X[i].lo - x[i])^2
        elseif  X[i].hi < x[i]
            d = d + (x[i] - X[i].hi)^2
        end
    end
    return d
end


function convex_hull(vec::Vector{IntervalBox{N,T}}) where{N, T}

    x_big = Interval{T}[]
    num_variables = N
    num_boxes = length(vec)

    for i in 1:num_variables
        lower_bound = Inf
        upper_bound = -Inf
        for j in 1:num_boxes
            if vec[j][i].lo < lower_bound
                lower_bound = vec[j][i].lo
            end
            if upper_bound < vec[j][i].hi
                upper_bound = vec[j][i].hi
            end
        end
        append!(x_big, Interval(lower_bound, upper_bound))
    end
    return IntervalBox(x_big...)
end
