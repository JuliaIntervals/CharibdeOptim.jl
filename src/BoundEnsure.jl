"""There is the possibility of the created vector existing outside of the bounds specified at initialization.
 Because of this, we created a separate function i.e BoundEnsure that checks and corrects for this. In the event
 that one of these rogue points are found, we’ll simply move it to the nearest boundary """

function bound_ensure(m::SVector{N, T}, u::SVector{N, T}, X::IntervalBox{N, T}) where {N, T}

    for i in 1:length(X)
      if m[i] < X[i].lo
        m = setindex(m, u[i]+rand()*(X[i].lo-u[i]), i)
      elseif X[i].hi < m[i]
        m = setindex(m, u[i]+rand()*(X[i].hi-u[i]), i)
      end
    end
    return m
end
