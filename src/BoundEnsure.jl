"""There is the possibility of the created vector existing outside of the bounds specified at initialization.
 Because of this, we created a separate function i.e BoundEnsure that checks and corrects for this. In the event
 that one of these rogue points are found, we’ll simply move it to the nearest boundary """
 
function BoundEnsure(V::Vector{T}, X::IntervalBox{N, T}) where{N, T}

    for i in length(X)
      if V[i] < X[i].lo
        V[i] = X[i].lo
      elseif X[i].hi < V[i]
        V[i] = X[i].hi
      end
    end

    return V
end
