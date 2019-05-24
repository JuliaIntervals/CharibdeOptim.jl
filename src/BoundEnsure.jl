"""There is the possibility of the created vector existing outside of the bounds specified at initialization.
 Because of this, we created a separate function i.e BoundEnsure that checks and corrects for this. In the event
 that one of these rogue points are found, we’ll simply move it to the nearest boundary """

function BoundEnsure(M, U, X)

    for i in 1:length(X)
      if M[i] < X[i].lo
        M[i] = U[i] + rand()*(X[i].lo - U[i])
      elseif X[i].hi < M[i]
        M[i] = U[i] + rand()*(X[i].hi - U[i])
      end
    end
    return M
end
