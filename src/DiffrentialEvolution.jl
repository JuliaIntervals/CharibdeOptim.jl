function DiffEvolution_min(f::Function, X::T, maxiter = 30 ) where {T}

   n = length(X)
   np = 10*n

   Pop = []                          #Initialsing Population
   for i in 1:np
      indv = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:n]
      push!(Pop, indv)
   end

   global_min = Inf
   X_best = last(Pop)

   for iter in 1:maxiter

      F = 2*rand()
      I = rand(1:n)
      CR = rand()
      PopNew = []

      for i in 1:np

         u = GenerateRandom(1, np, i)
         v = GenerateRandom(1, np, i, u)
         w = GenerateRandom(1, np, i, u, v)   # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         M = BoundEnsure(Pop[u] + F*(Pop[v] - Pop[w]), Pop[u], X)
                                             # Mutation : Mutant Vector is created
         for j in 1:n                        # Recombination or CrossOver :  Mutant vector is itself is modified by Crossover rate (CR)
            if j != I
               if rand() > CR
                  M[j] = Pop[i][j]
               end
            end
         end

         if f(M...) <= f(Pop[i]...)          # Selection : Best Indivual is selected among the Modified Mutant Individual and in the originl Individual
            push!(PopNew, M)
         else
            push!(PopNew, Pop[i])
         end

         if f(PopNew[i]...) < global_min
            global_min = f(PopNew[i]...)
            X_best = PopNew[i]
         end
      end
      Pop = PopNew
   end
   return global_min, X_best                # best individual is output
end

function DiffEvolution_max(f::Function, X::T, maxiter = 30) where {T}
    maxima, maximiser=  DiffEvolution_min(x -> -f(x), X, maxiter)
    return -maxima, maximiser
end
