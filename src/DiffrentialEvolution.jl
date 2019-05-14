function DiffEvolution(f::Funnction, X::IntervalBox{N, T}; maxiter = 30 ) where {N, T}

   n = length(X)
   np = 10*n

   Pop = []                          #Initialsing Population
   FuncVal = []
   for i in 1:np
      indv = []
      for j in 1:n
         r = 1 - rand()
         x = X[j].lo + r*(X[j].hi - X[j].lo)
         appned!(indv, x)
      end
      append!(Pop, indv)
      append!(FuncVal, f(indv...))
   end

   for iter in 1:maxiter

      F = 2*rand()
      I = rand(1:n)
      CR = rand()
      PopNew = []

      for i in 1:np

         u = generateRandom(1, np, i)
         v = generateRandom(1, np, i, u)
         w = generateRandom(1, np, i, u, v)   # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         M = Pop[u] + F*(Pop[v] - Pop[w])
         M = BoundEnsure(M, X)               # Mutation : Mutant Vector is created

         for j in 1:d                        # Recombination or CrossOver :  Mutant vector is itself is modified by Crossover rate (CR)
            if j != I
               if rand() > CR
                  M[j] = Pop[i][j]
               end
            end
         end

         if f(M...) < f(Pop[i]...)          # Selection : Best Indivual is selected among the Modified Mutant Individual and in the originl Individual
            append!(PopNew, M)
         else
            append!(popNew,Pop[i])
         end
      end

      for i in 1:np
         Pop[i] = PopNew[i]
         FuncVal[i] = f(Pop[i]...)
      end

   end

   return min(FuncVal...)                # best individual is output

end
