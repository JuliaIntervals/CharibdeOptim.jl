struct ConstraintCond{T}
   f ::Operation
   c ::Interval{T}
end


function DiffEvolution(f::Funnction, X::IntervalBox{N, T}, constraints::Vector{ConstraintCond{T}}; maxiter = 30 ) where {N, T}

   n = length(X)
   np = 10*n
   nc = length(constraints)

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

   for iter in maxiter

      F = 2*rand()
      I = rand(1:n)
      CR = rand()
      PopNew = []

      for i in 1:np

         u = GenerateRandom(1, np, i)
         v = GenerateRandom(1, np, i, u)
         w = GenerateRandom(1, np, i, u, v)    # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         M = Pop[u] + F*(Pop[v] - Pop[w])
         M = BoundEnsure(M, X)                 # Mutatation : Mutant Vector is created

         for j in 1:d                          # Recombination or CrossOver :  Mutant vector is itself modified by Crossover rate (CR)
            if j != I
               if rand() > CR
                  M[j] = Pop[i][j]
               end
            end
         end

         C1, C2 =0, 0                                # Selection : Best individual (among modified Mutant and Original individual) is being selected by their response toward Constraints

         for constraint in constraints
            if constraint.f(M...) ∈ constraint.c:
               C1 = C1 + 1
            end
            if constraint.f(Pop[i]...) ∈ constraint.c:
               C2 = C2 + 1
            end
         end

         if C1, C2 == nc, nc :
            if f(M...) < f(Pop[i]...)
               append!(PopNew, M)
            else
               append!(popNew, Pop[i])
            end
         else
            if C2 <= C1:
               append!(Popnew, M)
            else
               append!(popNew, Pop[i])
            end
         end

      for i in 1:np
         Pop[i] = PopNew[i]
         FuncVal[i] = f(Pop[i]...)
      end

   end                                                    #To do : Termination crieteria have to figure out

   return min(FuncVal...)

end
