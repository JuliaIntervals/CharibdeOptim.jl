struct ConstraintCond{T}
   f ::Operation
   c ::Interval{T}
end


function diff_minimise(f::Function, X::T, constraints::Vector{ConstraintCond{T}}; maxiter = 30 ) where{T}

   n = length(X)
   np = 10*n
   nc = length(constraints)

   Pop = Array{Float64,1}[]                          #Initialsing Population
   for i in 1:np
      indv = [X[j].lo + r*(X[j].hi - X[j].lo) for j in 1:n]
      push!(Pop, indv)
   end

   for iter in maxiter

      F = 2*rand()
      I = rand(1:n)
      CR = rand()
      PopNew = Array{Float64,1}[]

      for i in 1:np

         u = generate_random(1, np, i)
         v = generate_random(1, np, i, u)
         w = generate_random(1, np, i, u, v)    # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         M = BoundEnsure(Pop[u] + F*(Pop[v] - Pop[w]), Pop[u], X)                # Mutatation : Mutant Vector is created

         for j in 1:n                          # Recombination or CrossOver :  Mutant vector is itself modified by Crossover rate (CR)
            if j != I
               if rand() > CR
                  M[j] = Pop[i][j]
               end
            end
         end

         C1, C2 =0, 0                                # Selection : Best individual (among modified Mutant and Original individual) is being selected by their response toward Constraints

         for constraint in constraints
            if constraint.f(M...) ∈ constraint.c
               C1 = C1 + 1
            end
            if constraint.f(Pop[i]...) ∈ constraint.c
               C2 = C2 + 1
            end
         end

         if (C1, C2) == (nc, nc)
            if f(M...) < f(Pop[i]...)
               push!(PopNew, M)
            else
               push!(popNew, Pop[i])
            end
         else
            if C2 <= C1
               push!(Popnew, M)
            else
               push!(popNew, Pop[i])
            end
         end

        if f(PopNew[i]...) < global_min
           global_min = f(PopNew[i]...)
           X_best = PopNew[i]
        end
      end
      Pop = PopNew
  end
  return global_min, X_best
end


function diff_maximise(f::Function, X::T, constraints::Vector{ConstraintCond{T}}, maxiter = 30 ) where{T}
   maxima, maximiser = DiffEvolution_min(f, X, constraints, maxiter)
   return -maxima, maximiser
end
