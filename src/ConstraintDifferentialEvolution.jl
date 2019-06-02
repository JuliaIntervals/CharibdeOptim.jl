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

      fac = 2*rand()
      ind = rand(1:n)
      cr = rand()
      PopNew = Array{Float64,1}[]

      for i in 1:np

         u = generate_random(1, np, i)
         v = generate_random(1, np, i, u)
         w = generate_random(1, np, i, u, v)    # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         M = BoundEnsure(Pop[u] + fac*(Pop[v] - Pop[w]), Pop[u], X)                # Mutatation : Mutant Vector is created

         for j in 1:n                          # Recombination or CrossOver :  Mutant vector is itself modified by Crossover rate (CR)
            if j != ind
               if rand() > cr
                  M[j] = Pop[i][j]
               end
            end
         end

         c1, c2 =0, 0                                # Selection : Best individual (among modified Mutant and Original individual) is being selected by their response toward Constraints

         for constraint in constraints
            if constraint.f(M...) ∈ constraint.c
               c1 = c1 + 1
            end
            if constraint.f(Pop[i]...) ∈ constraint.c
               c2 = c2 + 1
            end
         end

         if (c1, c2) == (nc, nc)
            if f(M...) < f(Pop[i]...)
               push!(PopNew, M)
            else
               push!(popNew, Pop[i])
            end
         else
            if c2 <= c1
               push!(Popnew, M)
            else
               push!(popNew, Pop[i])
            end
         end

        if f(PopNew[i]...) < global_min
           global_min = f(PopNew[i]...)
           x_best = PopNew[i]
        end
      end
      Pop = PopNew
  end
  return global_min, x_best
end


function diff_maximise(f::Function, X::T, constraints::Vector{ConstraintCond{T}}, maxiter = 30 ) where{T}
   maxima, maximiser = DiffEvolution_min(f, X, constraints, maxiter)
   return -maxima, maximiser
end
