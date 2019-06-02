function diffevol_minimise(f::Function, X::T, maxiter = 30 ) where {T}

   n = length(X)
   np = 10*n

   pop = Array{Float64,1}[]                          #Initialsing Population
   for i in 1:np
      indv = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:n]
      push!(pop, indv)
   end

   global_min = Inf
   x_best = last(pop)

   for iter in 1:maxiter

      fac = 2*rand()
      ind = rand(1:n)
      cr = rand()
      pop_new = Array{Float64,1}[]

      for i in 1:np

         u = generate_random(1, np, i)
         v = generate_random(1, np, i, u)
         w = generate_random(1, np, i, u, v)   # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         m = bound_ensure(pop[u] + fac*(pop[v] - pop[w]), pop[u], X)
                                             # Mutation : Mutant Vector is created
         for j in 1:n                        # Recombination or CrossOver :  Mutant vector is itself is modified by Crossover rate (CR)
            if j != ind
               if rand() > cr
                  m[j] = pop[i][j]
               end
            end
         end

         if f(m...) <= f(pop[i]...)          # Selection : Best Indivual is selected among the Modified Mutant Individual and in the originl Individual
            push!(pop_new, m)
         else
            push!(pop_new, pop[i])
         end

         if f(pop_new[i]...) < global_min
            global_min = f(pop_new[i]...)
            x_best = pop_new[i]
         end
      end
      pop = pop_new
   end
   return global_min, x_best                # best individual is output
end

function diffevol_maximise(f::Function, X::T, maxiter = 30) where {T}
    maxima, maximiser=  DiffEvolution_min(x -> -f(x), X, maxiter)
    return -maxima, maximiser
end
