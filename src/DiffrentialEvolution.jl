function diffevol_minimise(f::Function, X::T, ibc_chnl::RemoteChannel{Channel{Tuple{T,Float64}}}, diffevol_chnl::RemoteChannel{Channel{Tuple{Array{Float64,1},Float64}}} ) where {T}

   n = length(X)
   np = 10*n

   pop = Vector{Float64}[]                          #Initialsing Population
   for i in 1:np
      indv = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:n]
      push!(pop, indv)
   end

   global_min = Inf
   x_best = last(pop)

   while true

      fac = 2*rand()
      ind = rand(1:n)
      cr = rand()
      pop_new = Vector{Float64}[]

      temp = global_min

      if isready(diffevol_chnl)
         (x_best, temp) = take!(diffevol_chnl)  # Receiveing best individual from diffevol_minimise
         if temp == Inf
            break
         end
         push!(pop, x_best)
         np = np + 1
      end
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

      if global_min < temp
         put!(ibc_chnl, (IntervalBox(Interval.(x_best)), global_min))   #sending the best individual to ibc_minimise
      end

      pop = pop_new
   end
   return
end
