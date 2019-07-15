function diffevol_minimise(f::Function, X::IntervalBox{N, T}, constraints::Vector{ConstraintCond{T}}, ibc_chnl::Union{Channel{Tuple{IntervalBox{N,T}, Float64}}, RemoteChannel{Channel{Tuple{IntervalBox{N,T},Float64}}} },
               diffevol_chnl::Union{Channel{Tuple{SArray{Tuple{N},Float64,1,N},Float64}}, RemoteChannel{Channel{Tuple{SArray{Tuple{N},Float64,1,N}, Float64}}}} ) where{N, T}

   n = length(X)
   np = 10*n

   pop = SArray{Tuple{n},Float64,1,n}[]

   for i in 1:np
      indv = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:n]                        #Initialsing Population
      push!(pop, SVector{n, Float64}(indv))
   end

   global_min = Inf
   x_best = pop[np]

   while true

      fac = 2*rand()
      ind = rand(1:n)
      cr = rand()
      pop_new = SArray{Tuple{n},Float64,1,n}[]

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

         m = bound_ensure(pop[u]+fac*(pop[v]-pop[w]), pop[u], X) # Mutation : Mutant Vector is created

         for j in 1:n                        # Recombination or CrossOver :  Mutant vector is itself is modified by Crossover rate (CR)
            if j != ind
               if rand() > cr
                  m = setindex(m, pop[i][j], j)
               end
            end
         end

         (c1, c2) = (0, 0)

         for constraint in constraints
            if constraint.C(m...) ∈ constraint.bound
               c1 = c1 + 1
            end
            if constraint.C(pop[i]...) ∈ constraint.bound
               c2 = c2 + 1
            end
         end

         if (c1, c2) == (length(constraints), length(constraints))
            if f(m...) < f(pop[i]...)
               push!(pop_new, m)
            else
               push!(pop_new, pop[i])
            end
         else
            if c2 <= c1
               push!(pop_new, m)
            else
               push!(pop_new, pop[i])
            end
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

end
