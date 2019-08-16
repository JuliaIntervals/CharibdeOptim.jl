"""
diffevol_minimise is a function which the implementation of DifferentialEvolution algorithm is designed in such a way so that it can only be used by charibde_min on different processor.
It can be use for both constrained as well as for unconstrained optimisation
"""
function diffevol_minimise(f::Function, X::IntervalBox{N, T}, ibc_chnl::Union{Channel{Tuple{IntervalBox{N,T}, T}}, RemoteChannel{Channel{Tuple{IntervalBox{N, T}, T}}} },
               diffevol_chnl::Union{Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}, RemoteChannel{Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}}}; np = 10*N, debug = false ) where{N, T}

   pop = SVector{N, T}[]

   for i in 1:np
      indv = [X[j].lo + (1-rand())*(X[j].hi - X[j].lo) for j in 1:N]                        #Initialsing Population
      push!(pop, SVector{N, T}(indv))
   end

   global_min = Inf
   x_best = pop[np]

   while true

      fac = 2*rand()
      ind = rand(1:N)
      cr = rand()
      pop_new = SVector{N, T}[]

      temp = global_min

      if isready(diffevol_chnl)
         (x_best, temp, new_box) = take!(diffevol_chnl)  # Receiveing best individual from diffevol_minimise
         if debug
            println("Box recevied from IBC: ", (x_best, temp))
         end
         if temp == Inf
            break
         end
         if new_box != nothing
            X = new_box
         end
         push!(pop, x_best)
         np = np + 1
      end

      for i in 1:np

         u = generate_random(1, np, i)
         v = generate_random(1, np, i, u)
         w = generate_random(1, np, i, u, v)   # Choosing index of three different individuals, different from the index of that individual whose mutant vector is going to form.

         m = bound_ensure(pop[u] + fac * (pop[v] - pop[w]), pop[u], X) # Mutation : Mutant Vector is created

         for j in 1:N                        # Recombination or CrossOver :  Mutant vector is itself is modified by Crossover rate (CR)
            if j != ind
               if rand() > cr
                  m = setindex(m, pop[i][j], j)
               end
            end
         end

         if f(m) <= f(pop[i])          # Selection : Best Indivual is selected among the Modified Mutant Individual and in the originl Individual
            push!(pop_new, m)
         else
            push!(pop_new, pop[i])
         end

         if f(pop_new[i]) < global_min
            global_min = f(pop_new[i])
            x_best = pop_new[i]
         end
      end

      if global_min < temp
         best_box = (IntervalBox(Interval.(x_best)), global_min)
         put!(ibc_chnl, best_box)   #sending the best individual to ibc_minimise
         if debug
            println("Box send to IBC: ", best_box)
         end
      end

      pop = pop_new
   end
   return
end

function diffevol_maximise(f::Function, X::IntervalBox{N, T}, ibc_chnl::Union{Channel{Tuple{IntervalBox{N,T}, T}}, RemoteChannel{Channel{Tuple{IntervalBox{N,T}, T}}} },
               diffevol_chnl::Union{Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}, RemoteChannel{Channel{Tuple{SVector{N, T}, T, Union{Nothing, IntervalBox{N, T}}}}}}; np = 10*N, debug = false ) where{N, T}

            diffevol_minimise(x -> -f(x), X, ibc_chnl, diffevol_chnl, np = np, debug = debug)
         end
