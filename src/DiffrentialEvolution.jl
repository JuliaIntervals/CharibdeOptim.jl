function DiffEvolution(f::Funnction, X::IntervalBox{N, T} ) where {N, T}

   n = length(X)
   np = 10*n

   Pop = []                          #Initialsing Population
   ran = []
   FuncVal = []
   for i in 1:np
      indv = []
      ranv = []
      for j in 1:n
         r = 1 - rand()
         x = X[j].lo + r*(X[j].hi - X[j].lo)
         appned!(indv, x)
         append!(ranv, r)
      end
      append!(Pop, indv)
      append!(ran, ranv)
      append!(FuncVal, f(indv...))
   end

   while 1
      F = 2*rand()
      I = rand(1:n)
      CR = rand()
      PopNew = []
      for i in 1:np
         u = generateRandom(1, np, i)
         v = generateRandom(1, np, i, u)
         w = generateRandom(1, np, i, u, v)
         M = Pop[u] + F*(Pop[v] - Pop[w])
         O = []
         for j in 1:d
            if j != I
               if ran[i][j] <= CR
                  append!(O,M[j])
               else
                  append!(O, Pop[i][j])
               end
            else
               append!(O, M[j])
            end
         end
         if f(O...) < f(Pop[i]...)
            append!(PopNew, O)
         else
            append!(popNew,Pop[i])
         end
      end

      for i in 1:np
         Pop[i] = PopNew[i]
         FuncVal[i] = f(Pop[i]...)
       end
   end                                                #To do : Termination crieteria have to figure out

return min(FuncVal...)

end

function generateRandom(min, max, x...)
   r = rand(min:max)
   if r in Set(x)
      return generateRandom(min, max, x...)
   else
      return r
end
