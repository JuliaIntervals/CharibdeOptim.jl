struct ConstraintCond{T}
   f ::Operation
   c ::Interval{T}
end


function DiffEvolution(f::Funnction, X::IntervalBox{N, T}, constraints::Vector{ConstraintCond{T}}) where {N, T}

   n = length(X)
   np = 10*n
   nc = length(constraints)

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
         u = generateRandom(1, np, i)         #To do: have to modify the strategy to choose u(base individual)
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

         C1, C2 =0, 0                          #Constraint Handeling
         for constraint in constraints
            if constraint.f(O...) ∈ constraint.c:
               C1 = C1 + 1
            end
            if constraint.f(Pop[i]...) ∈ constraint.c:
               C2 = C2 + 1
            end
         end

         if C1 == nc & C2 == nc :
            if f(O...) < f(Pop[i]...)
               append!(PopNew, O)
            else
               append!(popNew,Pop[i])
            end
         else
            if C2 <= C1:
               append!(Popnew, O)
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
