using Distributed
addprocs(2)

@everyhwere using CharibdeOptim
@everyhwere using ModelingToolkit
@everywhere using IntervalArithmetic

# --- Michalewicz Problem ---

@everywhere function michalewicz(X)
      sum = 0
      for i in 1:length(X)
            sum = sum + sin(X[i])*(sin(i*(X[i]^2)/π))^20
      end
      return -sum
end

(global_min, minimisers, information) = charibde_min(michalewicz, IntervalBox(0..π, 50))


# --- Sine Envolope Problem ---

@everyhwere function sine_envolope(X)
      sum = 0
      for i in 1:length(X)-1
            sum = sum + ( 0.5 + ( sin(√(X[i+1]^2 + X[i]^2) -0.5)/(0.001*(X[i+1]^2 + X[i]^2) + 1) )^2)
      end
      return -sum
end

(global_min, minimisers, information) = charibde_min(sine_envolope, IntervalBox(-100..100, 5))


# --- Shekel Problem ---

@everywhere function shekel(X)
      sum = 0
      for i in 1:30
            s = 0
            for j in 1:length((X))
                  s = s + (X[j] - i*j)^2
            end
            sum = sum + (1 / (i + s))
      end
      return -sum
end

(global_min, minimisers, information) = charibde_min(shekel, IntervalBox(0..10, 5))


# --- Egg-Holder Problem ---

@everywhere function egg_holder(X)
      sum = 0
      for i in 1:length(X)-1
            sum = sum + ((X[i+1] + 47)sin(√abs(X[i+1] + 47 + X[i]/2)) + X[i]sin(√abs(X[i] - X[i+1] -47)))
      end
      return -sum
end

(global_min, minimisers, information) = charibde_min(egg_holder, IntervalBox(-512..512, 5))


# --- Rana Problem ---

@everywhere function rana(X)
      sum = 0
      for i in 1: length(X)-1
            sum = sum + (X[i]*cos(√(abs(X[i+1] + X[i] + 1))) * sin(√(abs(X[i+1] - X[i] + 1))) + (1+ X[i+1]) * sin(√(abs(X[i+1] + X[i] + 1))) * cos(√(abs(X[i+1] - X[i] + 1))))
      end
      return -sum
end

(global_min, minimisers, information) = charibde_min(rana , IntervalBox(-512..512, 5))


# --- Keane Problem ---

@everywhere function keane(X)
      temp_sum = 0
      temp_pro = 1
      temp_sum1 = 0
      for i in 1:length(X)
            temp_sum = temp_sum + cos(X[i])^4
            temp_pro = temp_pro * cos(X[i])^2
            temp_sum1 = temp_sum1 + i * X[i]^2
      end
      sum = abs(temp_sum - 2 * temp_pro) / √temp_sum1
      return -sum
end

@everyhwere vars = ModelingToolkit.@variables x y w z
@everywhere C1 = constraint(vars , x*y*w*z, Interval(0.75, Inf))
@everyhwere C2 = constriant(vars , x+y+z+w, Interval(-Inf, 30))

(global_min, minimisers, information) = charibde_min(keane, IntervalBox(0..10, 4), [C1, C2])
