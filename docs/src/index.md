# CharibdeOptim.jl

This package gives you the global optimum value (max or min) of any mathematical objective function in given domain by running two different algorithms (Interval Branch & Contract and Differential evolution) in parallel, this approach of solving mathematical optimum value problems of constrained and unconstrained type is called *Charibde* which can be used through functions `charibde_min` or `charibde_max` (to get minimum or maximum value of the objective function) .

To run Charibde on two workers(processors) there are few requirements

1. The Julia session should have atleast 2 more workers other than the master processor with package `CharibdeOptim` loaded on them, which can we have through like this
  ```
  julia> using Distributed
  julia> addprocs(2)
  julia> @everywhere using CharibdeOptim
  ```  
2. Objective function and Search - Space (in form of `IntervalBox`) should be defined on both of the workers by using macro `@everywhere`, like this
  ```
  julia> @everywhere using IntervalArithmetic
  julia> @everywhere X = IntervalBox(4..5, 4..5)
  julia> @everyhwere f = X->((x,y) = X; x^2 + y^2)
  ```
And finally you can call `charibde_min` or `charibde_max` according to the required optimum value
 `julia> charibde_min(f, X)`

While using the package through JuMP, you really not have to care about all these, you just have to tell the number of workers (through `workers` keyword argument which is already set to 2) on them you want to run *Charibde* and the package automatically add required number of workers and load and define all the required packages and stuff on them.  

`charibde_min` supports following keyword arguments,

1. `structure`
   - Can be set to `IntervalOptimisation.SortedVector` or `IntervalOptimisation.Heapedvector`.
   - Defines the exploration strategy of search-space in IBC.
   - Default is set to `IntervalOptimisation.SortedVector`.
2. `workers`
   - Tells on how many workers Charibde will work.
   - Should be set to 1 or 2.
   - Default is set to 2.
3. `tol`
   - Defines the accuracy level of the solution.
   - Default is set to *1e-6*.
4. `np`
   - Defines number of individuals in initial population of Differential Evolution.
   - Advised to keep it around 10 to 20 times the dimension of search-space.
   - Default is set to 10 times the dimension of search-space.
5. `debug`
   - You have to set it to `true` to access debug mode of .

## Acknowledge constraints

In case of Constrained Optimisation you also need to define constraints and input them in *Charibde*.
To define each constraint the package provides the function `constraint` which takes three arguments
  - array of variables (whose length should be equal to dimension of the search space) defined using macro `@ModelingToolkit.variables`
  - expression constraint in form of `Operation` using `ModelingToolkit`
  - Interval box defining bound of that constraint
  and a keyword argument `epsilon` which define tolerance in case of equality constraint (default is set to *1e-4*)
and returns an object of type `Constraint`.

Note: All these `Constraint` objects for every constraint should be define on each worker using macro `@everywhere`.

After defining every constraint using function `constraint` you need to pass all returned objects (of type `Constraint`) in an array along with objective function and search-space in `charibde_max` or `charibde_min` according to the required optimum value (minimum or maximum).
