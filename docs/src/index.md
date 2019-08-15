# CharibdeOptim.jl

This package gives you the global optimum value (max or min) of any mathematical objective function in a given domain by running two different algorithms (Interval Branch & Contract and Differential evolution) in parallel, this approach of solving mathematical optimum value problems of constrained and unconstrained type is called *Charibde* which can be used through functions `charibde_min` or `charibde_max` (to get minimum or maximum value of the objective function) .

The package also provides the function `ibc_minimise` and `ibc_maximise` to use *Interval Branch & Contract* algorithm to solve the problems.

To run Charibde on two workers(processors) there are few requirements

1. The Julia session should have at least 2 more workers other than the master processor with package `CharibdeOptim` loaded on them, which can we have through like this
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
   - `ibc_minimise` and `ibc_maximise` also supports it.
   - Can be set to `IntervalOptimisation.SortedVector` or `IntervalOptimisation.Heapedvector`.
   - Defines the exploration strategy of search-space in IBC, either the next unexplored box should be popped out from a sorted vector (boxes are ranked upon the value of lower bound of the objective function on the box) or from the minimum heap data structure.
   - Default is set to `IntervalOptimisation.SortedVector`.
2. `workers`
   - Tells on how many workers Charibde will work.
   - Should be set to 1 or 2.
   - Default is set to 2.
3. `tol`
   - `ibc_minimise` and `ibc_maximise` also supports it.
   - Defines the accuracy level of the solution.
   - Default is set to *1e-6*.
4. `np`
   - Defines the number of individuals in the initial population of Differential Evolution.
   - Advised to keep it around 10 to 20 times the dimension of search-space.
   - Default is set to 10 times the dimension of search-space.
5. `debug`
   - You have to set it to `true` to access debug mode of `ibc_minimise`/`ibc_maximise` or `charibde_min` and `charibde_max`.

## Acknowledge constraints

In case of Constrained Optimisation, you also need to define constraints and input them in *Charibde* or *IBC*.
To define each constraint the package provides the function `constraint` which takes three arguments
  - an array of variables (whose length should be equal to the dimension of the search space) defined using macro `@ModelingToolkit.variables`
  - expression of constraint in the form of `Operation` using `ModelingToolkit`
  - Interval box defining bound of that constraint
  and a keyword argument `epsilon` which define tolerance in case of equality constraint (default is set to *1e-4*)
and returns an object of type `Constraint`.

All these `Constraint` objects for every constraint should be defined on each worker using macro `@everywhere` while using *Charibde*.
Note: No need to use `@everywhere` anywhere while using *IBC*.

After defining every constraint using function `constraint` you need to pass all returned objects (of type `Constraint`) in an array along with objective function and search-space in `charibde_max`/`charibde_min` or `ibc_maximise`/`ibc_minimise`  according to the required optimum value (minimum or maximum).


## Usage

Functions `charibde_min`/`charibde_max` and `ibc_maximise`/`ibc_minimise` are provided to find the global minimum or maximum, respectively, of a standard Julia function f of one or several variables.

They return an Interval that is guaranteed to contain the global minimum (maximum), a Vector of Intervals or IntervalBoxes whose union contains all the minimisers and an object of type `Information` which have three fields (`de_to_ibc`, `ibc_to_de` and `iterations`) giving information about how many updates were sent from DE to IBC , from IBC to DE and how many total iterations occurred in IBC .

### Examples

#### Unconstrained Optimisation

##### Using IBC

```julia
julia> using CharibdeOptim

julia> using IntervalArithmetic

julia> (global_min, minimisers, info)= ibc_minimise(X->((x,y)= X;x^2 + y^2), IntervalBox(2..3, 3..4))
([13, 13.0001], IntervalBox{2,Float64}[[2, 2.00001] × [3, 3.00001], [2, 2.00001] × [3, 3.00001]], CharibdeOptim.Information(0, 0, 42))

julia> global_min
[13, 13.0001]

julia> minimisers
2-element Array{IntervalBox{2,Float64},1}:
 [2, 2.00001] × [3, 3.00001]
 [2, 2.00001] × [3, 3.00001]

julia> info
CharibdeOptim.Information(0, 0, 42)
```
##### Using Charibde

```julia
julia> using Distributed

julia> addprocs(2)
2-element Array{Int64,1}:
 2
 3

julia> @everywhere using CharibdeOptim

julia> @everywhere using IntervalArithmetic

julia> charibde_min(X->((x,y)=X;x^3 + 2y + 5), IntervalBox(2..4, 2..3))
([17, 17.0001], IntervalBox{2,Float64}[[2, 2.00001] × [2, 2.00001], [2, 2.00001] × [2, 2.00001]], CharibdeOptim.Information(23, 22, 23))

julia> (global_min, minimisers, info) = charibde_min(X->((x,y)=X;x^3 + 2y + 5), IntervalBox(2..4, 2..3))
([17, 17.0001], IntervalBox{2,Float64}[[2, 2.00001] × [2, 2.00001], [2, 2.00001] × [2, 2.00001], [2, 2.00001] × [2, 2.00001]], CharibdeOptim.Information(26, 26, 26))

julia> global_min
[17, 17.0001]

julia> minimisers
3-element Array{IntervalBox{2,Float64},1}:
 [2, 2.00001] × [2, 2.00001]
 [2, 2.00001] × [2, 2.00001]
 [2, 2.00001] × [2, 2.00001]

julia> info
CharibdeOptim.Information(26, 26, 26)
```

#### Constrained Optimisation

##### Using IBC

```julia
julia> using ModelingToolkit

julia> vars = ModelingToolkit.@variables x y
(x(), y())

julia> C1 = constraint(vars, x+y, -Inf..4)
CharibdeOptim.Constraint{Float64}([-∞, 4], Contractor in 2 dimensions:
  - forward pass contracts to 1 dimensions
  - variables: Symbol[:x, :y]
  - expression: x() + y())

julia> C2 = constraint(vars, x+3y, -Inf..9)
CharibdeOptim.Constraint{Float64}([-∞, 9], Contractor in 2 dimensions:
  - forward pass contracts to 1 dimensions
  - variables: Symbol[:x, :y]
  - expression: x() + 3 * y())

julia> (maxima, maximisers, info) = ibc_maximise(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4),[C1, C2])
([-8.00001, -7.99999], IntervalBox{2,Float64}[[1.99762, 1.99763] × [2.00237, 2.00238], [1.99768, 1.99769] × [2.00231, 2.00232], [2.00227, 2.00228] × [1.99772, 1.99773], [2.00221, 2.00222] × [1.99778, 1.99779], [1.99791, 1.99792] × [2.00208, 2.00209], [1.99797, 1.99798] × [2.00202, 2.00203], [1.9972, 1.99721] × [2.00279, 2.0028], [1.99739, 1.9974] × [2.0026, 2.00261], [1.99727, 1.99728] × [2.00272, 2.00273], [2.00212, 2.00213] × [1.99787, 1.99788]  …  [1.9993, 1.99931] × [2.00069, 2.0007], [1.99931, 1.99932] × [2.00068, 2.00069], [2.00068, 2.00069] × [1.99931, 1.99932], [2.00069, 2.00071] × [1.99929, 1.99931], [1.99928, 1.99929] × [2.00071, 2.00072], [2.00071, 2.00072] × [1.99928, 1.99929], [2.00075, 2.00076] × [1.99924, 1.99925], [2.00072, 2.00073] × [1.99927, 1.99928], [1.99932, 1.99933] × [2.00067, 2.00068], [2.0007, 2.00071] × [1.99929, 1.9993]], CharibdeOptim.Information(0, 0, 5827))
```

##### Using Charibde

```julia
julia> @everywhere using ModelingToolkit

julia> @everywhere vars = ModelingToolkit.@variables x y

julia> @everywhere C1 = constraint(vars, x+y, -Inf..4)

julia> @everywhere C2 = constraint(vars, x+3y, -Inf..9)

julia> (maxima, maximisers, info) = charibde_max(X->((x,y)=X;-(x-4)^2-(y-4)^2), IntervalBox(-4..4, -4..4), [C1, C2])
([-8, -7.99999], IntervalBox{2,Float64}[[2.00243, 2.00245] × [1.99755, 1.99757], [2.00237, 2.00238] × [1.99762, 1.99763], [1.99692, 1.99693] × [2.00307, 2.00308], [1.99699, 1.997] × [2.003, 2.00301], [2.00231, 2.00232] × [1.99768, 1.99769], [2.00224, 2.00225] × [1.99775, 1.99776], [1.99706, 1.99707] × [2.00293, 2.00294], [2.00291, 2.00292] × [1.99708, 1.99709], [2.00209, 2.0021] × [1.9979, 1.99791], [2.00216, 2.00217] × [1.99783, 1.99784]  …  [2.0008, 2.00081] × [1.99919, 1.9992], [1.9992, 1.99921] × [2.00079, 2.0008], [2.00077, 2.00078] × [1.99922, 1.99923], [1.99923, 1.99924] × [2.00076, 2.00077], [1.99922, 1.99923] × [2.00077, 2.00078], [2.00071, 2.00072] × [1.99928, 1.99929], [2.00073, 2.00074] × [1.99926, 1.99927], [1.99924, 1.99925] × [2.00075, 2.00076], [1.99925, 1.99926] × [2.00074, 2.00075], [2.00076, 2.00077] × [1.99923, 1.99924]], CharibdeOptim.Information(50, 3, 5777))

julia> maxima
[-8, -7.99999]
```
