#  Optimal Synthesis of four-bar mechanisms using Evolutionary Centers Algorithm

This work employs ECA for generating optimal four-bar mechanisms.

## Requirements

* [Julia 0.6](https://julialang.org/)
* [Metaheuristics package](https://github.com/jmejia8/Metaheuristics.jl)
* [Plots](https://github.com/JuliaPlots/Plots.jl)

## Example

```julia
include("findStructure.jl")
include("animation.jl")

ncase   = 1
nframes = 50
# case information
D, bounds, precision_points, error_func = case_info(ncase)

# solution
p, my_error = findStructure(ncase)

C, X0, X1, X2, X3 = generateTrayectory(p, nframes)

# create an animation
getAnimation(C, X0, X1, X2, X3, nframes, my_error, precision_points)
```


## Study Case 1
![Mechanism](https://www.candaana.com/eca/case1.gif "Mechanism")
## Study Case 2
![Mechanism 2](https://www.candaana.com/eca/case2.gif "Mechanism 2")
## Study Case 3
![Mechanis 2m](https://www.candaana.com/eca/case3.gif "Mechanism 2")
