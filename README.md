#  Optimal Synthesis of four-bar mechanisms using Evolutionary Centers Algorithm

This works employs ECA for generate optimal four-bar mechanisms.

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