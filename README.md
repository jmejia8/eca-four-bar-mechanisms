#  Optimal Synthesis of four-bar mechanisms using Evolutionary Centers Algorithm

This work employs ECA for generating optimal four-bar mechanisms.

## Requirements

* [Julia 0.6](https://julialang.org/)
* [Metaheuristics package](https://github.com/jmejia8/Metaheuristics.jl)
* [Plots](https://github.com/JuliaPlots/Plots.jl)

## Example

```julia
include("findStructure.jl")

ncase = 3

println("Optimize...")
p, er, precisionpts = findStructure(ncase)


println("Generate animation...")
img = animateFourBarMechanism(p;
                precision_points=precisionpts,
                title   = @sprintf("Error: %e", er),
                xlimits = (-5, 5),
                ylimits = (-5, 5)
                )
```


## Study Case 1
![Mechanism](https://www.candaana.com/eca/case1.gif "Mechanism")
## Study Case 2
![Mechanism 2](https://www.candaana.com/eca/case2.gif "Mechanism 2")
## Study Case 3
![Mechanis 2m](https://www.candaana.com/eca/case3.gif "Mechanism 2")
