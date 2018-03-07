using Metaheuristics
include("tools.jl")


function myerror(x, error_func, precision_points)
    f, g = error_func(x, grashof(x), precision_points)

    return f + 1000g
end

function findStructure(ncase)
    # case information
    D, bounds, precision_points, error_func = case_info(ncase)

    # fitness
    fitness(x) = myerror(x, error_func, precision_points)
    
    # optimize
    return eca(fitness, D; saveLast = "last.csv",
                           saveConvergence = "conv.csv",
                           max_evals=20000D,
                           limits=bounds)


end