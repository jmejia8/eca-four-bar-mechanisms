include("findStructure.jl")
include("animation.jl")

function mytest()
    ncase   = 1
    nframes = 50
    # case information
    D, bounds, precision_points, error_func = case_info(ncase)

    # solution
    p, my_error = findStructure(ncase)

    C, X0, X1, X2, X3 = generateTrayectory(p, nframes)

    getAnimation(C, X0, X1, X2, X3, nframes, my_error, precision_points)
end