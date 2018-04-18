using Metaheuristics
include("meritobar.jl")

function objFunc(x::Vector{Float64})
    F, g = merito1(x)
    if g < 0.0
        g = 0.0
    end

  return F + 100g
end 

function getKC(bounds::Vector{Float64} = [0.1, 50])
    # return optimal K
    D = 3

    # optimize using ECA
    K, F = eca(objFunc, D;
                   limits = bounds,
                 saveLast = "last_control.csv",
                 saveConvergence = "conv_control.csv",
                 # showIter = false,
                 # p_bin = 0.03,
                 η_max=4,
                 K = 3,
                 N = 10,
                 # N = 40,
                max_evals = 10000)

    # K, F = diffEvolution(objFunc, D;
    #                limits = bounds,
    #                N = 100,
    #              # saveLast = "gen.csv",
    #              # showIter = false,
    #              # p_bin = 0.03,
    #              CR_min = 0.8,
    #             CR_max = 1,
    #             F_min  = 0.3,
    #             F_max  = 0.9,
    #             individual = Metaheuristics.xfg_indiv,
    #             # showResults=false,
    #             max_evals = 100*200)

    return K, F
end

function getθ(KC)
    t0 = 0.0
    tf = 2
    TSPAN=(t0, tf)

    θ2  = 0.0
    θ2p = 0.0
    x30 = 0.0 # DC motor
    x40 = 0.0 # ∫error
    x50 = 0.0 # ∫v (desired)
    x60 = 0.0 # ∫I^2

    # initial conditions
    x0 = [θ2, θ2p, x30, x40, x50, x60]

    # dynamic system solution
    return myode45(mod4b,TSPAN,x0,KC)

end

# statistics (ECA)
#   minimum         median             mean         maximum          std
# ECA
# 2.055877e-01   2.340816e-01    2.292241e-01    2.345938e-01    9.315045e-03
# ED
# 2.120037e-01   2.232873e-01    2.206308e-01    2.340816e-01    8.093537e-03