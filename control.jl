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
                 saveLast = "gen.csv",
                 showIter = false,
                max_evals = 5000)

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
#   minimum         mean             median         maximum          std
# 2.149911e-01   2.307464e-01    2.340816e-01    2.345938e-01    6.313997e-03