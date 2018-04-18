using Metaheuristics
using Mechanisms

function findStructure(ncase)
    if ncase == 1
        precisionpts = PTS_VERTICAL_LINE

        D, f, g = getOptimizationProblem(precisionpts)

        bounds =
            [ 0 60;
              0 60;
              0 60;
              0 60;
            -60 60;
            -60 60.0;
              0 2π;
            -60 60;
            -60 60;
             repmat([0 2π], size(precisionpts, 1))
            ]'

        K = 7; N = K*D; η_max = 2
    elseif ncase == 2
        precisionpts = PTS_ELLIPTIC

        # synchronization
        X0 = [0,0.0,0]
        θ2 = [π/6, π/4, π/3, 10π/24, π/2]
        D, f, g = getOptimizationProblem(precisionpts, X0, θ2)

        bounds =
            [ 0 60;
              0 60;
              0 60;
              0 60;
            -60 60;
            -60 60.0;
            ]'

        K = 7; N = K*D; η_max = 2
    else

        D, f, g = getOptimizationProblem(PTS_PAIR_1, PTS_PAIR_2)
        precisionpts = PTS_PAIR_1

        bounds =
            [ 0 60;
              0 60;
              0 60;
              0 60;
            -60 60;
            -60 60.0;
              0 2π;
            -60 60;
            -60 60;
             repmat([0 2π], size(precisionpts, 1))
            ]'

        precisionpts = [PTS_PAIR_1 ; PTS_PAIR_2]
        K = 7; N = K*D; η_max = 2
    end

    objFunction(p, α = 1000.0) = f(p) + α*sum(g(p))
    
    p, err = eca(objFunction, D;
                    η_max = η_max,
                    K     = K,
                    N     = N,
                    limits= bounds,
                    p_bin= 0.03,
                    max_evals=20000D)
    if sum(g(p)) > 0
        warn("No feasible.")
    end

    p, err, precisionpts
end