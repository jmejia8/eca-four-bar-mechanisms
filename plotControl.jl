using PyPlot
include("control.jl")

function plotTuning(KC, K_m)
    t, u, θ2       = modelovolt(KC)
    t_m, u_m, θ2_m = modelovolt(K_m)

    subplot(1,2,1)
    grid("on")
    plot(  t,  θ2, "b")
    plot(t_m, θ2_m, "r--")
    legend(["Optimum tuning", "Trial and error tuning"], loc=4)
    ylim([25, 31])
    xlabel("t [s]")
    ylabel("\$\\theta_2\$ [rad / s]")

    subplot(1,2,2)
    grid("on")
    plot(  t,  u, "b")
    plot(t_m, u_m, "r--" )
    legend(["Optimum tuning", "Trial and error tuning"], loc=4)
    xlabel("t [s]")
    ylabel("\$ u(t) \$ [V]")

end

function plotConvergence()
    figure(figsize=(16, 4), dpi=80)

    fname = "output/conv_control_median.csv"
    case1 = readcsv(fname)

    # subplot(1,3,1)
    # title("Study Case 1")
    plot(case1[:,1], log.(1+case1[:,2]), "k", lw=2)
    xlabel("Evaluations")
    ylabel("Log Error")
    grid("on")
end

function main()
    K, F = getKC()
    K_m = [45.55, 5.25, 1]

    # while F > 2.05588e-01
    #     K, F = getKC()
    #     @printf("%e ==== ", F)
    #     println(K)
    # end

    plotTuning(K, K_m)
end

# main()
plotConvergence()