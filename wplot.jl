using PyPlot
include("findStructure.jl")

function readStructure(ncase)
    D, bounds, precision_points, error_func = case_info(ncase)

    # read saved cases (last generation of mechanism)
    if ncase == -1
        data = readcsv("output/case1_best.csv")
    elseif ncase == -3
        data = readcsv("output/case3_best.csv")
    else
        return findStructure(ncase)
    end

    errs = zeros(size(data,1))
    
    for i = 1:length(errs) 
        errs[i] = myerror(data[i,:], error_func, precision_points)
    end
    
    best = indmin(errs)
    
    return data[best,:], errs[best]
end

function plotConvergence()
    figure(figsize=(16, 4), dpi=80)

    fname = "output/conv_median_case1_final.csv"
    case1 = readcsv(fname)

    fname = "output/conv_median_case2_final.csv"
    case2 = readcsv(fname)

    fname = "output/conv_median_case3_final.csv"
    case3 = readcsv(fname)

    subplot(1,3,1)
    title("Study Case 1")
    plot(case1[:,1], log.(1+case1[:,2]), "k", lw=2)
    xlabel("Evaluations")
    ylabel("Log Error")
    grid("on")

    subplot(1,3,2)
    title("Study Case 2")
    plot(case2[:,1], log.(1+case2[:,2]), "k", lw=2)
    xlabel("Evaluations")
    ylabel("Log Error")
    xlim([0,60000])
    grid("on")

    subplot(1,3,3)
    title("Study Case 3")
    plot(case3[:,1], log.(1+case3[:,2]), "k", lw=2)
    xlabel("Evaluations")
    ylabel("Log Error")
    xlim([0,380000])
    grid("on")
end

# plotConvergence()

function wplot()
    ncase   = 1
    nframes = 50
    
    # case information
    D, bounds, precision_points, error_func = case_info(ncase)

    # get optimal mechanism
    p, my_error = readStructure(ncase)
    println(p)


    C, X0, X1, X2, X3 = generateTrayectory(p, nframes)

    t = rand(1:nframes, 1)[1]


    # pyplot settings
    figure(figsize=(8, 8), dpi=80)
    title(@sprintf("Error: %e", my_error))

    # limits
    if ncase==1
        t = div(nframes, 4)
        xlim([-30, 25])
        ylim([10, 65])
    elseif ncase==2
        t = div(nframes, 5)
        xlim([-5, 15])
        ylim([-5, 15])
    elseif ncase==3
        t = div(nframes, 4)
        xlim([-2, 3])
        ylim([-2, 3])
    end

    # axis
    plot(zeros(401), -200:200, color="gray")
    plot(-200:200, zeros(401), color="gray")

    # paths
    plot(C[:,1], C[:,2], "b--")
    plot(X2[:,1], X2[:,2], "g--")
    
    # points
    plot(precision_points[:,1], precision_points[:,2], "ko")
    
    # bars
    plot([X0[t:t, 1] , X2[t:t, 1] ], [X0[t:t, 2] , X2[t:t, 2] ], "g", lw=2)
    plot([X1[t:t, 1] , X0[t:t, 1] ], [X1[t:t, 2] , X0[t:t, 2] ], "b", lw=2)
    plot([X1[t:t, 1] , X3[t:t, 1] ], [X1[t:t, 2] , X3[t:t, 2] ], "b", lw=2)
    plot([X2[t:t, 1] , X3[t:t, 1] ], [X2[t:t, 2] , X3[t:t, 2] ], "b", lw=2)
    plot([X2[t:t, 1] , C[t:t, 1]], [X2[t:t, 2] , C[t:t, 2]], "r", lw=2)
    plot([X3[t:t, 1] , C[t:t, 1]], [X3[t:t, 2] , C[t:t, 2]], "r", lw=2)

    # links
    plot(X0[t, 1:1], X0[t, 2:2], "bo")
    plot(X1[t, 1:1], X1[t, 2:2], "bo")
    plot(X2[t, 1:1], X2[t, 2:2], "go")
    plot(X3[t, 1:1], X3[t, 2:2], "bo")

    # coupler
    plot(C[t, 1:1], C[t, 2:2], "ro")

    grid("on")


end

wplot()