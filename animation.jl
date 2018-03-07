using Plots
gr(size=(500, 500), legend=false)

function getAnimation(C, X0, X1, X2, X3, nframes, my_error, precision_points)
    @gif for t = 1:nframes
        plot(title=@sprintf("Error: %e", my_error))
        scatter!(precision_points[:,1], precision_points[:,2])
        plot!(C[:,1], C[:,2], linestyle=:dot, linecolor=:green)
        
        scatter!(X0[t, 1:1], X0[t, 2:2])
        scatter!(X1[t, 1:1], X1[t, 2:2])
        scatter!(X2[t, 1:1], X2[t, 2:2])
        scatter!(X3[t, 1:1], X3[t, 2:2])
        
        scatter!(C[t, 1:1], C[t, 2:2], markercolor=:red)
       

        plot!([X0[t, 1] , X2[t, 1] ], [X0[t, 2] , X2[t, 2] ])
        plot!([X1[t, 1] , X0[t, 1] ], [X1[t, 2] , X0[t, 2] ])
        plot!([X1[t, 1] , X3[t, 1] ], [X1[t, 2] , X3[t, 2] ])
        
        plot!([X2[t, 1] , X3[t, 1] ], [X2[t, 2] , X3[t, 2] ], linecolor=:red)
        plot!([X2[t, 1] , C[t, 1]], [X2[t, 2] , C[t, 2]], linecolor=:red)
        plot!([X3[t, 1] , C[t, 1]], [X3[t, 2] , C[t, 2]], linecolor=:red) 
    end
end