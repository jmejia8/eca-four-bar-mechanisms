using Plots
plotly( size=(1000, 500) )

function plotTuning(KC, K_m)
	t, u, θ2       = modelovolt(KC)
	t_m, u_m, θ2_m = modelovolt(K_m)

	p1 = plot(  t,  θ2, ylimits=(25, 31), linecolor=:blue, label="Optimum tuning")
	plot!(t_m, θ2_m, linecolor=:red, linestyle=:dot , label="Trial and error tuning")

	p2 = plot(  t,  u, linecolor=:blue, label="Optimum tuning")
	plot!(t_m, u_m, linecolor=:red, linestyle=:dot , label="Trial and error tuning")

	plot(p1, p2)
end