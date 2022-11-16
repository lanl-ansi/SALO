using PyPlot

"""
	plot_SALO(Ls::Matrix{Float64}, L::Float64, l::Int64, k::Int64, ks::Vector{Int64}, T::Union{Float64,Int64}, n::Int64)

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the SALO approach.

_INPUT_:
`Ls`: Values of the objective for the various values of l (rows) and k (columns). 
`L`: Minimal objective value.
`l0`: Corresponding value of l.
`k0`: Corresponding value of k.
`ks`: Values of k assessed.
`T`: Observation time.
`n`: Number of agents.
"""
function plot_SALO(Ls::Matrix{Float64}, L::Float64, l0::Int64, k0::Int64, ks::Vector{Int64}, T::Union{Float64,Int64}, n::Int64)
	fs = ks/T
	gra = (.8,.8,.8,1.)

	if n > 5
		l1 = findmax(Ls)[2][1]
		l2 = findmax([Ls[1:l1-1,:];-1000*ones(1,length(ks));Ls[l1+1:n,:]])[2][1]
		if l1 > l2
			lt = l1
			l1 = l2
			l2 = lt
		end
		l3 = findmax([Ls[1:l1-1,:];-1e6*ones(1,length(ks));Ls[l1+1:l2-1,:];-1e6*ones(1,length(ks));Ls[l2+1:n,:]])[2][1]
		ll = sort([l1,l2,l3])
	
		Lred = Ls[setdiff(1:n,ll),:]
		Lrmax = [maximum(Lred[:,i]) for i in 1:length(ks)]
		Lrmin = [minimum(Lred[:,i]) for i in 1:length(ks)]
	
		PyPlot.fill([ks;ks[end:-1:1]]./T,[Lrmax;Lrmin[end:-1:1]],color=gra)
		
		nn = 3
	else
		nn = n
	end

	for l in 1:nn
		PyPlot.plot(fs,Ls[ll[l],:],label="l = $l")
	end
	xlabel("f")
	ylabel("Log-likelihood")
	title("Best SALO: l = $(l0), f = $(round(k0/T,sigdigits=5)), L = $(round(L,sigdigits=5))")
end

"""
	plot_SALOrelax_1(Ls::Array{Float64,1}, L::Float64, γ::Array{Float64,1}, l0::Int64, k0::Int64, ks::Array{Int64,1}, T::Union{Float64,Int64})

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the SALO-relax approach.

_INPUT_:
`Ls`: Values of the objective for the various values of k. 
`L`: Minimal objective value.
`k0`: Corresponding value of k.
`ks`: Values of assessed.
"""
function plot_SALOrelax_1(Ls::Vector{Float64}, L::Float64, k0::Int64, ks::Array{Int64,1}, T::Union{Float64,Int64})
	fs = ks/T

	PyPlot.plot(fs,Ls)
	xlabel("f")
	ylabel("Log-likelihood")
	title("Best SALO-relax: f = $(round(k0/T,sigdigits=5)), L = $(round(L,sigdigits=5))")
end

"""
	plot_SALOrelax_2(Ls::Vector{Float64}, L::Float64, γ::Array{Float64,1}, l0::Int64)

Plots the value of the estimated forcing amplitude vs. the agent index, obtained by the SALO-relax approach.

_INPUT_:
`Ls`: Values of the objective for the various values of k. 
`L`: Minimal objective value.
`γ`: Corresponding forcing amplitudes.
`l0`: Index with largest amplitude.
"""
function plot_SALOrelax_2(Ls::Vector{Float64}, L::Float64, γ::Array{Float64,1}, l0::Int64)
	PyPlot.plot(1:length(γ),γ,"o")
	xlabel("node index")
	ylabel("Estimated amplitude")
	title("Source: l = $(l0), γ_l = $(round(γ[l0],sigdigits=4))")
end


