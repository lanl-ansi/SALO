using Distributed

# Define the number of parallel threads.
n_thr = 3
@info "How many parallel threads are required? (Default is 3)"
try
	global n_thr = parse(Int64,readline())
catch eee
	@info "Default number of threads used: 3."
end

@info "Loading optimizer..."

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end

# If the directory "data" does not exists, then create it.
if !isdir("data/")
	mkdir("data/")
end

@everywhere include("SALO.jl")

#=
"""
	run_SALO_par(id::String, Xs::Matrix{Float64}, τ::Float64, ls::Vector{Int64}, ks::Vector{Int64}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized verions of "run_SALO".

The parameter "id" identifies the system to be identified.
"""
=#
@everywhere function run_SALO_par(id::String, Xs::Matrix{Float64}, τ::Float64, ls::Vector{Int64}, ks::Vector{Int64}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/τ
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

	A1h = zeros(n,n)
	a2h = ones(n)

	args = Array{Tuple{String,Matrix{Float64},Matrix{Float64},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Matrix{Float64},Vector{Float64},Float64,Float64,Float64},1}()

	for l in ls
		for k in ks
			push!(args,(id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp))
		end
	end

	pmap(Lmax_SALO_par,args)
end


#=
"""
	Lmax_SALO_par(id::String, x::Matrix{Float64}, Dx::Matrix{Float64}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Matrix{Float64}, a1h::Vector{Float64}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of "Lmax_SALO", i.e., stores the data in files.
"""
=#
@everywhere function Lmax_SALO_par(tups::Tuple{String,Matrix{Float64},Matrix{Float64},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Int64,Matrix{Float64},Vector{Float64},Float64,Float64,Float64})
	id,x,Dx,xt,Dxt,l,k,A1h,a2h,b,μ,bp = tups

	@info "===================================================================================="
	@info "Computing "*id*" with SALO optimization: ℓ = $l, k = $k."
	@info "===================================================================================="
	
	nn,N = size(x)
	n = Int(nn/2)

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxt[l,k+1]*xtk')

	glk = norm(Dxt[l,k+1])^2

# Definition of the optimization problem.
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => μ, "bound_push" => bp))
	@variable(system_id, A1[i = 1:n, j = 1:n])
	for i = 1:n
		for j = 1:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end
	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@constraint(system_id, c1[i=1:n], a2[i] >= 0.)
	@variable(system_id, γ >= 0.)
	set_start_value(γ,1.)

# tr(L^2*θ*θ^T)
	@NLexpression(system_id, AtAS01, 
		      sum(sum(sum(A1[j,i]*A1[j,m]*S0[m,i] for i = 1:n) for j = 1:n) for m = 1:n))

# tr(-L*D*ω*θ^T) = tr(-D*L*θ*ω^T)
	@NLexpression(system_id, AtAS02a, 
		      sum(sum(A1[j,i]*a2[j]*S0[n+j,i] for i = 1:n) for j = 1:n))
	@NLexpression(system_id, AtAS02b, 
		      sum(sum(a2[i]*A1[i,j]*S0[j,n+i] for i = 1:n) for j = 1:n))

# tr(D^2*ω*ω^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + AtAS02a + AtAS02b + AtAS03)


# tr(-L*θ*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(sum(A1[i,j]*S1[j,i] for i = 1:n) for j = 1:n))

# tr(-D*ω*(Dω)^T)
	@NLexpression(system_id, AS12, sum(a2[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1, 
		      sum(sum(A1[l,i]*A1[l,j]*Fk[j,i] for i = 1:n) for j = 1:n))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2a, 
		      sum(A1[l,i]*a2[l]*Fk[n+l,i] for i = 1:n))
	@NLexpression(system_id, AAF2b, 
		      sum(a2[l]*A1[l,j]*Fk[j,n+l] for j = 1:n))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3, a2[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3, AAF1 + AAF2a + AAF2b + AAF3)

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4, sum(flk[i]*A1[l,i] for i = 1:n) + flk[n+l]*a2[l])


	@NLexpression(system_id, γ2, γ^2)

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2*γ/sqrt(N)*sqrt(T3 + 2*T4 + glk) + b*sum(sum(abs(A1[i,j])+abs(A1[j,i]) for j = i+1:n) for i = 1:n-1))

	optimize!(system_id)

	mL = zeros(n,n)
	for i in 1:n
		for j in 1:n
			mL[i,j] = value(A1[i,j])
		end
	end

	writedlm("data/"*id*"_SALO_$(l).$(k)_obj.csv",-objective_value(system_id),',')
	writedlm("data/"*id*"_SALO_$(l).$(k)_A1.csv",mL,',')
	writedlm("data/"*id*"_SALO_$(l).$(k)_a2.csv",value.(a2),',')
	writedlm("data/"*id*"_SALO_$(l).$(k)_g.csv",value(γ),',')
end

#=
"""
	run_SALOrelax_par(Xs::Matrix{Float64}, τ::Float64, ks::Vector{Int64}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of "run_SALOrelax".

The parameter "id" identifies the system to be identified.
"""
=#
@everywhere function run_SALOrelax_par(id::String, Xs::Matrix{Float64}, τ::Float64, ks::Vector{Int64}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

# Computing the needed inputs (time series, discrete derivative, and their Fourier transforms).
	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/τ
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

	A1h = zeros(n,n)
	a2h = ones(n)

# Run the optimizations
	args = Array{Tuple{String,Matrix{Float64},Matrix{Float64},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Matrix{Float64},Vector{Float64},Tuple{Float64,Float64},Float64,Float64},1}()

	for k in ks
		push!(args,(id,x,Dx,xt,Dxt,k,A1h,a2h,b,μ,bp))
	end

	pmap(Lmax_SALOrelax_par,args)
end

#=
"""
	Lmax_SALOrelax_par(id::String, x::Matrix{Float64}, Dx::Matrix{Float64}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Matrix{Float64}, a2h::Vector{Float64}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)

Parallelized version of Lmax_SALOrelax.
"""
=#
@everywhere function Lmax_SALOrelax_par(tups::Tuple{String,Matrix{Float64},Matrix{Float64},Array{Complex{Float64},2},Array{Complex{Float64},2},Int64,Matrix{Float64},Vector{Float64},Tuple{Float64,Float64},Float64,Float64})
	id,x,Dx,xt,Dxt,k,A1h,a2h,b,μ,bp = tups
	
	@info "===================================================================================="
	@info "Computing "*id*" with SALOrelax optimization: k = $k."
	@info "===================================================================================="
	
	nn,N = size(x)
	n = Int(nn/2)
	b1,b2 = b

# Definition of the needed parameters.
	S0 = (x*x')/N
	S1 = (x*Dx')/N

	xtk = xt[:,k+1]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k+1]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k+1])^2 for l = 1:n]

# Definition of the optimization problem. 
	system_id = Model(optimizer_with_attributes(Ipopt.Optimizer, "mu_init" => μ, "bound_push" => bp))

	@variable(system_id, A1[i = 1:n, j = 1:n])
	for i in 1:n
		for j in 1:n
			set_start_value(A1[i,j],A1h[i,j])
		end
	end

	@variable(system_id, a2[i = 1:n])
	set_start_value.(a2,a2h)
	@constraint(system_id, c1[i=1:n], a2[i] >= 0.)

	@variable(system_id, γ[l = 1:n])
	for l in 1:n
		@constraint(system_id, γ[l] >= 0.)
	end
	set_start_value.(γ,ones(n))

# tr(L^2*θ*θ^T)
	@NLexpression(system_id, AtAS01, 
		      sum(sum(sum(A1[j,i]*A1[j,m]*S0[m,i] for i = 1:n) for j = 1:n) for m = 1:n))

# tr(-L*D*ω*θ^T) = tr(-D*L*θ*ω^T)
	@NLexpression(system_id, AtAS02a, 
		      sum(sum(A1[j,i]*a2[j]*S0[n+j,i] for i = 1:n) for j = 1:n))
	@NLexpression(system_id, AtAS02b, 
		      sum(sum(a2[i]*A1[i,j]*S0[j,n+i] for i = 1:n) for j = 1:n))

# tr(D^2*ω*ω^T)
	@NLexpression(system_id, AtAS03, sum(a2[i]^2*S0[n+i,n+i] for i = 1:n))
	
# tr(A^T*A*S0)
	@NLexpression(system_id, T1, AtAS01 + AtAS02a + AtAS02b + AtAS03)


# tr(-L*θ*(Dω)^T)
	@NLexpression(system_id, AS11, 
		      sum(sum(A1[i,j]*S1[j,i] for i = 1:n) for j = 1:n))

# tr(-D*ω*(Dω)^T)
	@NLexpression(system_id, AS12, sum(a2[i]*S1[n+i,i] for i = 1:n))

# tr(A*S1)
	@NLexpression(system_id, T2, AS11 + AS12)


# tr((L_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF1[l = 1:n], 
		      sum(sum(A1[l,i]*A1[l,j]*Fk[j,i] for i = 1:n) for j = 1:n))

# tr((L_{l,:})^T*D_{l,:}*F(k)) = tr((D_{l,:})^T*L_{l,:}*F(k))
	@NLexpression(system_id, AAF2a[l = 1:n], 
		      sum(A1[l,i]*a2[l]*Fk[n+l,i] for i = 1:n))
	@NLexpression(system_id, AAF2b[l = 1:n], 
		      sum(a2[l]*A1[l,j]*Fk[j,n+l] for j = 1:n))

# tr((D_{l,:})^T*D_{l,:}*F(k))
	@NLexpression(system_id, AAF3[l = 1:n], a2[l]^2*Fk[n+l,n+l])

# tr((A_{l,:})^T*A_{l,:}*F(k))
	@NLexpression(system_id, T3[l = 1:n], AAF1[l] + AAF2a[l] + AAF2b[l] + AAF3[l])

	
# (f_l(k))^T*A_{l,:}
	@NLexpression(system_id, T4[l = 1:n], sum(flk[l][i]*A1[l,i] for i = 1:n) + flk[l][n+l]*a2[l])


	@NLexpression(system_id, γ2, sum(γ[l]^2 for l = 1:n))

	
	@NLobjective(system_id, Min, T1 + 2*T2 + .5*γ2 - 2/sqrt(N)*sum(γ[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n) + b1*sum(sum(abs(A1[i,j])+abs(A1[j,i]) for j in i+1:n) for i = 1:n-1) + b2*sum(γ[l] for l in 1:n))

	optimize!(system_id)

	mL = zeros(n,n)
	for i in 1:n
		for j in 1:n
			mL[i,j] = value(A1[i,j])
		end
	end
	
	writedlm("data/"*id*"_SALOr_$(k)_obj.csv",objective_value(system_id),',')
	writedlm("data/"*id*"_SALOr_$(k)_A1.csv",mL,',')
	writedlm("data/"*id*"_SALOr_$(k)_a2.csv",value.(a2),',')
	writedlm("data/"*id*"_SALOr_$(k)_g.csv",value.(γ),',')
end

@info "Loaded."


