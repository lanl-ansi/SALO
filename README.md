# System Agnostic Localization of Oscillations (SALO)
Implementation of the *System Agnostic Localization of Oscillations (SALO)* algorithm, presented in the manuscript [Delabays, Lokhov, Tyloo, and Vuffray (2022)](https://arxiv.org). 
The SALO algorithm identifies the source and frequency of a forced oscillation in a complex network of coupled dynamical agents, based on the position and velocity time series of each agent. 


## System requirements
The code has been developped on Julia 1.6 and should work on any subsequent version, with up-to-date packages. The packages requires are:

- DelimitedFiles
- Distributed
- FFTW
- Ipopt
- JuMP
- LinearAlgebra
- PyPlot


## Summary of the files and folders
- **example\_data**: Contains the data necessary to run the two examples. 
- **plot\_SALO.jl**: Loads the scripts plotting the outcome of the SALO algorithm.
- **run\_examples.jl**: Loads the functions running the examples for the SALO algorithms. 
- **SALO.jl**: Loads the SALO and SALO-relaxed algorithms.
- **SALO\_parallel.jl** Loads the parallelized versions of the SALO and SALO-relaxed algorithms.


## List of functions
- [Lmax\_SALO (Lmax\_SALO\_par)](#LmaxSALO)
- [Lmax\_SALOrelax (Lmax\_SALOrelax\_par)](#LmaxSALOr)
- [plot\_SALO](#plot\_SALO)
- [plot\_SALOrelax\_1](#plot\_SALOrelax\_1)
- [plot\_SALOrelax\_2](#plot\_SALOrelax\_2)
- [run\_example\_ntw20](#run\_example\_ntw20)
- [run\_example\_ntw3](#run\_example\_ntw3) 
- [run\_SALO (run\_SALO\_par)](#runSALO)
- [run\_SALOrelax (run\_SALOrelax\_par)](#runSALOr)


## Detailed documentation for the functions


### <a id="LmaxSALO"></a> Lmax\_SALO (Lmax\_SALO\_par)
*./SALO.jl (./SALO\_parallel.jl)*

- `Lmax_SALO(x::Matrix{Float64}, Dx::Matrix{Float64}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Matrix{Float64}, a1h::Vector{Float64}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)`
- `Lmax_SALO_par(id::String, x::Matrix{Float64}, Dx::Matrix{Float64}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, l::Int64, k::Int64, A1h::Matrix{Float64}, a1h::Vector{Float64}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)`

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) and location (l) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (γ). 

**INPUT**:\\
(`ntw`: For `Lmax_SALO_par` only. Name of the system under investigation, for data labelling purpose.)\
`x`: Time series of the phase angles.\
`Dx`: Time series of the phase frequencies. \
`xt`: (Inverse) Fourier transform of x.\
`Dxt`: (Inverse) Fourier transform of Dx.\
`l`: Fixed location of the forcing.\
`k`: Fixed index of the forcing frequency (f = k/T).\
`A1h`: Warm start for A1.\
`a2h`: Warm start for a2.\
`b`: Regularization parameter to avoid overfitting.\
`μ`: Initial value of the barrier parameter (in IPOPT).\
`bp`: Initial value of the bound_push parameter (in IPOPT).


**OUTPUT**:\
`objective`: Value of the optimized objective function.\
`A1`: Best estimate of the dynamcis matrix.\
`a2`: Best estimate of the dampings.\
`γ`: Best estimate of the forcing amplitude.


---

### <a id="LmaxSALOr"></a> Lmax\_SALOrelax (Lmax\_SALOrelax\_par)
*./SALO.jl (./SALO\_parallel.jl)*


 - `Lmax_SALOrelax(x::Matrix{Float64}, Dx::Matrix{Float64}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Matrix{Float64}, a2h::Vector{Float64}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)`
- `Lmax_SALOrelax_par(id::String, x::Matrix{Float64}, Dx::Matrix{Float64}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64,2}, k::Int64, A1h::Matrix{Float64}, a2h::Vector{Float64}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)`

Minimizes the quadratic error in the estimation of the forced trajectory, for a fixed frequency (k) of the forcing. The optimization parameters are the dynamics matrix (A1), the damings (a2), and the forcing amplitude (γ). 


**INPUT**:
(`ntw`: For `Lmax_SALOrelax_par` only. Name of the system under investigation, for data labelling purpose.)\
`x`: Time series of the phase angles.\
`Dx`: Time series of the phase frequencies.\
`xt`: (Inverse) Fourier transform of x.\
`Dxt`: (Inverse) Fourier transform of Dx.\
`k`: Fixed index of the forcing frequency (f = k/T).\
`A1h`: Warm start for A1.\
`a2h`: Warm starrt for a2.\
`b`: Regularizaion parameter to avoid overfitting.\
`μ`: Initial value of the barrier parameter (in IPOPT).\
`bp`: Initial value of the bound_push parameter (in IPOPT).


**OUTPUT**:\
`objective`: Value of the optimized objective function.\
`A1`: Best estimate of the dynamcis matrix.\
`a2`: Best estimate of the dampings.\
`γ`: Best estimate of the forcing amplitude at each possible location.


---


### plot\_SALO
*./plot\_SALO.jl*

- `plot_SALO(Ls::Matrix{Float64}, L::Float64, l::Int64, k::Int64, ks::Vector{Int64}, T::Union{Float64,Int64}, n::Int64)`

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the SALO approach.

**INPUT**:\
`Ls`: Values of the objective for the various values of l (rows) and k (columns).\
`L`: Minimal objective value.\
`l0`: Corresponding value of l.\
`k0`: Corresponding value of k.\
`ks`: Values of k assessed.\
`T`: Observation time.\
`n`: Number of agents.


---

### plot\_SALOrelax\_1
*./plot\_SALO.jl*

- `plot_SALOrelax_1(Ls::Array{Float64,1}, L::Float64, γ::Array{Float64,1}, l0::Int64, k0::Int64, ks::Array{Int64,1}, T::Union{Float64,Int64})`

Plots the value of the objective vs. the estimated (discretized) forcing frequencies by the SALO-relax approach.

**INPUT**:\
`Ls`: Values of the objective for the various values of k.\
`L`: Minimal objective value.\
`k0`: Corresponding value of k.\
`ks`: Values of assessed.


---

### plot\_SALOrelax\_2
*./plot\_SALO.jl*

- `plot_SALOrelax_2(Ls::Vector{Float64}, L::Float64, γ::Array{Float64,1}, l0::Int64)`

Plots the value of the estimated forcing amplitude vs. the agent index, obtained by the SALO-relax approach.

**INPUT**:\
`Ls`: Values of the objective for the various values of k.\
`L`: Minimal objective value.\
`γ`: Corresponding forcing amplitudes.\
`l0`: Index with largest amplitude.


---

### run\_example\_ntw20
*./run\_examples.jl*

- `run_example_ntw20()`

Runs the SALO and SALO-relax algorithms on a 20-node system. All data supporting the examples can be found in the folder `example\_data`.


---

### run\_example\_ntw3
*./run\_examples.jl*

- `run_example_ntw3()`

Runs the SALO and SALO-relax algorithms on a 3-node system. All data supporting the examples can be found in the folder `example\_data`.

---

### <a id="runSALO"></a> run\_SALO (run\_SALO\_par)
*./SALO.jl (./SALO_parallel.jl)*

- `run_SALO(Xs::Matrix{Float64}, τ::Float64, ls::Vector{Int64}, ks::Vector{Int64}, plot::Bool=false, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)`
- `run_SALO_par(id::String, Xs::Matrix{Float64}, τ::Float64, ls::Vector{Int64}, ks::Vector{Int64}, b::Float64=0., μ::Float64=1e-1, bp::Float64=1e-1)`

Runs the SALO algorithm on the time series `Xs`, with candidate sources' indices in `ls` and forcing's candidate frequencies indices in `ks`. Returns the Maximum Likelihood for each possible pair (l,k) as well as the identified system parameters corresponding to the largest likelihood.

**INPUT**:\
(`id`: For the parallel version only. Name of the system under investigation, for data labelling purpose.)\
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).\
`τ`: Time step size.\
`ls`: Values of l to be tried.\
`ks`: Values of k to be tried.\
`plot`: For the non-parallel version only. If true, generates the plots of the objective function vs. frequency.\
`b`: Regularization parameter to avoid overfitting.\
`μ`: Initial value of the barrier parameter (in IPOPT).\
`bp`: Initial value of the bound_push parameter (in IPOPT).

**OUTPUT (written in folder ./data/ for the parallel version)**:\
`Ls_l0`: values of the objective function for the various values of k and l, for the SALO algorithm.\
`(L_l0, A_l0, d_l0, γ_l0, k_l0, l_l0)`: results under l0.\
        `L_l0`: Minimal value of the objective.\
        `A_l0`: Estimate of the dynamics matrix.\
        `d_l0`: Estimate of the dampings.\
        `γ_l0`: Estimate of the forcing amplitude.\
        `k_l0`: Estimate frequency index (see theory).\
        `l_l0`: Estimate of the forcing location.

---

### <a id="runSALOr"></a> run\_SALOrelax (run\_SALOrelax\_par)
*./SALO.jl (./SALO\_parallel.jl)*

- `run_SALOrelax(Xs::Matrix{Float64}, τ::Float64, ks::Vector{Int64}, plot::Bool=false, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)`
- `run_SALOrelax_par(Xs::Matrix{Float64}, τ::Float64, ks::Vector{Int64}, b::Tuple{Float64,Float64}=(0.,0.), μ::Float64=1e-1, bp::Float64=1e-1)`

Runs the SALO-relax algorithm on the time series `Xs`, with forcing's candidate frequencies indices in `ks`. Returns the Maximum Likelihood for each possible k as well as the identified system parameters corresponding to the largest likelihood.

**INPUT**:\
(`ntw`: For parallel version only. Name of the system under investigation, for data labelling purpose.)\
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (rows n+1:2*n).\
`τ`: Time step size.\
`ks`: Values of k to try.\
`plot`: For the non-parallel version only. If true, generates the plots of the objective function vs. k.\
`b`: Regularization parameters to avoid overfitting.\
`μ`: Initial value of the barrier parameter (in IPOPT).\
`bp`: Initial value of the bound\_push parameter (in IPOPT).

**OUTPUT (written in folder ./data/ for the parallel version)**:\
`Ls_l1`: values of the objective function for the various values of k and l, under l1.\
`(L_l1, A_l1, d_l1, γ_l1, k_l1, l_l1)`: results under l1.\
        `L_l1`: Minimal value of the objective.\
        `A_l1`: Estimate of the dynamics matrix.\
        `d_l1`: Estimate of the damings.\
        `γ_l1`: Estimate of the forcing amplitude.\
        `l_l1`: Estimate of the forcing location.\
        `k_l1`: Estimate frequency index (see theory).


---



