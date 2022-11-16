using DelimitedFiles

include("SALO.jl")


function run_example_ntw3()
	@info "Running the example script for the SALO algorithm on the 3-node example."
	
	ntw = "ntw3"
	
	###############################################################
	# Loading the time series
	###############################################################
	@info "Loading time series in 'Xs'..."
	Xs = readdlm("./example_data/"*ntw*"_Xs.csv",',')
	@info "Loaded."
		
	
	###############################################################
	# Defining the parameters
	###############################################################
	τ = readdlm("./example_data/"*ntw*"_tau.csv",',')[1]
	ks = vec(Int64.(readdlm("./example_data/"*ntw*"_ks.csv",',')))
	ls = vec(Int64.(readdlm("./example_data/"*ntw*"_ls.csv",',')))
	
	
	###############################################################
	# Ground truth
	###############################################################
	gtnode = Int64(readdlm("example_data/"*ntw*"_gtnode.csv",',')[1])
	gtfreq = readdlm("example_data/"*ntw*"_gtfreq.csv",',')[1]
	
	
	###############################################################
	# Running SALO
	###############################################################
	@info "Running SALO..."
	xxx = run_SALO(Xs,τ,ls,ks,true)
	
	@info "==========================================================================="
	@info "SALO identifed a forcing at node $(xxx[2][5]) with frequency $(xxx[2][6]/(τ*size(Xs)[2]))."
	@info "Ground truth is a forcing at node $gtnode with frequency $gtfreq."
	@info "==========================================================================="
	
	
	@info "Pause for 5 seconds."
	pause(5)
	
	###############################################################
	# Running SALO-relax
	###############################################################
	@info "Running SALO-relax..."
	xxx = run_SALOrelax(Xs,τ,ks,true)
	
	@info "==========================================================================="
	@info "SALO-relax identifed a forcing at node $(xxx[2][5]) with frequency $(xxx[2][6]/(τ*size(Xs)[2]))."
	@info "Ground truth is a forcing at node $gtnode with frequency $gtfreq."
	@info "==========================================================================="
	
	return nothing
end


function run_example_ntw20()
	@info "Running the example script for the SALO algorithm on the 20-node example."
	
	ntw = "ntw20"
	
	###############################################################
	# Loading the time series
	###############################################################
	@info "Loading time series in 'Xs'..."
	Xs = readdlm("./example_data/"*ntw*"_Xs.csv",',')
	@info "Loaded."
	
	
	###############################################################
	# Defining the parameters
	###############################################################
	τ = readdlm("./example_data/"*ntw*"_tau.csv",',')[1]
	ks = vec(Int64.(readdlm("./example_data/"*ntw*"_ks.csv",',')))
	ls = vec(Int64.(readdlm("./example_data/"*ntw*"_ls.csv",',')))
	
	
	###############################################################
	# Ground truth
	###############################################################
	gtnode = Int64(readdlm("example_data/"*ntw*"_gtnode.csv",',')[1])
	gtfreq = readdlm("example_data/"*ntw*"_gtfreq.csv",',')[1]
	
	
	###############################################################
	# Running SALO
	###############################################################
	@info "Running SALO..."
	xxx = run_SALO(Xs,τ,ls,ks,true)
	
	@info "==========================================================================="
	@info "SALO identifed a forcing at node $(xxx[2][5]) with frequency $(xxx[2][6]/(τ*size(Xs)[2]))."
	@info "Ground truth is a forcing at node $gtnode with frequency $gtfreq."
	@info "==========================================================================="
	
	
	@info "Pause for 5 seconds."
	pause(5)
	
	###############################################################
	# Running SALO-relax
	###############################################################
	@info "Running SALO-relax..."
	
	xxx = run_SALOrelax(Xs,τ,ks,true)
	
	@info "==========================================================================="
	@info "SALO-relax identifed a forcing at node $(xxx[2][5]) with frequency $(xxx[2][6]/(τ*size(Xs)[2]))."
	@info "Ground truth is a forcing at node $gtnode with frequency $gtfreq."
	@info "==========================================================================="

	return nothing
end


