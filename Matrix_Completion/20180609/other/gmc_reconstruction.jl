# Animation library
include("../lib/animation.jl")

# Matrix completion technique
include("proposed.jl")

# The solution of subproblem 2 of ADMM is obtained via 
# a call to the conjugate gradient algorithm.
include("../lib/cg.jl")

### Reconstruct the input point cloud based on the proposed graph-based matrix completion algorithm
function gmc_reconstruction(reconstruction_path::String, anm::Animation, anm_matrix_rank::Int64,
						 under_sample_ratio::Float64, tau::Float64, 
						 max_admm_iters::Int64, gamma::Float64, noise::String)
			 

	# Parameters setup
	anm_matrix = arrange(anm.data, [[:coors, :frames], :vertices])
	m,n = size(anm_matrix)

	# Initialization
	convergence_curve_gmc = zeros(max_admm_iters)

	# Ground-trouth
	U,S,V = svd(anm_matrix)
	thr_S = zeros(size(S))
	thr_S[1:anm_matrix_rank] = S[1:anm_matrix_rank]
	anm_matrix_thr = U*diagm(thr_S)*V'
	
	# Subsampling matrix
	Omega=zeros(m,n)
	for f=1:3:m
		srand(f)
		s = Int(round(under_sample_ratio*n))
		Omega_indices=randperm(n)[1:s]
		Omega[f, Omega_indices]=ones(s)
		Omega[f+1, :] = Omega[f, :]
		Omega[f+2, :] = Omega[f, :]
	end
	
	
	if(noise=="noisy")
		filename = reconstruction_path*"average_graph_"*string(round(Int,under_sample_ratio*100))*".jld"		
		file = load(filename)
		L = file["data"]
	else
		filename = reconstruction_path*string("average_graph_100.txt")
		L = sparse(readdlm(filename))
	end

	anm_matrix_gmc, convergence_curve_gmc = proposed_gmc(anm_matrix, Omega, tau, max_admm_iters, gamma, L, true)

	
	# Performance evaluations
	NMSE_gmc=vecnorm(anm_matrix_gmc - anm_matrix, 2)/vecnorm(anm_matrix,2) 	
	KGerror_gmc=100*vecnorm(anm_matrix-anm_matrix_gmc, 2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	KGerror_thr=100*vecnorm(anm_matrix_thr-anm_matrix,2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	
	# I/0
	anm_gmc = deepcopy(anm)
	anm_gmc.data = LArray(reshape(anm_matrix_gmc, 3, round(Int, m/3), n), [:coors, :frames, :vertices])
	
	println("Threshold KG error:"*string(KGerror_thr))
	writedlm(reconstruction_path*string(under_sample_ratio)*"_convergence_curve_gmc.txt", convergence_curve_gmc)
	writedlm(reconstruction_path*string(under_sample_ratio)*"_KGerror_gmc.txt", KGerror_gmc)
	writedlm(reconstruction_path*string(under_sample_ratio)*"_gmc_parameters.txt", [gamma tau])
	writedlm(reconstruction_path*string(under_sample_ratio)*"_NMSE_gmc.txt", NMSE_gmc)
	export_animation_obj(anm_gmc, reconstruction_path, "gmc_"*string(under_sample_ratio)*"_frame_")
	println("KG error GMC (Subsampling ratio:"*string(round(under_sample_ratio*100))*"%):"*string(round(KGerror_gmc,4)))
	println("NMSE GMC:"*string(NMSE_gmc))
end