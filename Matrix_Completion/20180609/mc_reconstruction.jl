
# Animation library
include("../lib/animation.jl")

# Matrix completion technique
include("proposed.jl")
include("graphLaplacians.jl")

# The solution of subproblem 2 of ADMM is obtained via 
# a call to the conjugate gradient algorithm.
include("../lib/cg.jl")

### Reconstruct the input point cloud based on the proposed graph-based matrix completion algorithm
function mc_reconstruction(reconstruction_path::String, anm::Animation, anm_matrix_rank::Int64,
						 under_sample_ratio::Float64, tau::Float64, 
						 max_admm_iters::Int64)
			 

	# Parameters setup
	anm_matrix = arrange(anm.data, [[:coors, :frames], :vertices])
	m,n = size(anm_matrix)

	# Initialization
	convergence_curve_mc = zeros(max_admm_iters)

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
	
	anm_matrix_mc, convergence_curve_mc = proposed_mc(anm_matrix, Omega, tau, max_admm_iters, true)

	
	# Performance evaluations
	NMSE_mc=vecnorm(anm_matrix_mc - anm_matrix, 2)/vecnorm(anm_matrix,2) 	
	KGerror_mc=100*vecnorm(anm_matrix-anm_matrix_mc, 2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	KGerror_thr = 100*vecnorm(anm_matrix_thr-anm_matrix,2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	
	# I/0
	anm_mc = deepcopy(anm)
	anm_mc.data = LArray(reshape(anm_matrix_mc, 3, round(Int, m/3), n), [:coors, :frames, :vertices])
	
	println("Threshold KG error:"*string(KGerror_thr))
	writedlm(reconstruction_path*string(under_sample_ratio)*"_thresh_mc.txt", [anm_matrix_rank KGerror_thr])
	writedlm(reconstruction_path*string(under_sample_ratio)*"_convergence_curve_mc.txt", convergence_curve_mc)
	writedlm(reconstruction_path*string(under_sample_ratio)*"_KGerror_mc.txt", KGerror_mc)
	writedlm(reconstruction_path*string(under_sample_ratio)*"_parameters_mc.txt", [tau])
	writedlm(reconstruction_path*string(under_sample_ratio)*"_NMSE_mc.txt", NMSE_mc)
	export_animation_obj(anm_mc, reconstruction_path, "mc_"*string(under_sample_ratio)*"_frame_")
	println("KG error MC (Subsampling ratio:"*string(round(under_sample_ratio*100))*"%):"*string(round(KGerror_mc,4)))
	println("NMSE MC:"*string(NMSE_mc))
end