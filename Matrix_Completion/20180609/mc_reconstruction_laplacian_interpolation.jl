
# Animation library
include("../lib/animation.jl")

# Matrix completion technique
include("proposed.jl")
include("graphLaplacians.jl")

# The solution of subproblem 2 of ADMM is obtained via
# a call to the conjugate gradient algorithm.
include("../lib/cg.jl")

### Reconstruct the input point cloud based on the proposed graph-based matrix completion algorithm
function mc_reconstruction_laplacian_interpolation(
	matrix::Matrix{Float64},
	Omega::Matrix{Float64},
	L1::SparseMatrixCSC{Float64,Int64},
	#L2::SparseMatrixCSC{Float64,Int64},
	matrix_rank::Int64,
	under_sample_ratio::Float64,
	tau::Float64,
	max_admm_iters::Int64,
	verbose::Bool,
	gamma::Float64
	)
	m, n = size(matrix)

	# Initialization
	convergence_curve_mc = zeros(max_admm_iters)

	# Ground-trouth
	U,S,V = svd(matrix)
	thr_S = zeros(size(S))
	thr_S[1:matrix_rank] = S[1:matrix_rank]
	matrix_thr = U*diagm(thr_S)*V'



    #matrix_mc, convergence_curve_mc,current_iter,error_between_iterations = proposed_gmc_v2(matrix, Omega, tau, max_admm_iters,gamma,L1,L2,verbose)
	matrix_mc, convergence_curve_mc,current_iter,error_between_iterations = proposed_gmc_v2(matrix, Omega, tau, max_admm_iters,gamma,L1,verbose)

    # println(current_iter)
	# Performance evaluations
	NMSE_mc=vecnorm(matrix_mc - matrix, 2)/vecnorm(matrix,2)
	KGerror_mc=100*vecnorm(matrix-matrix_mc, 2)/vecnorm(matrix.-mean(matrix, 2), 2)
	KGerror_thr = 100*vecnorm(matrix_thr-matrix,2)/vecnorm(matrix.-mean(matrix, 2), 2)

	# I/0
	#mc = deepcopy(anm)
	#mc.data = LArray(reshape(matrix_mc, 3, round(Int, m/3), n), [:coors, :frames, :vertices])



	#writedlm("thresh_mc.txt", [matrix_rank KGerror_thr])
	#writedlm("convergence_curve_mc.txt", convergence_curve_mc)
	#writedlm("KGerror_mc.txt", KGerror_mc)
	#writedlm("parameters_mc.txt", [tau])
	#writedlm("NMSE_mc.txt", NMSE_mc)

	#export_animation_obj(mc, reconstruction_path, "mc_"*string(under_sample_ratio)*"_frame_")
	if verbose
	println("Threshold KG error:"*string(KGerror_thr))
	println("KG error MC (Subsampling ratio:"*string(round(under_sample_ratio*100))*"%):"*string(round(KGerror_mc,4)))
	println("NMSE MC:"*string(NMSE_mc))
	end
	return matrix_mc, convergence_curve_mc, current_iter,error_between_iterations
end
