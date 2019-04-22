
# Animation library
include("../lib/animation.jl")

# Matrix completion technique
include("proposed.jl")
include("graphLaplacians.jl")

# The solution of subproblem 2 of ADMM is obtained via
# a call to the conjugate gradient algorithm.
include("../lib/cg.jl")






### Reconstruct the input point cloud based on the proposed graph-based matrix completion algorithm
function mc_reconstruction_candes(
	matrix::Matrix{Float64},
	Omega::Matrix{Float64},
	matrix_rank::Int64,
	under_sample_ratio::Float64,
	tau::Float64,
	max_admm_iters::Int64,
	verbose::Bool
	)
	# Initialization
	convergence_curve_mc = zeros(max_admm_iters)



	matrix_mc, convergence_curve_mc, current_iter = proposed_mc_candes(matrix, Omega, tau, max_admm_iters, verbose)


	NMSE_mc=vecnorm(matrix_mc - matrix, 2)/vecnorm(matrix,2)







	return matrix_mc, convergence_curve_mc, current_iter
end
