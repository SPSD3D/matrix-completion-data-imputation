
# Animation library
include("../lib/animation.jl")

# Matrix completion technique
include("proposed.jl")
include("graphLaplacians.jl")

# The solution of subproblem 2 of ADMM is obtained via
# a call to the conjugate gradient algorithm.
include("../lib/cg.jl")

function proposed_mc_candes(
	M::Array{Float64, 2},
	Omega::Array{Float64, 2},
	tau::Float64,
	max_admm_iters::Int64,
	verbose=true)
	epsilon=1e-6
	rho=1/100
	# Initializations
	norm_error = zeros(max_admm_iters)
	m,n=size(M)
	X = zeros(m,n)
	Y = zeros(m,n)
    current_iter=0;

	# ADMM iterations
	for admm_iter=1:max_admm_iters
		current_iter=1.0*admm_iter


		X = svt(Y,tau/rho)
		Y=Y+0.2*Omega.*(M-X)


		## Print verbose information
		norm_error[admm_iter]=vecnorm(X - M, 2)/vecnorm(M,2)
		if(admm_iter>1)
		if( abs(norm_error[admm_iter]-norm_error[admm_iter-1])<= epsilon)
			break
		end
		end

		if(verbose)
			println("Iteration: "*string(admm_iter)*":"*string(norm_error[admm_iter]))
		end
	end

	return X, norm_error,current_iter
end
