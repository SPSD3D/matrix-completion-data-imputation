function proposed_mc(
	M::Array{Float64, 2},
	Omega::Array{Float64, 2},
	tau::Float64,
	max_admm_iters::Int64,
	verbose=true)

	# Parameter setup
	m,n=size(M)
	rho=1/10
	overrelax_parameter=1.5
	epsilon=1e-6

	# Initializations
	norm_error = zeros(max_admm_iters)
	OM = Omega.*M
	norm_factor = vecnorm(M,2)
	X = zeros(m,n)
	Y = zeros(m,n)
	Z = zeros(m,n)
	A=spdiagm(vec(Omega),0)+rho*speye(m*n,m*n)

    current_iter=0;

	# ADMM iterations
	for admm_iter=1:max_admm_iters
		current_iter=1.0*admm_iter

		## Sub-optimization problem 1 (the soft-thresholding)

		X = svt(Y - 1/rho*Z,tau/rho)

		# if( res_error[admm_iter]<= epsilon)
		# 	break
		# end

		## Sub-optimization problem 2 (the generalized Sylvester equation)

		B=sparse(vec(OM+rho*(overrelax_parameter*X+(1-overrelax_parameter)*Y) + Z))
		Y = mysylvester(A,B,m,n)

		## Sub-optimization problem 3 (the dual variable update)
		Z = Z + rho*(X-Y)


		## Print verbose information
		norm_error[admm_iter]=vecnorm(X - M, 2)/vecnorm(M, 2)

		#norm_error[admm_iter]=(mean((X-M).^2)/mean((M-mean(M)).^2))

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
