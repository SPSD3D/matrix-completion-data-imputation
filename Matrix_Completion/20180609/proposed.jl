# Proposed algorithms for graph-based rank minimization



function proposed_gmc_v2(
	M::Array{Float64, 2},
	Omega::Array{Float64, 2},
	tau::Float64,
	max_admm_iters::Int64,
	gamma::Float64,
	L1::SparseMatrixCSC{Float64,Int64},
	#L2::SparseMatrixCSC{Float64,Int64},
				 verbose=true)

	# Parameter setup
	m,n=size(M)
	rho=1/100
	#rho=1/10
	overrelax_parameter=1.5
	epsilon=1e-8

	# Initializations
	norm_error = zeros(max_admm_iters)
	error_between_iterations = zeros(max_admm_iters-1)
	OM = Omega.*M
	norm_factor = vecnorm(M,2)
	X = zeros(m,n)
	Y = zeros(m,n)
	Z = zeros(m,n)


	#mlen=size(M,1);
    #mlen=size(M,direction);



	#A=spdiagm(vec(Omega),0)+rho*speye(m*n,m*n)+gamma*kron(L'*L, speye(mlen,mlen))

	current_iter=0;

	# ADMM iterations
	for admm_iter=1:max_admm_iters
        current_iter=1.0*admm_iter
		## Sub-optimization problem 1 (the soft-thresholding)
		X = svt(Y - 1/rho*Z,tau/rho)






		C=gamma*L1'*L1
		#C=gamma*L1'*L1+gamma*L2'*L2

		A=spdiagm(vec(Omega),0)+rho*speye(m*n,m*n)+C
		B=sparse(vec(OM+rho*(overrelax_parameter*X+(1-overrelax_parameter)*Y) + Z))
		Y = mysylvester(A,B,m,n)

		## Sub-optimization problem 3 (the dual variable update)
		Z = Z + rho*(X-Y)

		## Print verbose information
		#norm_error[admm_iter]=vecnorm(X - M, 2)/norm_factor
		norm_error[admm_iter]=(mean((X-M).^2)/mean((M-mean(M)).^2))

		if(admm_iter>1)

			error_between_iterations[admm_iter-1]=abs(norm_error[admm_iter]-norm_error[admm_iter-1])


		if( error_between_iterations[admm_iter-1]<= epsilon)
			break
		end


		end

		if(verbose)
			println("Iteration: "*string(admm_iter)*":"*string(norm_error[admm_iter]))
		end
	end
	#println(current_iter)
	return X, norm_error,current_iter,error_between_iterations
end











function proposed_gmc(M::Array{Float64, 2}, Omega::Array{Float64, 2}, tau::Float64, max_admm_iters::Int64,
				 gamma::Float64, L::SparseMatrixCSC{Float64,Int64},direction::Int64,
				 verbose=true)

	# Parameter setup
	m,n=size(M)
	rho=1/100
	#rho=1/10
	overrelax_parameter=1.5
	#epsilon=1e-6

	# Initializations
	norm_error = zeros(max_admm_iters)
	OM = Omega.*M
	norm_factor = vecnorm(M,2)
	X = zeros(m,n)
	Y = zeros(m,n)
	Z = zeros(m,n)


	#mlen=size(M,1);
    mlen=size(M,direction);



	#A=spdiagm(vec(Omega),0)+rho*speye(m*n,m*n)+gamma*kron(L'*L, speye(mlen,mlen))

	C=gamma*(kron(L'*L, speye(mlen,mlen)))

    A=spdiagm(vec(Omega),0)+rho*speye(m*n,m*n)+C

	# ADMM iterations
	for admm_iter=1:max_admm_iters

		## Sub-optimization problem 1 (the soft-thresholding)
		X = svt(Y - 1/rho*Z,tau/rho)



		## Sub-optimization problem 2 (the generalized Sylvester equation)
		B=sparse(vec(OM+rho*(overrelax_parameter*X+(1-overrelax_parameter)*Y) + Z))
		Y = mysylvester(A,B,m,n)

		## Sub-optimization problem 3 (the dual variable update)
		Z = Z + rho*(X-Y)

		## Print verbose information
		norm_error[admm_iter]=vecnorm(X - M, 2)/norm_factor
		if(verbose)
			println("Iteration: "*string(admm_iter)*":"*string(norm_error[admm_iter]))
		end
	end

	return X, norm_error
end





function proposed_cgmc(M::Array{Float64, 2}, Omega::Array{Float64, 2}, tau::Float64, max_admm_iters::Int64,
				 gamma::Float64, L::SparseMatrixCSC{Float64, Int64}, cluster_size::Int64, Omega_cl::Array{Float64, 2},
				 xu::Float64, verbose=true)

	# Parameter setup
	m,n=size(M)
	rho=1/10000
	overrelax_parameter=1.5
	#epsilon=1e-6

	# Initializations
	norm_error = zeros(max_admm_iters)
	OM = Omega.*M
	norm_factor = vecnorm(M,2)
	X = zeros(m,n)
	Y = zeros(m,n)
	Z = zeros(m,n)
	V = zeros(m,n)
	Zcl = zeros(m,n)

	KLL = spzeros(round(Int,m/3)*n*cluster_size,n*cluster_size)
	for ki=1:cluster_size:m
		f = round(Int,ki/3)+1
		KLL[n*cluster_size*(f-1)+1:n*f*cluster_size, :] = kron(L[n*(f-1)+1:n*f,:]'*L[n*(f-1)+1:n*f,:], speye(cluster_size,cluster_size))
	end


	# ADMM iterations
	for admm_iter=1:max_admm_iters

		## Sub-optimization problem 1 (the soft-thresholding)
		for ki=1:cluster_size:m
			block_i = ki:ki+cluster_size-1
			X[block_i,:] = svt(Y[block_i,:] - 1/rho * V[block_i,:], tau/rho)
		end

		#if( res_error[admm_iter]<= epsilon)
		#	break
		#end

		## Sub-optimization problem 2 (the generalized Sylvester equation)
		for ki=1:cluster_size:m
			block_i = ki:ki+cluster_size-1
			f = round(Int,ki/3)+1
			Y[block_i,:] = mysylvester(spdiagm(vec(Omega[block_i,:]+xu*Omega_cl[block_i,:]),0) + gamma* KLL[n*cluster_size*(f-1)+1:n*f*cluster_size, :]+
			rho*speye(cluster_size*n,cluster_size*n),sparse(vec(xu*Zcl[block_i,:] + V[block_i,:] +
			rho*(overrelax_parameter*X[block_i,:] + (1-overrelax_parameter)*Y[block_i,:]) + OM[block_i,:])),cluster_size,n)
		end

		## Sub-optimization problem 3 (Adjacent vertices averaging)
		for ki=cluster_size+1:cluster_size:m-cluster_size+1
			Zcl[ki-3:ki-1,:] = (Omega_cl[ki-3:ki-1,:].*Y[ki-3:ki-1,:]+Omega_cl[ki:ki+2,:].*Y[ki:ki+2,:])/2
			Zcl[ki:ki+2,:] = Zcl[ki-3:ki-1,:]
			#Zcl[ki-3:ki-1,:] = Omega_cl[ki-3:ki-1,:].*M[ki-3:ki-1,:]
			#Zcl[ki:ki+2,:] = Omega_cl[ki:ki+2,:].*M[ki:ki+2,:]
		end

		## Sub-optimization problem 3 (the dual variable update)
		# Sub-optimization 4
		for ki=1:cluster_size:m
			block_i = ki:ki+cluster_size-1
			V[block_i,:] = V[block_i,:] + rho*(X[block_i,:]-Y[block_i,:])
		end


		## Print verbose information
		norm_error[admm_iter]=vecnorm(X - M, 2)/norm_factor
		if(verbose)
			println("Iteration: "*string(admm_iter)*":"*string(norm_error[admm_iter]))
		end
	end

	return X, norm_error
end

# Lib
function svt(Y,tau)

	#print(tau)
	U,S,V=svd(Y)
	soft_thrs_S = max.(0, S.-tau) - max.(0, -S.-tau)

	# Reconstruct the matrix
	X = U*diagm(soft_thrs_S)*V'

	return X
end

function mysylvester(A,b,m,n)

	return reshape(cg(A,b),m,n)
end
