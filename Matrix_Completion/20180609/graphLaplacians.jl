function create_graph(Data::Array{Float64, 2}, construction_type::String, 
						  construction_parameter, weight_type::String, verbose=true)

	m,n = size(Data)
	adj_matrix=spzeros(n,n)::SparseMatrixCSC{Float64, Int64}
	lpcn_matrix = spzeros(n,n)::SparseMatrixCSC{Float64, Int64}

	norm_factor = norm(Data)

	if(construction_type=="e-n")
	
		I = speye(n,n)
		D = zeros(n)
		for i=1:n
			if(verbose)
				println(i)
			end
		
			for j=1:n
				if(i != j)
					dist = norm(Data[:, i] - Data[:, j])/norm_factor
					if(dist<construction_parameter)
						if(weight_type=="heat")
							if(verbose)
								println(dist)
							end
							adj_matrix[i,j] = exp(-dist)
						elseif(weight_type=="binary")
							adj_matrix[i,j] =  1
						end
					end
				end
			end
			tmp = sum(adj_matrix[i,:])
			if(tmp>0)
				D[i]=1/tmp
				lpcn_matrix[i,:] = I[i,:]  - D[i]*adj_matrix[i,:]
			else
				lpcn_matrix[i,:] = spzeros(1,n)
			end
		end

		
		
	elseif(construction_type=="k-nn")
	
		for i=1:n
			if(verbose)
				println(i)
			end

			dist=spzeros(n)
			for j=1:n
				dist[j] = norm(Data[:, i] - Data[:, j])
			end

			index = sortperm(dist)
			if(weight_type=="heat")
				adj_matrix[i,index[1:construction_parameter]] = exp.(-dist[index[1:construction_parameter]]/n)
			elseif(weight_type=="binary")
				adj_matrix[i,index[1:construction_parameter]] = ones(1, construction_parameter)
			end
		end

		D = spdiagm(1./sum(adj_matrix, 2)[:,1])
		lpcn_matrix = speye(n,n) - D*adj_matrix

	end

 	return lpcn_matrix 	 
end

function compareEstimatedLaplacian(anm::Animation,tau::Float64,max_admm_iters::Int64,knn::Int64)

	anm_matrix = arrange(anm.data, [[:coors, :frames], :vertices])
	m,n=size(anm_matrix)

	L_gt = create_graph(anm_matrix, "k-nn", knn, "binary", true)
	
	undersample=[0.1 0.2 0.3 0.4]
	L_under_error = zeros(length(undersample))
	L_recon_error = zeros(length(undersample))
	iter = 1
	for u = undersample
		println(iter)
		Omega=zeros(m,n)
		for f=1:3:m
			s = Int(round(u*n))
			Omega_indices=randperm(n)[1:s]
			Omega[f, Omega_indices]=ones(s)
			Omega[f+1, :] = Omega[f, :]
			Omega[f+2, :] = Omega[f, :]
		end
		UndersampledData = Omega.*anm_matrix
		
		ReconData, ~ = proposed_mc(anm_matrix, Omega, tau, max_admm_iters, true)

		L_undersampled = create_graph(mean(UndersampledData, 1), "k-nn", knn, "binary", true)
		L_recon = create_graph(mean(ReconData, 1), "k-nn", knn, "binary", true)
		
		L_under_error[iter] = vecnorm(L_gt - L_recon, 2)/vecnorm(L_gt,2)
		L_recon_error[iter] = vecnorm(L_gt - L_undersampled, 2)/vecnorm(L_gt,2)
		iter = iter + 1
	end
	
	#using PyPlot
	#semilogy(missing, kgerror_lsm_hc, marker="o", label="LSM, spatial coherence ",linestyle="--")
	#([0.363483,0.23018,0.172262,0.137269],[0.576946,0.576698,0.576518,0.576345])
	
	return L_under_error,L_recon_error

end