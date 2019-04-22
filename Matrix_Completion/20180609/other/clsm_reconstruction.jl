# Animation library
include("../lib/animation.jl")

include("lsm.jl")

# The solution of subproblem 2 of ADMM is obtained via 
# a call to the conjugate gradient algorithm.
include("../lib/cg.jl")

using JLD

### Reconstruct the input point cloud based on the proposed graph-based matrix completion algorithm
function clsm_reconstruction(reconstruction_path::String, anm::Animation, anm_matrix_rank::Int64,
						 under_sample_ratio::Float64, cluster_size::Int, graph_type::String)
			 

	# Parameters setup
	anm_matrix = arrange(anm.data, [[:coors, :frames], :vertices])
	m,n = size(anm_matrix)

	# Ground-trouth
	U,S,V = svd(anm_matrix)
	thr_S = zeros(size(S))
	thr_S[1:anm_matrix_rank] = S[1:anm_matrix_rank]
	anm_matrix_thr = U*diagm(thr_S)*V'


	# anm_matrix_rank2 = round(Int, (anm_matrix_rank/m)*cluster_size)
	# anm_matrix_thr = zeros(m,n)
	# for ki=1:cluster_size:m
		# block_i = ki:ki+cluster_size-1
		# data_per_block = anm_matrix[block_i,:]
		# U,S,V = svd(data_per_block)
		# thr_S = zeros(size(S))
		# thr_S[1:anm_matrix_rank2] = S[1:anm_matrix_rank2]
		# anm_matrix_thr[block_i,:] = U*diagm(thr_S)*V'
	# end
	
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
	
	Omega_cl=zeros(m,n)
	epsilon = 0.001
	for k=cluster_size+1:cluster_size:m-cluster_size+1
		proximity_metric_vec = anm_matrix[k-3:k-1, :] - anm_matrix[k:k+2, :]
		proximity_metric_ind = find(sqrt(sum(abs(proximity_metric_vec).^2,1)).< epsilon)
		Omega_cl[k-3:k+2, proximity_metric_ind] = ones(6, length(proximity_metric_ind))
	end

	if(graph_type=="average")
		L = sparse(readdlm(reconstruction_path*string("average_graph.txt")))
	elseif(graph_type=="block_average")
		L = spzeros(round(Int,m/3)*n,n)
		for f=1:round(Int,m/3)-1
			data = load(reconstruction_path*string("L_"*string(f)*".jld"))
			data_label=collect(data)[1][1]
			L[n*(f-1)+1:n*f,:] = data[data_label]
		end
	else
		L = speye(n,n)
	end
		
	anm_matrix_clsm = zeros(m,n)
	for ki=1:cluster_size:m
		block_i = ki:ki+cluster_size-1
		f = round(Int,ki/3)+1		
		anm_matrix_clsm[block_i,:]  = lsm(anm_matrix[block_i,:], Omega[block_i,:], L[n*(f-1)+1:n*f,:])
	end

	# Performance evaluations
	NMSE_clsm=vecnorm(anm_matrix_clsm - anm_matrix, 2)/vecnorm(anm_matrix,2) 	
	KGerror_clsm=100*vecnorm(anm_matrix-anm_matrix_clsm, 2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	KGerror_thr=100*vecnorm(anm_matrix_thr-anm_matrix,2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	
	# I/0
	anm_clsm = deepcopy(anm)
	anm_clsm.data = LArray(reshape(anm_matrix_clsm, 3, round(Int, m/3), n), [:coors, :frames, :vertices])
	
	println("Threshold KG error:"*string(KGerror_thr))
	writedlm(reconstruction_path*string(under_sample_ratio)*"_KGerror_clsm.txt", KGerror_clsm)
	writedlm(reconstruction_path*string(under_sample_ratio)*"_NMSE_clsm.txt", NMSE_clsm)
	writedlm(reconstruction_path*string(under_sample_ratio)*"_clsm_parameters.txt", [cluster_size])	
	export_animation_obj(anm_clsm, reconstruction_path, "clsm_"*string(under_sample_ratio)*"_frame_")
	println("KG error CLSM (Subsampling ratio:"*string(round(under_sample_ratio*100))*"%):"*string(round(KGerror_clsm,4)))
	println("NMSE CLSM:"*string(NMSE_clsm))
end