# Animation library
include("../lib/animation.jl")


### Reconstruct the input point cloud based on the proposed graph-based matrix completion algorithm
function thr_reconstruction(reconstruction_path::String, anm::Animation, anm_matrix_rank::Int64)
			 

	# Parameters setup
	anm_matrix = arrange(anm.data, [[:coors, :frames], :vertices])
	m,n = size(anm_matrix)

	# Ground-trouth
	U,S,V = svd(anm_matrix)
	thr_S = zeros(size(S))
	thr_S[1:anm_matrix_rank] = S[1:anm_matrix_rank]
	anm_matrix_thr = U*diagm(thr_S)*V'

	# Performance evaluations
	KGerror_thr=100*vecnorm(anm_matrix_thr-anm_matrix,2)/vecnorm(anm_matrix.-mean(anm_matrix, 2), 2)
	
	# I/0
	anm_thr = deepcopy(anm) 
	anm_thr.data = LArray(reshape(anm_matrix_thr, 3, round(Int, m/3), n), [:coors, :frames, :vertices])

	writedlm(reconstruction_path*"thr_parameters.txt", [anm_matrix_rank])		
	export_animation_obj(anm_thr, reconstruction_path, "anm_thr_frame_")

	println("Threshold KG error:"*string(KGerror_thr))
	
	return KGerror_thr
end