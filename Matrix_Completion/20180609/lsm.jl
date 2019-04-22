function lsm(Data::Array{Float64, 2},Omega::Array{Float64, 2},Lc::SparseMatrixCSC{Float64, Int64})
	
	m,n = size(Data)
	Reconstructed = zeros(m,n)
	I=speye(n,n)
	
	wght=10
	for fi=1:m
	
		println(fi)
		#delta=Lc*Data[fi,:]
		ind = find(Omega[fi, :])
		sub=length(ind)
		Le=[Lc ; wght*I[ind, :]]
		Dframe_sub = Data[fi, ind]
		b=[spzeros(n) ; wght*Dframe_sub]
		Reconstructed[fi, :] = cg(Le'*Le,Le'*b)

	end
	
	return Reconstructed

end