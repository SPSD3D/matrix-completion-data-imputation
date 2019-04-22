include("utilities.jl")

# Application of vr_pca on matrix X for leading eigenvector estimation, using a
# fraction of X, with learning rate eta and initial estimation W
function vr_pca(X, iterations, fraction, eta; W = 0, sampling = :random)
    xr, xc = size(X)
    if (W==0) W = randn(xr, 1) end
    Wt = W
    # main iteration loop
    for s in 1:iterations
        U  = (1/xc) * X * X' * W
        samples = sample_columns(X, fraction, sampling)
        println("Using $(length(samples)) columns from input matrix")
        for (m, x) in enumerate(samples)
            Wt += eta * (x * (x'*Wt - x'*W) + U)
            Wt /= vecnorm(Wt)
        end
        W = Wt
    end
    return W
end

# Application of vr_pca on matrix X for k-dimensional subspace estimation, using
# a fraction of X, with learning rate eta and initial estimation W
function vr_pca_block(X, k, iterations, fraction, eta; W = 0, sampling = :random)
    # if no estimated basis W is supplied, a random k-dimensional one is used,
    # otherwise it's cropped to xr-by-min(k, xr) dimensions in order to be
    # multiplied with the covariance matrix
    xr, xc = size(X)
    if (W==0) W = randn_basis(xr, k)
    else W = W[1:xr, 1:min(k, xr)]
    end
    Wt = W
    # main iteration loop
    for s in 1:iterations
        U  = (1/xc) * X * X' * W
        samples = sample_columns(X, fraction, sampling)
        println("Using $(length(samples)) columns from input matrix.") 
        for (m, x) in enumerate(samples)
            BU, BS, BV = svd(Wt' * W)
            B = BV * BU'
            Wt += eta * (x * (x'*Wt - x'*W*B) + U*B)
            Wt *= inv(chol(Wt'*Wt))
        end
        W = Wt
    end
    return W
end
