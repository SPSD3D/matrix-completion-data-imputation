#= Programming =#
# result[1] = f0(seq[1]), result[n] = cf(result[n-1], seq[n])
function chain_reduce(f0, cf, seq)
    len = length(seq)
    if len == 0 return []
    else
        seed = f0(seq[1])
        result = [seed]
        for el in drop(seq, 1)
            seed = cf(seed, el)
            push!(result, seed)
        end
    end
    return result
end

# result[1] = f0(seq[1]), result[n] = cf(seq[n-1], seq[n])
function chain_map(f0, cf, seq)
    len = length(seq)
    if len == 0 return []
    else
        seed = f0(seq[1])
        result = [seed]
        for i in 2:length(seq)
            seed = cf(seq[i-1], seq[i])
            push!(result, seed)
        end
    end
    return result
end

remove_if(fn, collection) = filter(el -> !fn(el), collection)
rows{T <: Any}(A::AbstractArray{T, 2}) = [A[i, :] for i = 1:size(A, 1)]
cols{T <: Any}(A::AbstractArray{T, 2}) = [A[:, i] for i = 1:size(A, 2)]

#= Compression =#
# Returns a random orthonormal base for a c-dimensional subspace of R^r
function randn_basis(r, c, orthogonal=true)
    if orthogonal
        return qr(randn(r, c))[1]
    else
        W = randn(r, c)
        for i in 1:c
            W[:, i] = normalize(W[:, i])
        end
        return(W)
    end
end

# Returns the number of vectors (k) needed to capture the desired metric down to
# the specified orders of magnitude, according to the given singular values.
function log_capture(svalues; order = 4, metric = :variance)
    metric_expts = Dict(:variance => 1, :energy => 2)
    metric_values = svalues .^ metric_expts[metric]
    log_ratios = map(v -> log10(first(metric_values) / v), metric_values)
    return searchsortedlast(log_ratios, order)
end

# Sample a fraction of the columns of X
function sample_columns(X, fraction, sampling = :random)
    xc = size(X, 2)
    if sampling == :random
        return [X[:, ci] for ci in rand(1:xc, round(Int64, fraction*xc))]
    elseif sampling == :uniform
        return [X[:, ci] for ci in 1:round(Int64, 1/fraction):xc]
    end
end

# KG error and bpvf (temporal and spatial)
bpvf_temp(d, k, n, qf, qd) = (6/d) + ((3*k)/(n*d)) * (n*qf + d*qd)

bpvf_spat(d, k, n, qf, qd) = (6/d) + ((3*k)/(n*d)) * (d*qf + n*qd)

function kg_error(D_hat, D, dimensions = 3)
    function trajectory_means(D)
        frames = [D[dimensions*(i-1) + 1 : dimensions*i, :]
                  for i = 1:div(size(D, 1), dimensions)]
        return  repmat(mean(frames), size(frames, 1))
    end
    numerator = vecnorm(D - D_hat)
    denominator = vecnorm(D - trajectory_means(D_hat))
    error = 100 * numerator / denominator
    return error, numerator, denominator
end

# Calculation of variance (or energy) captured by an estimated k-dimensional
# est_basis relative to those captured by the corresponding svd_basis and the
# totals of the singular values (svalues). Metric can be either :variance or
# :energy.
function basis_dist(est_basis::Array{Float64, 2}, svd_basis::Array{Float64, 2},
                    svalues::Array{Float64, 1}, metric = :energy)
    metric_expts = Dict(:variance => 1, :energy => 2)
    weights = svalues .^ metric_expts[metric]
    total = sum(weights)

    projs       = abs(est_basis' * svd_basis)
    k, n        = size(projs)
    svd_capture = sum(weights[1:k])

    # filter the maximum projection for each vector in svd_basis, which equals
    # the portion of the captured value corresponding to that vector.
    max_projs           = zeros(k, n)
    max_proj_ind_tuples = [findmax(projs[:, c]) for c in 1:n]
    for (c, (max_proj, ind)) in enumerate(max_proj_ind_tuples)
        max_projs[ind, c] = max_proj
    end
    pca_capture = sum(max_projs[:, 1:k] * weights[1:k])

    return pca_capture, svd_capture, total
end
