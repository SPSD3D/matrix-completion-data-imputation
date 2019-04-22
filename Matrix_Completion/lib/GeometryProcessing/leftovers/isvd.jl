vecsvd(vector) = svd(reshape(vector, size(vector, 1), 1))

isvd{T<:Number}(vector::Array{T, 1}) = vecsvd(vector)

function isvd_update(U, S, V, vec)
    w = U' * vec
    p = U * w
    r = vec - p
    vr, vc = size(V)
    Uh, Sh, Vh = svd([diagm(S) w; zeros(size(S))' vecnorm(r)])
   
    return [U (normalize(r))] * Uh, Sh, [V zeros(vr); zeros(vc)' 1] * Vh
end

function isvd{T<:Number}(vectors::Array{Array{T, 1}, 1})
    U, S, V = vecsvd(first(vectors))
    for vec in drop(vectors, 1)
        U, S, V = isvd_update(U, S, V, vec)
    end
    return U, S, V
end

function isvd{T<:Number}(matrix::Array{T, 2})
    U, S, V = vecsvd(matrix[:, 1])
    for ci in drop(indices(matrix, 2), 1)
        U, S, V = isvd_update(U, S, V, matrix[:, ci])
    end
    return U, S, V
end
