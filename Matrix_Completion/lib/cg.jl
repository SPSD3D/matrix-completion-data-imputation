function cg(A, b)
    K = length(b)
    x = zeros(K)
    # Initalization for the CG algorith for each input sample
    r = b

    p = r

    rsold = r'*r
    for iter = 1:K
        #println(iter)
        # Compute the gradient vector
        Ap = A * p
        alpha =  rsold / (p' * Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = r' * r

        sqrt(rsnew) < 1e-12 && break
        p = r + (rsnew / rsold) * p

        #println(rsold)
        #println(rsnew)
        #rsold[1,1] = rsnew[1,1]

        rsold = rsnew
        
    end
    return x
end
