# X: data matrix (of size m x n) to be decomposed.
# lambda: regularization parameter, default = 1 / sqrt(max(m, n)).
# mu: augmented lagrangian parameter, default = 10 * lambda.
# tolerance - reconstruction error tolerance, default = 1e-4.
# max_iterations - maximum number of iterations, default = 1000.
function rpca(X;
              L = zeros(X),
              S = zeros(X),
              Y = zeros(X),
              lambda = 1/sqrt(max(size(X)...)),
              mu = 10 * lambda,
              tolerance = 1e-14,
              max_iterations = 4000,
              verbose = false,
              return_full_state = false)



              X_norm = vecnorm(X)

    for iter = (1:max_iterations)
        # ADMM step: update L and S
        L = Do(1/mu, X - S + (1/mu)*Y)
        S = So(lambda/mu, X - L + (1/mu)*Y)
        # and augmented lagrangian multiplier
        Z = X - L - S
        Y = Y + mu*Z

        err = vecnorm(Z) / X_norm
        if (iter == 1 || mod(iter, 10) == 0 || err < tolerance) && verbose
            println("iter: $iter err: $err rank(L): $(rank(L))")
        end
        if err < tolerance break; end
    end
    if return_full_state
        return L, S, Y, lambda, mu
    else
        return L, S
    end
end

# X: data matrix (of size m x n) to be decomposed.
# lambda: regularization parameter, default = 1 / sqrt(max(m, n)).
# mu: augmented lagrangian parameter, default = 10 * lambda.
# tolerance - reconstruction error tolerance, default = 1e-4.
# max_iterations - maximum number of iterations, default = 1000.
function rpca_no_Y(X;
                   lambda = 1/sqrt(max(size(X)...)),
                   mu = 10 * lambda,
                   tolerance = 1e-6,
                   max_iterations = 1000,
                   verbose = false,
                   return_full_state = false)
    m, n = size(X)
    X_norm = vecnorm(X)

    # initial solution
    L = zeros(m, n); S = zeros(m, n)

    for iter = (1:max_iterations)
        # ADMM step: update L and S
        L = Do(1/mu, X - S)
        S = So(lambda/mu, X - L)
        Z = X - L - S

        err = vecnorm(Z) / X_norm
        if (iter == 1 || mod(iter, 10) == 0 || err < tolerance) && verbose
            println("iter: $iter err: $err rank(L): $(rank(L))")
        end
        if err < tolerance break; end
    end
    if return_full_state
        return L, S, lambda, mu
    else
        return L, S
    end
end

# shrinkage operator
function So(tau, X)
    sign.(X) .* max.(abs.(X) - tau, 0);
end

# shrinkage operator for singular values
function Do(tau, X)
    U, S, V = svd(X)
    U * So(tau, diagm(S)) * V'
end

# X: data matrix (of size m x n) to be decomposed.
# lambda: regularization parameter, default = 1 / sqrt(max(m, n)).
# mu: augmented lagrangian parameter, default = 10 * lambda.
# tolerance - reconstruction error tolerance, default = 1e-4.
# max_iterations - maximum number of iterations, default = 1000.
function rpca_fixed_iterations(X;
                               L = zeros(X),
                               S = zeros(X),
                               Y = zeros(X),
                               lambda = 1/sqrt(max(size(X)...)),
                               mu = 10 * lambda,
                               tolerance = 1e-3,
                               iterations = 1000,
                               verbose = false,
                               return_full_state = false)
    X_norm = vecnorm(X)

    for iter = 1:iterations
        # ADMM step: update L and S
        L = Do(1/mu, X - S + (1/mu)*Y)
        S = So(lambda/mu, X - L + (1/mu)*Y)
        # and augmented lagrangian multiplier
        Z = X - L - S
        Y = Y + mu*Z

        err = vecnorm(Z) / X_norm
        if (iter == 1 || mod(iter, 10) == 0 || err < tolerance) && verbose
            println("iter: $iter err: $err rank(L): $(rank(L))")
        end
        # if err < tolerance break; end
    end
    if return_full_state
        return L, S, Y, lambda, mu
    else
        return L, S
    end
end
