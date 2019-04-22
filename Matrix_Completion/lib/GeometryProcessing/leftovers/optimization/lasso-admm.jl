using JLD
using CPUTime
using ParallelDataTransfer

function make_sparse_vector(length::Int, sparsity::Number; t=Float)
    # if given in percent, adjust sparsity to fraction
    if sparsity > 1 sparsity /= 100 end
    entry_count = round(Int, sparsity * length)
    positions = zeros(Int, entry_count)
    for entry in 1:entry_count
        entry_pos = rand(1:length)
        # ensure unique positions
        while entry_pos in positions
            entry_pos = rand(1:length)
        end
        positions[entry] = entry_pos
    end
    sparse_vector = zeros(t, length)
    for pos in positions
        sparse_vector[pos] = randn(t) end
    sparse_vector
end

# Accept a vector y and a fat matrix A as input, and solve the minimization
# problem argmin_x(1/2*||y - Ax||_2^2 + lambda||x||_1) by performing LASSO using
# ADMM, returning the minimizing vector x and the vector of costs at each
# iteration.
function lucey_lasso_admm(A::Matrix{Float}, y::Vector{Float},
                          lambda::Float)
    Ac = size(A, 2)
    # parameter initialization
    max_iter = 500; rho = 1e-4; max_rho = 5;
    I = speye(Ac); lngr = zeros(Ac); c = randn(Ac)

    # fast soft thresholding and norm functions
    fast_sthresh(x, th) = sign(x) .* max(abs(x) - th, 0)
    norm_1(x) = sum(abs(x[:]))
    norm_2(x) = sum(x .^ 2)

    x = Vector{Float}(Ac)
    cost = Vector{Float}(max_iter)
    for n = 1:max_iter
        # solve sub-problem to get x
        x = (A'*A + rho*I) \ (A'*y + rho*c - lngr)
        # solve sub-problem to get c
        c = fast_sthresh(x + lngr/rho, lambda/rho)
        # update the lagrangian
        lngr = lngr + rho*(x - c)
        # section 3.3 in Boyd's book
        rho = min(max_rho, rho * 1.1)
        # get the current cost
        cost[n] = 0.5 * norm_2(y - A*x) + lambda * norm_1(x)
    end
    x
end

# Solve lasso problem via ADMM (adapted from http://www.stanford.edu/~boyd):
# argmin_x 1/2*||Ax - b||_2^2 + lambda*||x||_1 . Rho is the augmented Lagrangian
# parameter, alpha is the over-relaxation parameter (typical values for alpha
# are between 1.0 and 1.8).
function boyd_lasso_admm(A::Matrix{Float}, y::Vector{Float};
                         lambda = 1.0, rho = 1.0, alpha = 1.2,
                         eps = 1e-6)
    # objective, shrinkage and cholesky (+ rho*I) factorization
    # objective(A, y, lambda, x, z) = 0.5 * sum((A*x - y).^2) + lambda*norm(z, 1)
    shrinkage(x, kappa) = max(0, x-kappa) - max(0, -x-kappa)
    chol_rho(A, rho) = chol(eye(size(A, 1)) + 1/rho*(A*A'))

    max_iter = 1000; Ac = size(A, 2); Aty = A'*y
    x = zeros(Ac); z = zeros(Ac); u = zeros(Ac); x_old = zeros(Ac); error = Inf;
    U = chol_rho(A, rho); L = U'

    for k = 1:max_iter
        q     = Aty + rho*(z - u)
        x_old = x
        x     = q/rho - (A'*(U \ ( L \ (A*q)))) / rho^2
        x_hat = alpha*x + (1 - alpha)*z;
        z     = shrinkage(x_hat + u, lambda/rho);
        u     = u + (x_hat - z);
        error = min(error, vecnorm(x - x_old) / vecnorm(x_old))
        error < eps && break
    end
    x, z, error
end

function remote_eval(pid, form)
    remotecall_fetch(eval, pid, form)
end

macro reval(pids, exp)
    quote
        @sync for pid in $pids
            @async remotecall_fetch(eval, pid, $exp)
        end
    end
end

function block_limits(length, block_size)
    i = 1; t = block_size; bs = [];
    while i <= length
        push!(bs, i:min(t, length))
        i += block_size
        t += block_size
    end
    return bs
end

function save_blocks(filename::String, A::Array, block_count::Int)
    m = size(A, 1)
    blims = block_limits(m, ceil(Int, m / block_count))
    blocks = Dict{String, Matrix}()
    for (bid, blim) in enumerate(blims)
        blocks[string(bid)] = A[blim, :]
    end
    save(filename, blocks)
end

function distributed_load(filename::String, var::Symbol, workers::Vector{Int})
    for (pid_id, pid) in enumerate(workers)
        remote_eval(pid, :($var = load($filename, $(string(pid_id)))))
    end
end

function admm_parallel(A::Matrix{Float}, y::Vector{Float}; N = 0,
                       lambda = 1.0, rho = 1.0, abstol = 1e-4, reltol = 1e-2,
                       distribute::Bool = false)

    distribute && nprocs() == 1 && addprocs(N != 0 ? N : Sys.CPU_CORES)
    worker_pids = workers(); N = length(worker_pids); n = size(A, 2)
    blims = block_limits(size(A, 1), ceil(Int, size(A, 1) / N))

    # send data and parameters to worker processes
    @sync for (bid, pid) in enumerate(worker_pids)
        @async begin
            remote_eval(pid, :(A = $(A[blims[bid], :])))
            remote_eval(pid, :(y = $(y[blims[bid]])))
            remote_eval(pid, :(lambda = $lambda; rho = $rho; n = $n))
        end
    end

    # define cholesky decomposition function
    @everywhere function chol_rho(A)
        chol(eye(size(A, 1)) + 1/rho*(A*A'))
    end

    @sync for pid in worker_pids
        @async begin
            remote_eval(pid, :(Aty = A'*y))
            remote_eval(pid, :(m = size(A, 1)))
            remote_eval(pid, :(U = chol_rho(A)))
            remote_eval(pid, :(L = U'))
            remote_eval(pid, :(u = zeros(n)))
            remote_eval(pid, :(x = zeros(n)))
            remote_eval(pid, :(r = zeros(n)))
        end
    end

    # admm_iteration returns w (= x + u) after one ADMM iteration at each
    # chunk. Relying on presence of predefined and cached variables (Aty, U, u,
    # x).
    @everywhere function admm_iteration(z)
        r = x - z
        global u += x - z
        q = Aty + rho*(z - u)
        global x = q/rho - (A'*(U \ ( L \ (A*q)))) / rho^2
        w = x + u
        # return  w, r^2, x^2 and u^2/rho^2
        return x + u, dot(r, r), dot(x, x), dot(u, u) / rho^2
    end

    function objective(A, y, x)
        0.5 * sum((A*x - y).^2) + lambda*norm(x, 1)
    end

    shrinkage(x) = max(0, x - lambda/(N*rho)) - max(0, -x - lambda/(N*rho))
    remote_admm_iteration(pid, z) = remotecall_fetch(admm_iteration, pid, z)
    function aggregate(outputs)
        z_mean = zeros(n); prires = 0; nxstack = 0; nystack = 0;
        for (z, r2, x2, u2rho2) in outputs
            z_mean += z; prires += r2; nxstack += x2; nystack += u2rho2;
        end
        return shrinkage(z_mean/N), sqrt(prires), sqrt(nxstack), sqrt(nystack)
    end

    max_iterations = 1000
    zprev = zeros(n); z = zeros(n);
    iteration_outputs = Vector{Any}(N)

    @time @CPUtime @profile for iter = 1:max_iterations
        # execute iterations remotely with available z
        @sync for (id, pid) in enumerate(worker_pids)
            @async iteration_outputs[id] = remote_admm_iteration(pid, z)
        end

        # aggregate the results and compute residuals
        z, prires, nxstack, nystack = aggregate(iteration_outputs)
        dualres = sqrt(N) * rho * norm(z - zprev, 2)
        eps_pri  = sqrt(n*N)*abstol + reltol * max(nxstack, sqrt(N)*norm(z, 2))
        eps_dual = sqrt(n*N)*abstol + reltol * nystack

        # check for termination conditions
        prires <= eps_pri && dualres <= eps_dual && break

        # log output, cache z
        println("Iteration $iter, residuals: primary = $prires, dual = $dualres")
        zprev = z
    end
    rmprocs(worker_pids)
    return z
end

function async_admm_parallel(A::Matrix{Float}, y::Vector{Float}; N = 0,
                             lambda = 1.0, rho = 1.0, abstol = 1e-4,
                             reltol = 1e-2, distribute::Bool = false)

    distribute && nprocs() == 1 && addprocs(N != 0 ? N : Sys.CPU_CORES)
    worker_pids = workers(); N = length(worker_pids); n = size(A, 2)
    blims = block_limits(size(A, 1), ceil(Int, size(A, 1) / N))

    # send data and parameters to worker processes
    @sync for (bid, pid) in enumerate(worker_pids)
        @async begin
            remote_eval(pid, :(A = $(A[blims[bid], :])))
            remote_eval(pid, :(y = $(y[blims[bid]])))
            remote_eval(pid, :(lambda = $lambda; rho = $rho; n = $n))
        end
    end

    # define cholesky decomposition function
    @everywhere function chol_rho(A)
         chol(eye(size(A, 1)) + 1/rho*(A*A'))
    end

    @sync for pid in worker_pids
        @async begin
            remote_eval(pid, :(Aty = A'*y))
            remote_eval(pid, :(m = size(A, 1)))
            remote_eval(pid, :(U = chol_rho(A)))
            remote_eval(pid, :(L = U'))
            remote_eval(pid, :(u = zeros(n)))
            remote_eval(pid, :(x = zeros(n)))
            remote_eval(pid, :(r = zeros(n)))
        end
    end

    # admm_iteration returns w (= x + u) after one ADMM iteration at each
    # chunk. Relying on presence of predefined and cached variables (Aty, U, u,
    # x).
    @everywhere function admm_iteration(z)
        r = x - z
        global u += x - z
        q = Aty + rho*(z - u)
        global x = q/rho - (A'*(U \ ( L \ (A*q)))) / rho^2
        w = x + u
        # return  w, r^2, x^2 and u^2/rho^2
        return x + u, dot(r, r), dot(x, x), dot(u, u) / rho^2
    end

    function objective(A, y, x)
        0.5 * sum((A*x - y).^2) + lambda*norm(x, 1)
    end

    shrinkage(x) = max(0, x - lambda/(N*rho)) - max(0, -x - lambda/(N*rho))
    remote_admm_iteration(pid, z) = remotecall_fetch(admm_iteration, pid, z)
    function aggregate(outputs)
        z_mean = zeros(n); prires = 0; nxstack = 0; nystack = 0;
        for (z, r2, x2, u2rho2) in outputs
            z_mean += z; prires += r2; nxstack += x2; nystack += u2rho2;
        end
        return shrinkage(z_mean/N), sqrt(prires), sqrt(nxstack), sqrt(nystack)
    end

    max_iterations = 1000
    zprev = zeros(n); z = zeros(n);
    iteration_outputs = Vector{Any}(N)

    @time @CPUtime @profile for iter = 1:max_iterations
        # execute iterations remotely with available z
        @sync for (id, pid) in enumerate(worker_pids)
            @async iteration_outputs[id] = remote_admm_iteration(pid, z)
        end

        # aggregate the results and compute residuals
        z, prires, nxstack, nystack = aggregate(iteration_outputs)
        dualres = sqrt(N) * rho * norm(z - zprev, 2)
        eps_pri  = sqrt(n*N)*abstol + reltol * max(nxstack, sqrt(N)*norm(z, 2))
        eps_dual = sqrt(n*N)*abstol + reltol * nystack

        # check for termination conditions
        prires <= eps_pri && dualres <= eps_dual && break

        # log output, cache z
        println("Iteration $iter, residuals: primary = $prires, dual = $dualres")
        zprev = z
    end
    rmprocs(workers)
    return z
end
