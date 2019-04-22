using PyPlot

include("lasso_admm.jl")

function make_sparse_vectors(length, sparsity_range)
    map(sparsity -> make_sparse_vector(length, sparsity),
        collect(sparsity_range))
end

norm_1(x) = sum(abs(x))
norm_2(x) = sum(x .^ 2)

norm_1_diff(xh, x0) = norm_1(xh-x0)
norm_1_frac(xh, x0) = norm_1(xh-x0) / norm_1(x0)
norm_2_diff(xh, x0) = norm_2(xh-x0)
norm_2_frac(xh, x0) = norm_2(xh-x0) / norm_2(x0)

function test_boyd(m, n, sparsity_range; lambda = 1.0, rho = 1.2, alpha = 1.4)
    A = randn(m, n)
    x0s = make_sparse_vectors(n, sparsity_range)
    ys = map(x -> A*x, x0s)
    x_z_hs = map(y -> boyd_lasso_admm(A, y, lambda, rho, alpha), ys)
    xhs = map(s -> s[1], x_z_hs)
    zhs = map(s -> s[2], x_z_hs)
    sparsities = collect(sparsity_range)
    xn2fs = map(norm_2_frac, xhs, x0s)
    zn2fs = map(norm_2_frac, zhs, x0s)
    for i in 1:length(sparsities)
        plot_diff(x0s[i],
                  xhs[i],
                  "x, sparsity = $(sparsities[i]), lambda = $lambda, rho = $rho, alpha = $alpha",
                  xn2fs[i])
    end
    for i in 1:length(sparsities)
        plot_diff(x0s[i],
                  zhs[i],
                  "z, sparsity = $(sparsities[i]), lambda = $lambda, rho = $rho, alpha = $alpha",
                  zn2fs[i])
    end
    # println("sparsity     x dist    z dist")
    # for i in 1:length(sparsities)
    #     print("  ", rpad(sparsities[i], 10, ' '))
    #     @printf("%03.6f  ", xn2fs[i])
    #     @printf("%03.6f  ", zn2fs[i])
    #     println()
    # end
    x0s, xhs, zhs
end

function plot_diff(original, estimation, title_supplement = "",
                   distance = 0.0, save_dir = "~/data/figures/admm/")
    wd = pwd()
    cd(expanduser(save_dir))
    PyPlot.close()
    n = length(original)
    scatter(collect(1:n), original, color = "green")
    scatter(collect(1:n), estimation, color = "red",
            label = string("distance = $distance"))
    xlabel("Element index"); ylabel("Value")
    legend()
    title_string = "Original (green), estimated (red)"
    if title_supplement != ""
        title_string = title_string * ", " * title_supplement
    end
    title(title_string)
    savefig("$title_supplement.svg")
    PyPlot.close()
    cd(wd)
end
