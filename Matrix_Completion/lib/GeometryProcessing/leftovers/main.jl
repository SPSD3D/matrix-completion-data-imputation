include("animation.jl")
include("io.jl")
include("vr_pca.jl")
include("compression.jl")

function vertex_distances(D, dimensions = 3)
    dr, dc = size(D)
    distances = []
    for frm_i in 1:div(dr, dimensions)
        frm_init = (frm_i - 1) * dimensions + 1
        frm_term = frm_init + dimensions - 1
        push!(distances, mapslices(vecnorm, D[frm_init:frm_term, :], 1))
    end
    return reduce(vcat, distances)
end

# VR_PCA
function compress_vr_pca(orig, output_name, dir = "~/data/exported";
                         k = 500, x_qbits = 12, d_qbits = 12, vr_pca_iters = 1,
                         vr_pca_fraction = 0.002, vr_pca_eta = 1000,
                         vr_pca_sampling = :uniform)
    comp = deepcopy(orig)
    comp.data, time =
        @timed batch_enc(orig.data,
                         method = :vr_pca,
                         k = k,
                         x_qbits = x_qbits,
                         d_qbits = d_qbits,
                         vr_pca_iters = vr_pca_iters,
                         vr_pca_fraction = vr_pca_fraction,
                         vr_pca_eta = vr_pca_eta,
                         vr_pca_sampling = vr_pca_sampling)
    add_data!(comp, :vertex, "distance",
              vertex_distances(comp.data - orig.data), 1)
    println("VR_PCA KG error = $(kg_error(comp.data, orig.data))")
    export_animation_pvd(comp, dir, output_name)
end

function test_vr_pca(anm, k_range, x_qbits_range, d_qbits_range,
                     vr_pca_iters_range, vr_pca_fraction_range, vr_pca_eta_range,
                     output_filename, output_directory = "~/data/logs")
    original = anm.data

    open(joinpath(ensure_dir(output_directory), output_filename), "w") do csv
        write(csv,
              "K, BLOCK_QBITS, DICT_QBITS, BPVF, ITERS, FRACTION, ",
              "ETA, KG, TIME\n")
        for
            k in collect(k_range),
            x_qbits in collect(x_qbits_range),
            d_qbits in collect(d_qbits_range),
            vr_pca_iters in collect(vr_pca_iters_range),
            vr_pca_fraction in collect(vr_pca_fraction_range),
            vr_pca_eta in collect(vr_pca_eta_range)
            encoded, time =
                @timed batch_enc(original,
                                 method = :vr_pca,
                                 k = k,
                                 x_qbits = x_qbits,
                                 d_qbits = d_qbits,
                                 vr_pca_iters = vr_pca_iters,
                                 vr_pca_fraction = vr_pca_fraction,
                                 vr_pca_eta = vr_pca_eta,
                                 vr_pca_sampling = :uniform)
            # error computation
            vr_pca_kg = kg_error(encoded, original)
            bpvf = bpvf_temp(size(anm, :frames), k,anm.vertex_count,
                             x_qbits, d_qbits)
            write(csv,
                  "$k, $x_qbits, $d_qbits, $bpvf, $vr_pca_iters, ",
                  "$vr_pca_fraction, $vr_pca_eta, $vr_pca_kg, $time\n")
        end
    end
end

# SVD
function compress_svd(orig, output_name, dir = "~/data/exported";
                      k = 500, x_qbits = 12, d_qbits = 12)
    comp = deepcopy(orig)
    comp.data, time =
        @timed batch_enc(orig.data,
                         method = :svd,
                         k = k,
                         x_qbits = x_qbits,
                         d_qbits = d_qbits)
    add_data!(comp, :vertex, "distance",
              vertex_distances(comp.data - orig.data), 1)
    println("SVD KG error = $(kg_error(comp.data, orig.data))")
    println("SVD time = $time")
    export_animation_pvd(comp, dir, output_name)
end

function test_svd(anm, k_range, x_qbits_range, d_qbits_range,
                  output_filename, output_directory = "~/data/logs")
    original = anm.data

    open(joinpath(ensure_dir(output_directory), output_filename), "w") do csv
        write(csv,
              "K, BLOCK_QBITS, DICT_QBITS, BPVF, SVD_KG, SVD_TIME\n")
        for
            k in collect(k_range),
            x_qbits in collect(x_qbits_range),
            d_qbits in collect(d_qbits_range)
            encoded, time =
                @timed batch_enc(original,
                                 method = :svd,
                                 k = k,
                                 x_qbits = x_qbits,
                                 d_qbits = d_qbits)
            # error computation
            svd_kg = kg_error(encoded, original)
            bpvf = bpvf_temp(size(anm, :frames), k, anm.vertex_count,
                             x_qbits, d_qbits)
            write(csv, "$k, $x_qbits, $d_qbits, $bpvf, $svd_kg, $time\n")
        end
    end
end

# VR_PCA on X with partial knowledge of it
function compress_vr_pca_X_part(orig, output_name, dir = "~/data/exported";
                                k = 500, x_qbits = 12, d_qbits = 12,
                                vr_pca_iters = 1, X_fraction = 1.0,
                                X_sampling = :random, vr_pca_fraction = 0.002,
                                vr_pca_eta = 1000, vr_pca_sampling = :uniform)
    comp = deepcopy(orig)
    comp.data, time =
        @timed batch_enc_partial(orig.data,
                                 method = :vr_pca,
                                 k = k,
                                 x_qbits = x_qbits,
                                 d_qbits = d_qbits,
                                 X_fraction = X_fraction,
                                 X_sampling = X_sampling,
                                 vr_pca_iters = vr_pca_iters,
                                 vr_pca_fraction = vr_pca_fraction,
                                 vr_pca_eta = vr_pca_eta,
                                 vr_pca_sampling = vr_pca_sampling)
    add_data!(comp, :vertex, "distance",
              vertex_distances(comp.data - orig.data), 1)
    println("VR_PCA KG error = $(kg_error(comp.data, orig.data))")
    export_animation_pvd(comp, dir, output_name)
end

function test_X_partial(anm, k_range, x_qbits_range, d_qbits_range,
                        X_fraction_range, vr_pca_iters_range,
                        vr_pca_fraction_range, vr_pca_eta_range,
                        output_filename, output_directory = "~/data/logs")
    original = anm.data

    open(joinpath(ensure_dir(output_directory), output_filename), "w") do csv
        write(csv,
              "K, BLOCK_QBITS, DICT_QBITS, BPVF, ITERS, FRACTION, ",
              "ETA, KG, TIME\n")
        for k in collect(k_range),
            x_qbits in collect(x_qbits_range),
            d_qbits in collect(d_qbits_range),
            X_fraction in collect(X_fraction_range),
            vr_pca_iters in collect(vr_pca_iters_range),
            vr_pca_fraction in collect(vr_pca_fraction_range),
            vr_pca_eta in collect(vr_pca_eta_range)
            encoded, time =
                @timed batch_enc_partial(original,
                                         method = :vr_pca,
                                         k = k,
                                         x_qbits = x_qbits,
                                         d_qbits = d_qbits,
                                         X_fraction = X_fraction,
                                         X_sampling = :random,
                                         vr_pca_iters = vr_pca_iters,
                                         vr_pca_fraction = vr_pca_fraction,
                                         vr_pca_eta = vr_pca_eta,
                                         vr_pca_sampling = :uniform)
            # error computation
            vr_pca_kg = kg_error(encoded, original)
            bpvf = bpvf_temp(size(anm, :frames), k, anm.vertex_count,
                             x_qbits, d_qbits)
            write(csv,
                  "$k, $x_qbits, $d_qbits, $bpvf, $vr_pca_iters, ",
                  "$vr_pca_fraction, $vr_pca_eta, $vr_pca_kg, $time\n")
        end
    end
end

function two_block_enc(X; method = :svd, k = 1, x_qbits = 12, d_qbits = 12,
                       block_fraction = 0.5)
    xr = size(X, 1)
    ci = round(Int, block_fraction*xr/3)
    blocks = [X[1:ci*3, :], X[ci*3+1:end, :]]
    bcount = size(blocks, 1)

    if method == :svd
        dicts = map(b -> svd(b*b')[1][:, 1:k], blocks)
    end

    recon_blocks = map((b, dict) -> compress(b, dict, x_qbits, d_qbits),
                       blocks, dicts)
    return reduce(vcat, recon_blocks)
end

function test_variable_block_size(anm, k_range,
                                  block_fraction_range, output_filename,
                                  output_directory = "~/data/logs")

    original = anm.data
    open(joinpath(ensure_dir(output_directory), output_filename), "w") do csv
        write(csv, "K, B0, B1, BPVF, ITERS, FRACTION, KG, TIME\n")
        x_qbits = 12
        d_qbits = 12
        for k in collect(k_range),
            block_fraction in collect(block_fraction_range)
            encoded, time = @timed two_block_enc(original,
                                                 method = :svd,
                                                 k = k,
                                                 x_qbits = x_qbits,
                                                 d_qbits = d_qbits,
                                                 block_fraction = block_fraction)
            # error computation
            b0 = round(Int, block_fraction * size(original, 1) / 3)
            b1 = round(Int, (1-block_fraction) * size(original, 1) / 3)
            kg = kg_error(encoded, original)[1]
            bpvf = bpvf_temp(size(anm, :frames), k, anm.vertex_count,
                             x_qbits, d_qbits)
            write(csv, "$k, $b0, $b1, $x_qbits, $d_qbits, $bpvf, $kg, $time\n")
        end
    end
end
