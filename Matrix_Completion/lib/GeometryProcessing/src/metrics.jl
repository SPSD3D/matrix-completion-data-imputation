"Compute the KG error. Arguments original and reconstructed must have labels
:coors, :vertices and :frames."
function kg_error(original::LArray, reconstructed::LArray; db = true)
    D_orig = arrange(original, [[:coors, :vertices], :frames])
    D_recon = arrange(reconstructed, [[:coors, :vertices], :frames])
    numerator = vecnorm(D_orig - D_recon)
    recon_mean = sum(D_recon, 2) ./ size(recon, :frames)
    denominator = vecnorm(D_orig .- recon_mean)
    error = 100 * numerator / denominator
    if db error = 10 * log10(error) end
    error
end

"Compute hybrid (including laplacian) NMSVE. Arguments original and
reconstructed must be LArrays with labels :vertices and :coors. "
function nmsve_hybrid_error(original::LArray, reconstructed::LArray,
                            conn::Connectivity; db = true, meshwise = false,
                            cartesian_weight = 0.5, normalized = true)
    function GL(data::LArray, vi::Int)
        vertex = arrange(data[(:vertices, vi)], [:coors])
        neighbours = arrange(data[(:vertices, conn[vi])], [:coors, :vertices])
        inv_dists = mapslices(v -> 1/vecnorm(v), neighbours .- vertex, 1)
        inv_scaled_neighbours = neighbours .* inv_dists
        vertex, vertex - sum(inv_scaled_neighbours, 2) / sum(inv_dists)
    end
    cerr = 0.0; lerr = 0.0      # cartesian and laplacian error
    for vi in 1:size(original, :vertices)
        ov, oGL = GL(original, vi); rv, rGL = GL(reconstructed, vi)
        cdist = vecnorm(ov - rv); ldist = vecnorm(oGL - rGL)
        cerr += meshwise ? cdist^2 : cdist
        lerr += meshwise ? ldist^2 : ldist
    end
    if meshwise cerr = sqrt(cerr); lerr = sqrt(lerr) end
    error = 2 * (cartesian_weight * cerr + (1-cartesian_weight) * lerr)
    if normalized error /= 2 * size(original, :vertices) end
    if db error = 10 * log10(error) end
    error
end

function mean_nmsve_hybrid_error(original::Animation, reconstructed::Animation;
                                 frames = :, db = true, meshwise = false,
                                 cartesian_weight = 0.5, normalized = true)
    mean(map(frame -> nmsve_hybrid_error(original.data[(:frames, frame)],
                                         reconstructed.data[(:frames, frame)],
                                         original.connectivity,
                                         db = db,
                                         meshwise = meshwise,
                                         cartesian_weight = cartesian_weight,
                                         normalized = normalized),
             (1:size(original, :frames))[frames]))
end

"Compute the distance error between corresponding vertices"
function vertex_error(original::LArray, reconstructed::LArray)
    frame_count = size(original, :frames)
    vertex_count = size(original, :vertices)
    diff = zeros(frame_count, vertex_count)
    for fi in 1:frame_count
        for vi in 1:vertex_count
            orig_v = arrange(original[(:frames, fi), (:vertices, vi)],
                             [:coors])
            recon_v = arrange(reconstructed[(:frames, fi), (:vertices, vi)],
                              [:coors])
            diff[fi, vi] = vecnorm(orig_v - recon_v)
        end
    end
    diff
end
