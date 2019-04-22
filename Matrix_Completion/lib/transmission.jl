"Accepts :uniform and :lloyd_max as kinds of quantizers."
function make_quantizer(array::Array, bits::Int; kind = :uniform)
    if kind == :uniform
        make_uniform_quantizer(array, bits)
    elseif kind == :lloyd_max
        make_lloyd_max_quantizer(array, bits)
    end
end

# For a quantization scheme with n levels, the input levels array contains n+1
# values, of which the first and the last are the minimum and maximum input
# value to the quantizer respectively. Centroids, if not supplied are computed
# by compute_centroids.
function make_quantizer(levels::Vector,
                        centroids::Vector = compute_centroids(levels))
    q(array::Array) = map(q, array)
    q(val) = searchsortedfirst(levels[2:end], val) - 1
    dq(array::Array) = map(dq, array)
    dq(lvl) = centroids[lvl+1]
    q, dq
end

# Compute quantization distortion as ||restored - original|| / ||original||.
distortion(array, q, dq) = vecnorm(array .- dq(q(array))) / vecnorm(array)

# Compute centroids (representation values for levels) as the mean of adjacent
# level boundaries (length(centroids) = length(levels) - 1).
function compute_centroids(levels::Vector)
    interval_count = length(levels) - 1
    centroids = zeros(interval_count)
    for li in 1:interval_count
        centroids[li] = (levels[li] + levels[li+1]) / 2
    end
    centroids
end

function make_uniform_quantizer(array, bits)
    lvl_count = 2 ^ bits
    min_val, max_val = extrema(array)
    step = (max_val - min_val) / lvl_count
    # return a quantizer and a dequantizer
    q(array::Array) = map(q, array)
    q(val) = floor(Int, (clamp(val, min_val, max_val) - min_val) / step)
    dq(array::Array) = map(dq, array)
    dq(val) = min_val + (0.5 + round(Int, clamp(val, 0, lvl_count-1))) * step
    return q, dq
end

function make_lloyd_max_quantizer(array, bits;
                                  tolerance = 1e-7, kind = :centroid)
    min_val, max_val = extrema(array)
    drange = max_val - min_val
    level_count = 2 ^ bits
    adjusted_tolerance = max(abs(tolerance), abs(max_val) * eps())
    levels = (collect(0:level_count)/level_count) * drange + min_val
    centroids = compute_centroids(levels)
    q, dq = make_quantizer(levels, centroids)
    old_distortion = distortion(array, q, dq)
    distortion_gain = old_distortion

    while distortion_gain > adjusted_tolerance
        # compute levels from centroids
        for li in 2:level_count
            levels[li] = (centroids[li-1] + centroids[li]) / 2
        end
        # compute centroids from sample means
        sample_counts = zeros(level_count)
        sample_sums = zeros(level_count)
        for val in array
            li = q(val) + 1
            sample_sums[li] += val
            sample_counts[li] += 1
        end
        for li in 1:level_count
            if sample_counts[li] == 0
                centroids[li] = (levels[li] + levels[li+1]) / 2
            else centroids[li] = sample_sums[li] / sample_counts[li]
            end
        end
        # compute distortion gain
        new_distortion = distortion(array, q, dq)
        distortion_gain = abs(new_distortion - old_distortion)
        old_distortion = new_distortion
    end
    q, dq
end

function transmit(data, dict, data_qbits = 12, dict_qbits = 12;
                  encode = (data, dict) -> dict' * data,
                  decode = (data, dict) -> dict  * data)
    enc_data = encode(data, dict)
    data_q, data_dq = make_quantizer(enc_data, data_qbits)
    dict_q, dict_dq = make_quantizer(dict, dict_qbits)
    rec_data = data_dq(data_q(enc_data))
    rec_dict = dict_dq(dict_q(dict))
    decode(rec_data, rec_dict)
end
