using MDCT

include("utilities.jl")
include("animation.jl")
include("vr-pca.jl")
include("isvd.jl")

function frame_block_svd(anm::Animation; d = 200, k = 50,
                         data_qbits = 12, dict_qbits = 12)
    blocks = divide(anm.data, :frames, d)
    flat_blocks, inflate = arrange_blocks(blocks, :eg)
    dictionaries = map(b -> svd(b * b')[1][:, 1:k], flat_blocks)
    recon_blocks = map((b, d) -> compress(b, d, data_qbits, dict_qbits),
                       flat_blocks,
                       dictionaries)
    recon_anm = deepcopy(anm)
    recon_anm.data = cat(:frames, inflate(recon_blocks))
    return recon_anm
end

function frame_block_aoi(anm::Animation; d = 200, k = 50,
                         data_qbits = 12, dict_qbits = 12)
    blocks = divide(anm.data, :frames, d)
    flat_blocks, inflate = arrange_blocks(blocks, :eg)
    dictionaries = chain_reduce(b -> svd(b * b')[1][:, 1:k],
                                (d, b) -> qr(b*b'*d[1:size(b, 1), :])[1],
                                flat_blocks)
    recon_blocks = map((b, d) -> compress(b, d, data_qbits, dict_qbits),
                       flat_blocks,
                       dictionaries)
    recon_anm = deepcopy(anm)
    recon_anm.data = cat(:frames, inflate(recon_blocks))
    return recon_anm
end

function frame_block_diff_dict(anm::Animation; d = 200, k = 50,
                               data_qbits = 12, dict_qbits = 12)
    blocks = divide(anm.data, :frames, d)
    flat_blocks, inflate = arrange_blocks(blocks, :eg)
    dictionaries = chain_reduce(b -> svd(b * b')[1][:, 1:k],
                                (d, b) -> qr(b*b'*d[1:size(b, 1), :])[1],
                                flat_blocks)
    # compute differential from each dictionary to the next
    dict_diffs = chain_map(identity, (d0, d1) -> d0' * d1, dictionaries)
    println("Dictionary differential sizes: $(map(size, dict_diffs)).")
    recon_dicts = chain_reduce(identity, (d, b) -> qr(pinv(d*d')*d*b)[1],
                               dict_diffs)
    println("Distances: $(map(norm, map(-, dictionaries, recon_dicts)))")
    recon_blocks = map((b, d) -> compress(b, d, data_qbits, dict_qbits),
                       flat_blocks,
                       recon_dicts)
    recon_anm = deepcopy(anm)
    recon_anm.data = cat(:frames, inflate(recon_blocks))
    return recon_anm
end

# Divide the extended (in order to be evenly symmetric at the edges) sequence to
# blocks of length d, which overlap by d/2, then apply MDCT to those blocks.
function seq_mdct{T <: Number}(seq::Vector{T}, d::Integer)
    if mod(d, 4) != 0 error("Block size must be a multiple of 4") end
    function mdct_blocks(vec, d)
        d2 = div(d, 2)
        map((i, t) -> vec[i: t],
            collect(1:d2:length(vec)-d+1),
            collect(d:d2:length(vec)))
    end
    d2 = div(d, 2)
    extended_seq = vcat(reverse(seq[1:d2]), seq, reverse(seq[end-d2+1:end]))
    map(mdct, mdct_blocks(extended_seq, d))
end

# Inverse the MDCT transforms of the overlaping blocks of sequence, return the
# reconstructed sequence.
function seq_imdct{T <: Number}(mdcts::Vector{Vector{T}})
    d2 = length(mdcts[1])
    if any(mdct -> length(mdct) != d2, mdcts)
        error("MDCTs must be of equal lengths")
    end
    imdcts = map(imdct, mdcts)
    top = imdcts[2:2:end]
    bot = imdcts[1:2:end]

    top = vcat(zeros(d2), top..., zeros(d2))
    bot = vcat(bot...)
    return (top + bot)[d2+1:end-d2]
end

function compress(B, D, data_qbits = 12, dict_qbits = 12;
                  encode = (B, D) -> D' * B,
                  decode = (B, D) -> D * B)
    B_enc = encode(B, D)        # encoding
    # construct quantizer and dequantizer functions for block and dictionary
    b_qzer, b_dqzer = make_linear_quantizer(B_enc, data_qbits)
    d_qzer, d_dqzer = make_linear_quantizer(D, dict_qbits)
    B_enc_q = b_qzer(B_enc)     # quantization
    D_q = d_qzer(D)
    B_enc_dq = b_dqzer(B_enc_q) # dequantization
    D_dq = d_dqzer(D_q)
    decode(B_enc_dq, D_dq)      # decoding
end

function make_linear_quantizer(array, bits)
    lvl_count = 2 ^ bits
    amin, amax = extrema(array)
    step = (amax - amin) / (lvl_count - 1)
    # return a quantizer and a dequantizer
    q(array::Array) = map(q, array)
    q(val) = round(Int, (clamp(val, amin, amax) - amin) / step)
    dq(array::Array) = map(dq, array)
    dq(val) = clamp(amin + step * round(Int, val), amin, amax)
    return q, dq
end

function arrange_blocks(blocks, conf)
    flat_blocks = []
    inv_transforms = []
    for block in blocks
        flat, inv = arrange_with_inverse(block, conf)
        push!(flat_blocks, flat)
        push!(inv_transforms, inv)
    end
    return flat_blocks, fbs -> map((inv, fb) -> inv(fb), inv_transforms, fbs)
end
