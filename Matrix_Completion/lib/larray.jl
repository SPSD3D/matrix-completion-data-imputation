# A labeled array (LArray) is a simple API designed around a multidimensional
# data array, the purpose of which is to allow manipulation through the
# dimension labels.
if !isdefined(:LArray)
    type LArray{T, N}
        array::Array{T, N}
        labels::Vector{Symbol}
    end
end

# LArray constructors from the specified array, with the specified dimension
# sizes and labels (together as returned by size, or split up)
function LArray(a::Array, sizes_labels::Tuple{Tuple, Vector{Symbol}})
    LArray(reshape(a, sizes_labels[1]), sizes_labels[2])
end

function LArray(a::Array, sizes::Tuple, labels::Vector{Symbol})
    LArray(reshape(a, sizes), labels)
end

# Permute dimensions of a labeled array according to the configuration supplied
function Base.permutedims(la::LArray, conf::Vector{Symbol})
    # Compute the permutation indices of new labels relative to the old
    function scan_permutation(old::Vector{Symbol}, new::Vector{Symbol})
        # ensure that conf is a permutation of la.labels
        if intersect(old, new) == old
            return map(l -> findfirst(old, l), new)
        else error("Provided label set different from current")
        end
    end
    permutation = scan_permutation(la.labels, conf)
    LArray(permutedims(la.array, permutation), conf)
end

Base.permutedims(la::LArray, conf) = permutedims(la, flatten(conf))

# Extension of getindex and setindex! for LArray. Ranges are specified as tuple
# pairs of labels and indices/ranges as for the procedures for simple arrays.
IRCV = Union{Int64, Range, Colon, Vector{Int64}, Vector{Bool}}
function adjust_ranges(la::LArray, ranges::Vararg{Tuple{Symbol, IRCV}})
    la_ranges = []
    for label in la.labels
        matches = filter(r -> r[1] == label, [ranges...])
        if matches != [] push!(la_ranges, matches[1][2])
        else push!(la_ranges, 1:size(la, label))
        end
    end
    la_ranges
end

function Base.getindex(la::LArray, rngs::Vararg{Tuple{Symbol, IRCV}})
    la_rngs = adjust_ranges(la, rngs...)
    # drop singleton dimensions from labels
    labels = la.labels[find(r -> !(typeof(r) <: Number), la_rngs)]
    if labels == [] getindex(la.array, la_rngs...)
    else LArray(getindex(la.array, la_rngs...), labels)
    end
end

function Base.setindex!(la::LArray, x::Array, rngs::Vararg{Tuple{Symbol, IRCV}})
    la_rngs = adjust_ranges(la, rngs...)
    setindex!(la.array, x, la_rngs...)
end

# Transparent access to underlying array
Base.getindex(la::LArray, inds...) = getindex(la.array, inds...)
Base.setindex!(la::LArray, x, inds...) = setindex!(la.array, x, inds...) 

# Divide the given la::LArray across the dimension specified by the supplied
# label into blocks of the given size(s)
function divide(la::LArray, label::Symbol, block_sizes::Vector{Int64})
    map((i, t) -> la[(label, i:t)], range_bounds(block_sizes)...)
end

function divide(la::LArray, label::Symbol, block_size::Int64)
    block_count, last_block_size = divrem(size(la, label), block_size)
    block_sizes = fill(block_size, block_count)
    if (last_block_size != 0) push!(block_sizes, last_block_size) end
    divide(la, label, block_sizes)
end

# Return a lower-dimensional layout of the data array together with a function
# that reconstructs the given LArray from it. The conf argument is a vector of
# Symbol and/or Vector{Symbol} (which dictate the order and fusion of dimensions
# in the resulting array). For example the configuration [[:c, :a], [:b :d]] for
# an la::LArray with la.labels = [:a, :b, :c] is equivalent to reshape(la.array,
# size(la.array, 3)*size(la.array, 1), size(la.array, 2)).  The superfluous :d
# label is discarded, i.e. it's equivalent to the configurations [[:c, :a],
# [:b]] and [[:ca, :a], :b].
function arrange_with_inverse(la::LArray, conf::Vector)
    conf = nest(conf)
    map!(group -> filter(l -> l in la.labels, group), conf)
    conf_labels = vcat(conf...)
    permuted = permutedims(la, conf_labels)
    fused_dims = map(length, conf)
    perm_sizes = size(permuted.array)
    rs_sizes = map((i, t) -> prod(perm_sizes[i:t]), range_bounds(fused_dims)...)
    return (reshape(permuted.array, rs_sizes...),
            a -> permutedims(LArray(reshape(a, perm_sizes...), vcat(conf...)),
                             la.labels))
end

arrange(la::LArray, conf) = arrange_with_inverse(la, conf)[1]

# Extensions of other Base functions. Size returns the size of the data array
# along with the labels, or the size of one dimension, if specified.
function Base.show(io::IO, la::LArray)
    println(io, summary(la))
    println(io, summary(la.array))
    println(io, la.array)
    println(io, la.labels)
end

function Base.size(la::LArray, label::Symbol)
    ind = findfirst(la.labels, label)
    return ind != 0 ? size(la.array, ind) : 0
end

function Base.size(la::LArray, labels::Vararg{Symbol})
    mapreduce(l -> size(la, l), *, 1, [labels...])
end

Base.size(la::LArray) = return size(la.array), la.labels
Base.size(la::LArray, d) = size(la.array, d)

function Base.cat{T, N <: Any}(label::Symbol, las::Vector{LArray{T, N}})
    labels = las[1].labels
    LArray(cat(findfirst(labels, label), map(la -> la.array, las)...), labels)
end

function Base.cat{T, N <: Any}(label::Symbol, las::Vararg{LArray{T, N}})
    cat(label, [las...])
end

# Return index bounds for consecutive ranges of the given sizes
range_bounds(sizes) = (cumsum([1, sizes[1:end-1]...]), cumsum(sizes))
# nest([:a, [:b, :c]]) => [[:a], [:b, :c]]
nest(conf) = map(el -> isa(el, Symbol) ? [el] : el, conf)
# flatten([:a, [:b, :c]]) => [:a, :b, :c]
flatten(conf) = vcat(conf...)
