importall Base

const Label = Symbol
const LabelIndex = Tuple{Label, Any}

"A labeled array (LArray) is a simple type designed around a multidimensional
data array, whose purpose is to allow access through the dimension labels."
mutable struct LArray{T, N} <: DenseArray{T, N}
    array::Array{T, N}
    labels::Vector{Label}
end

"LArray constructors from the specified array, with the specified dimension
sizes and labels (together as returned by size, or split up)"
function LArray(a, sizes_labels::Tuple{Tuple, Vector{Label}})
    LArray(Array(reshape(a, sizes_labels[1])), sizes_labels[2])
end

function LArray(a, sizes::Tuple, labels::Vector{Label})
    LArray(Array(reshape(a, sizes)), labels)
end

function LArray(rv::RowVector, args...)
    LArray(Array(rv), args...)
end

# LArray utilities
"Return the number of the dimension named by the given label (zero if the label
is not present)."
dim_num(la::LArray, label::Label) = findfirst(la.labels, label)

function translate_labels(la::LArray, conf; signal_errors = false)
    map(conf) do elem
        if isa(elem, Label)
            dn = dim_num(la, elem)
            dn == 0 && signal_errors && error("Invalid label $elem")
            dn
        else translate_labels(la, elem, signal_errors = signal_errors)
        end
    end
end

function is_valid_conf(conf)
    if !all(el -> (isa(el, Label) || isa(el, Vector{Label})), conf)
        error("Invalid configuration: $conf.")
    end
end

function translate_indices(la::LArray, label_indices::Vararg{LabelIndex})
    la_indices = Vector(ndims(la)); fill!(la_indices, :)
    for (label, index) in label_indices
        label_dim = dim_num(la, label)
        if label_dim != 0 la_indices[label_dim] = index end
    end
    la_indices
end

# DenseArray interface: IndexStyle, size, getindex, setindex!
IndexStyle(::Type{<:LArray}) = IndexLinear()

function size(la::LArray, label::Label)
    dim = dim_num(la, label)
    dim == 0? 0 : size(la.array, dim)
end

size(la::LArray, dim::Int) = size(la.array, dim)
size(la::LArray) = size(la.array)

"Extension of getindex for LArray. Indices are supplied as tuples of labels and
index specifiers."
function getindex(la::LArray, label_indices::Vararg{LabelIndex})
    la_indices = translate_indices(la, label_indices...)
    # drop singleton dimensions from labels
    labels = la.labels[find(r -> !(typeof(r) <: Number), la_indices)]
    if labels == [] getindex(la.array, la_indices...)
    else LArray(getindex(la.array, la_indices...), labels)
    end
end

"Extension of setindex! for LArray. Indices are supplied as tuples of labels and
index specifiers."
function setindex!(la::LArray, x, label_indices::Vararg{LabelIndex})
    la_indices = translate_indices(la, label_indices...)
    setindex!(la.array, x, la_indices...)
end

# Transparent access to underlying array
getindex(la::LArray, inds...) = getindex(la.array, inds...)
setindex!(la::LArray, x, inds...) = setindex!(la.array, x, inds...)

"Permute dimensions of a labeled array according to the configuration supplied."
function permutedims(la::LArray, conf::Vector{Label})
    LArray(permutedims(la.array, translate_labels(la, conf)), conf)
end

permutedims(la::LArray, conf::Vector) = permutedims(la, vcat(conf...))

"Return a lower-dimensional layout of the data array together with a function
that reconstructs the given LArray from it. The conf argument is a vector of
Label and/or Vector{Label} (which dictate the order and fusion of dimensions in
the resulting array). For example the configuration [[:c, :a], [:b :d]] for an
la::LArray with la.labels = [:a, :b, :c] is equivalent to reshape(la.array,
size(la.array, 3)*size(la.array, 1), size(la.array, 2)). The superfluous :d
label is discarded, i.e. it's equivalent to the configurations [[:c, :a], [:b]]
and [[:c, :a], :b]."
function arrange_with_inverse(la::LArray, conf::Vector)
    # canonicalize configuration
    conf = map(el -> isa(el, Label) ? [el] : el, conf)
    # filter each configuration group to contain labels present in la
    map!(group -> filter(l -> l in la.labels, group), conf, conf)
    # store the flattened labels and permute la's dimensions accordingly
    conf_labels = vcat(conf...)
    permuted = permutedims(la, conf_labels)
    # number of la dimensions fused in each dimension of the result
    fused_dims = map(length, conf)
    perm_sizes = size(permuted.array)
    # arranged sizes are the products of the respective fused dimensions sizes
    arr_sizes = map(r -> prod(perm_sizes[r]), block_ranges(fused_dims))
    return (reshape(permuted.array, arr_sizes...),
            a -> permutedims(LArray(reshape(a, perm_sizes...), vcat(conf...)),
                             la.labels))
end

"Returns first value of arrange_with_inverse, called with the same arguments"
arrange(la::LArray, conf::Vector) = arrange_with_inverse(la, conf)[1]

# Extensions of other Base functions. Size returns the size of the data array
# along with the labels, or the size of one dimension, if specified.
function show(io::IO, la::LArray)
    println(io, summary(la))
    println(io, summary(la.array))
    println(io, la.array)
    println(io, la.labels)
end

function cat(label::Label, las::Vararg{LArray})
    labels = las[1].labels
    if !all(la -> la.labels == labels, las[2:end])
        error("LArrays to be concatenated should have the same labels.")
    end
    LArray(cat(findfirst(labels, label), map(la -> la.array, las)...), labels)
end

cat(label::Label, las::Vector{LArray}) = cat(label, las...)

# Return index bounds for consecutive ranges of the given sizes
function block_ranges(sizes)
    map((i, t) -> i:t, cumsum([1, sizes[1:end-1]...]), cumsum(sizes))
end

function block_ranges(la::LArray, label::Label, block_size::Int)
    block_count, last_block_size = divrem(size(la, label), block_size)
    block_sizes = fill(block_size, block_count)
    if (last_block_size != 0) push!(block_sizes, last_block_size) end
    block_ranges(block_sizes)
end

"Segmentation is an iterator which successively returns the segments of an
LArray according to the given range specifications."
mutable struct Segmentation
    la::LArray
    ranges::Vector
end

start(seg::Segmentation) = 1
next(seg::Segmentation, s::Int) = seg.la[seg.ranges[s]...], s+1
done(seg::Segmentation, s::Int) = s > length(seg.ranges)
length(seg::Segmentation) = length(seg.ranges)

function Segmentation(la::LArray, label_indices...)
    labels = map(r -> r[1], label_indices)
    # all possible combinations of dimension ranges (cartesian product)
    range_combs = collect(product(map(r -> r[2], label_indices)...))[:]
    # pair labels with each range combination
    Segmentation(la, map(rc -> collect(zip(labels, rc)), range_combs))
end

function Segmentation(la::LArray, label::Label, block_size::Int)
    Segmentation(la, (label, block_ranges(la, label, block_size)))
end
