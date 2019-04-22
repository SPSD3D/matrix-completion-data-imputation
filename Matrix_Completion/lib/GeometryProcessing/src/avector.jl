"An anchored vector (AVector) consists of a data vector together with a vector
of indices (anchors) which indicate the starting indices of segments on the data
vector. It is meant as a faster and more compact replacement for a hash table,
when data is static and segments are naturally addressed by integer indices. The
anchors vector must cover the entire length of vector, i.e. anchors[1] = 1 and
anchors[end] = length(vector)+1."
mutable struct AVector{T}
    anchors::Vector{Int}
    vector::Vector{T}
end

function show(io::IO, av::AVector)
    println(io, summary(av))
    println(io, "Anchors")
    println(io, av.anchors)
    println(io, "Vector")
    print(io, av.vector)
end

length(av::AVector) = length(av.anchors) - 1
getindex(av::AVector, i::Int) = av.vector[av.anchors[i] : av.anchors[i+1]-1]
getindex(av::AVector, spec) = map(ind -> av[ind], collect(1:length(av))[spec])
endof(av::AVector) = endof(av.anchors) - 1

start(av::AVector) = 1
next(av::AVector, s::Int) = av[s], s+1
done(av::AVector, s::Int) = s > length(av)

function AVector{T <: Any}(segmented_data::Vector{Vector{T}})
    anchors = cumsum([1, map(length, segmented_data)...])
    AVector(anchors, reduce(vcat, segmented_data))
end

unflatten(av::AVector) = [av[i] for i in 1:length(av)]

