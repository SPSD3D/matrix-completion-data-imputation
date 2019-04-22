# An anchored vector (AVector) consists of a data vector together with a vector
# of indices (anchors) which indicate the starting indices of segments on the
# data vector. It is meant as a faster and more compact replacement for a hash
# table, when data is static and segments are naturally addressed by integer
# indices. The anchors vector must cover the entire length of vector,
# i.e. anchors[1] = 1 and anchors[end] = length(vector)+1.
if !isdefined(:AVector)
    type AVector{T}
        anchors::Vector{Int64}
        vector::Vector{T}
    end
end

function Base.show(io::IO, av::AVector)
    println(io, summary(av))
    println(io, "Anchors")
    println(io, av.anchors)
    println(io, "Vector")
    print(io, av.vector)
end

function AVector{T <: Any}(segmented_data::Vector{Vector{T}})
    anchors = cumsum([1, map(length, segmented_data)...])
    AVector(anchors, reduce(vcat, segmented_data))
end

function Base.getindex(av::AVector, i::Int64)
    av.vector[av.anchors[i] : av.anchors[i+1]-1]
end
