include("avector.jl")
include("file-processing.jl")
include("laplacian.jl")

if !isdefined(:Partitioning)
    type Partitioning
        vp_inds::Vector{Int64}
        partitions::AVector{Int64}
    end
end

function Base.show(io::IO, prt::Partitioning)
    println(io, summary(prt))
    println(io, "Vertex partition indices")
    println(io, prt.vp_inds)
    println(io, "Partition vertices")
    print(io, prt.partitions)
end

vertex_partition(prt::Partitioning, vi::Int64) = prt.vp_inds[vi]

partition_vertices(prt::Partitioning, pi::Int64) = prt.partitions[pi]

function Base.size(prt::Partitioning, element::Symbol)
    if element == :partitions
        length(prt.partitions.anchors) - 1
    elseif element == :vertices
        map(pi -> length(partition_vertices(prt, pi)),
            1:size(prt, :partitions))
    end
end

function Partitioning(vp_inds::Vector{Int64})
    min_pi, max_pi  = extrema(vp_inds)
    partition_count = max_pi - min_pi + 1
    vp_inds        -= min_pi - 1 # adjust to 1-based indexing
    partitions      = [Vector{Int64}() for i = 1:partition_count]
    for (vi, pi) in enumerate(vp_inds) push!(partitions[pi], vi) end
    Partitioning(vp_inds, AVector(partitions))
end

function Partitioning(filename::String)
    Partitioning(parse_simple_file(filename))
end

function is_edge(cnct::Connectivity, prt::Partitioning, vi::Int64)
    vp = vertex_partition(prt, vi)
    any(vi -> vertex_partition(prt, vi) != vp, neighbours(cnct, vi))
end

function edge_vertices(cnct::Connectivity, prt::Partitioning, pi::Int64)
    filter(vi -> is_edge(cnct, prt, vi), partition_vertices(prt, pi))
end
           
