"A Partitioning describes a partitioning of a set of vertices into subsets."
mutable struct Partitioning
    vp_inds::Vector{Int}
    partitions::AVector{Int}
end

const Partition = Vector{Int}

function show(io::IO, prt::Partitioning)
    println(io, summary(prt))
    println(io, "Vertex partition indices")
    println(io, prt.vp_inds)
    println(io, "Partition vertices")
    print(io, prt.partitions)
end

length(prt::Partitioning) = length(prt.partitions)
getindex(prt::Partitioning, spec) = prt.partitions[spec]
endof(prt::Partitioning) = endof(prt.partitions)

start(prt::Partitioning) = start(prt.partitions)
next(prt::Partitioning, s::Int) = next(prt.partitions, s)
done(prt::Partitioning, s::Int) = done(prt.partitions, s)

# partition_vertices supplied for convenience, same with getindex
vertex_partition(prt::Partitioning, vi::Int) = prt.vp_inds[vi]
partition_vertices(prt::Partitioning, vi::Int) = prt[vi]

function size(prt::Partitioning, element::Symbol)
    if element == :partitions
        length(prt.partitions.anchors) - 1
    elseif element == :vertices
        map(pi -> length(prt[pi]),
            1:size(prt, :partitions))
    end
end

function Partitioning(vp_inds::Vector{Int})
    min_pi, max_pi  = extrema(vp_inds)
    partition_count = max_pi - min_pi + 1
    vp_inds        -= min_pi - 1 # adjust to 1-based indexing
    partitions      = [Vector{Int}() for i = 1:partition_count]
    for (vi, pi) in enumerate(vp_inds) push!(partitions[pi], vi) end
    Partitioning(vp_inds, AVector(partitions))
end

function Partitioning(filename::String)
    Partitioning(parse_simple_file(filename))
end

unflatten(prt::Partitioning) = unflatten(prt.partitions)

"Return the partition vertices that are connected to the outside."
function edge_vertices(part::Partition, conn::Connectivity)
    outside = trues(length(conn)); outside[part] = false
    # scan for partition vertices connected to outside vertices
    part_cut = sum(conn.adjacency[part, outside], 2)[:]
    part[findn(part_cut)]
end

"Identify the vertices outside of the partition that are connected to the
partition vertices. Depending on the value of the `whole` keyword (true by
default), returns the extended partition (true) or only the extension vertices
(false)"
function extend(part::Partition, conn::Connectivity; whole = true)
    outside = trues(length(conn)); outside[part] = false
    outside = collect(1:length(conn))[outside]
    # scan for outside vertices connected to partition vertices
    outside_cut = sum(conn.adjacency[outside, part], 2)[:]
    # filter and return these vertices
    extension = outside[findn(outside_cut)]
    if whole vcat(part, extension) else extension end
end

"Same as the simple extend function, but apply it repeatedly."
function extend(part::Partition, conn::Connectivity, times::Int; whole = true)
    extension = Partition()
    for i in 1:times
        append!(extension, extend(vcat(part, extension), conn, whole = false))
    end
    if whole vcat(part, extension) else extension end
end
