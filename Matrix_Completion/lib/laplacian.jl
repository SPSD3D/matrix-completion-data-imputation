include("avector.jl")
include("larray.jl")
include("utilities.jl")

if !isdefined(:Connectivity)
    typealias Connectivity AVector{Int64}
end

function adjacency_matrix(fv_inds::LArray)
    # assuming that vertices are indexed sequentially, compute the vertex count
    # as the maximum vertex index in fv_inds and call the more specific method.
    adjacency_matrix(fv_inds, maximum(fv_inds.array))
end

function adjacency_matrix(fv_inds::LArray, vertex_count::Int64)
    face_count = size(fv_inds, :faces)
    A = falses(vertex_count, vertex_count)
    for face in cols(arrange(fv_inds, [:v_inds, :faces]))
        i, j, k = face[1], face[2], face[3]
        A[i, j] = true; A[j, k] = true; A[i, k] = true
        A[j, i] = true; A[k, j] = true; A[k, i] = true
    end
    A
end

function Connectivity(fv_inds::LArray, vertex_count::Int64)
    AVector(map(find, cols(adjacency_matrix(fv_inds, vertex_count))))
end

neighbours(cnct::Connectivity, vi::Int64) = cnct[vi]
