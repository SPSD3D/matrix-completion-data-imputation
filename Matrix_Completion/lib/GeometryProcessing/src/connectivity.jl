"A Connectivity consists of adjacency and graph degree matrix for its vertices,
as well as a connectivity mapping of each vertex to its neighbours."
mutable struct Connectivity
    adjacency::SparseMatrixCSC{Int, Int}
    degree::Vector{Int}
    neighbours::AVector{Int}
end

function show(io::IO, conn::Connectivity)
    println(io, summary(conn))
    println(io, "Adjacency matrix")
    println(io, conn.adjacency)
    println(io, "Degree matrix")
    println(io, conn.degree)
    println(io, "Neighbours")
    println(io, conn.neighbours)
end

length(conn::Connectivity) = length(conn.degree)
getindex(conn::Connectivity, vi::Int) = conn.neighbours[vi]
function getindex(conn::Connectivity, spec)
    map(vi -> conn.neighbours[vi], collect(1:length(conn))[spec])
end
endof(conn::Connectivity) = endof(conn.degree)

function Connectivity(fv_inds::LArray, vcnt::Int = maximum(fv_inds.array))
    A = spzeros(Int, vcnt, vcnt)
    for face in cols(arrange(fv_inds, [:v_inds, :faces]))
        i, j, k = face[1], face[2], face[3]
        A[i, j] = 1; A[j, k] = 1; A[i, k] = 1
        A[j, i] = 1; A[k, j] = 1; A[k, i] = 1
    end
    D = sum(A, 2)[:]
    N = AVector(map(find, cols(A)))
    Connectivity(A, D, N)
end

"Return the sparse laplacian of the requested kind.

Available kinds are: :d_a for the (D - A), :i_da for the (I - D^-1 * A).

If the optional argument `vinds` is supplied, it signifies the vertices for
which the laplacian will be computed."
function laplacian(conn::Connectivity, kind::Symbol, vinds = :)
    A = conn.adjacency[vinds, vinds]
    D = sum(A, 2)[:]
    if     kind == :d_a  return spdiagm(D) - A
    elseif kind == :i_da return speye(A) - spdiagm(1 ./ D) * A
    end
end
