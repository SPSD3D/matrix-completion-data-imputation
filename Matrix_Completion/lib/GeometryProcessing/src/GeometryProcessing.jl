module GeometryProcessing

if Int === Int32
    const Float = Float32
else
    const Float = Float64
end

using WriteVTK

importall Base
importall Base.Iterators

include("utilities.jl")
include("transmission.jl")
include("file-processing.jl")
include("avector.jl")
include("larray.jl")
include("connectivity.jl")
include("partitioning.jl")
include("animation.jl")
include("metrics.jl")

export
    Int, Float,
    # utilities.jl
    chain_reduce, chain_map,
    # transmission.jl
    make_quantizer, transmit,
    # file-processing.jl
    ensure_dir, find_files,
    # avector.jl
    AVector, unflatten,
    # larray.jl
    LArray, Segmentation, divide, arrange_with_inverse, arrange,
    # connectivity.jl
    Connectivity, laplacian,
    # partitioning.jl
    Partition, Partitioning, vertex_partition, partition_vertices,
    edge_vertices, extend,
    # animation.jl
    Animation, anm_conf, add_data!, get_conf, import_animation, fix_data,
    export_animation_pvd, export_animation_obj, add_vertex_data!, add_face_data!,
    # metrics.jl
    kg_error, nmsve_hybrid_error, mean_nmsve_hybrid_error, vertex_error
end
