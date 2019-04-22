using Meshes
using WriteVTK

include("utilities.jl")
include("file-processing.jl")
include("larray.jl")
include("laplacian.jl")
include("partitioning.jl")

if !isdefined(:Animation)
    type Animation
        data::LArray{Float64, 3}          # [:coors, :vertices, :frames]
        vertex_data::Dict{Symbol, LArray} # [:data, :vertices, :frames]
        face_data::Dict{Symbol, LArray}   # [:data, :faces, :frames]
        fv_inds::LArray{Int64, 2}         # [:v_inds, :faces]
        connectivity::Connectivity
    end
end

# Preset configurations for animation data matrices, and relevant extensions for
# LArray functions permutedims, arrange and size.
function anm_conf(conf_name::Symbol)
    confs = Dict(:eg           => [[:coors, :frames], [:vertices]],
                 :eg_transpose => [[:vertices], [:coors, :frames]],
                 :frame_column => [[:coors, :vertices], [:frames]])
    return confs[conf_name]
end

function Base.show(io::IO, anm::Animation)
    println(io, summary(anm))
    println(io, "Animation data")
    println(io, anm.data)
    println(io, "Vertex data")
    println(io, anm.vertex_data)
    println(io, "\nFace data")
    println(io, anm.face_data)
    println(io, "\nFace vertex indices")
    println(io, anm.fv_inds)
    println(io, "Connectivity")
    print(io, anm.connectivity)
end

Base.permutedims(la::LArray, cname::Symbol) = permutedims(la, anm_conf(cname))

arrange(la::LArray, conf::Symbol) = arrange(la, anm_conf(conf))

Base.size(anm::Animation) = size(anm.data)
function Base.size(anm::Animation, label::Symbol)
    max(size(anm.fv_inds, label), size(anm.data, label))
end

# Connectivity and partitioning
function neighbours(anm::Animation, vi::Int64)
    neighbours(anm.connectivity, vi)
end

function adjacency_matrix(anm::Animation)
    adjacency_matrix(anm.fv_inds, size(anm, :vertices))
end

function edge_vertices(anm::Animation, prt::Partitioning, pi::Int64)
    edge_vertices(anm.connectivity, prt, pi)
end

function edge_vertices(anm::Animation, prt::Partitioning)
    map(pi -> edge_vertices(anm, prt, pi), 1:size(prt, :partitions))
end

# Backend function to add data to the specified field of animation.
function add_data!(anm::Animation, field::Symbol, name::Symbol, data::LArray)
    if haskey(getfield(anm, field), name)
        println("Replacing data with name \"$name\". under $field")
        delete!(getfield(anm, field), name)
    end
    get!(getfield(anm, field), name, data)
end

# Fix the data supplied in an array according to the specified element count and
# frame count of animation. 1D arrays are treated as frame-invariant scalar
# data. 2D arrays are treated as frame-invariant vector data, except for the
# case where the array dimensions match element and frame count of animation,
# where they are treated as frame-changing scalar data. 3D arrays are treated as
# frame-changing vector data.
function fix_data(anm::Animation, ar::Array, element::Symbol)
    ar_size = size(ar)
    elem_count = size(anm, element)
    frame_count = size(anm, :frames)
    # dimensions indices with sizes elem_count, frame_count and data dimensions
    elem_dim = findfirst(ar_size, elem_count)
    frame_dim = findfirst(ar_size, frame_count)
    data_dim = findfirst(x -> x!=elem_count && x!=frame_count, ar_size)
    # check that at least one dimension size matches element count in animation
    if (elem_dim == 0) error("Invalid animation element count") end
    # dispatch on input array dimensions
    if ndims(ar) == 1           # scalar data, frame-invariant
        fixed = ar'
    elseif ndims(ar) == 2
        if frame_dim != 0       # scalar data, frame-changing
            fixed = permutedims(ar, (elem_dim, frame_dim))
            fixed = reshape(fixed, 1, elem_count, frame_count)
        else                    # vector data, frame-invariant
            fixed = permutedims(ar, (data_dim, elem_dim))
        end
    elseif ndims(ar) == 3       # vector data, frame-changing
        if (frame_dim == 0) error("Invalid animation frame count") end
        fixed = permutedims(ar, (data_dim, elem_dim, frame_dim))
    end
    if (ndims(fixed) == 2) return LArray(fixed, [:data, element])
    else return LArray(fixed, [:data, element, :frames])
    end
end

function add_vertex_data!(anm::Animation, name::Symbol, data::Array)
    add_data!(anm, :vertex_data, name, fix_data(anm, data, :vertices))
end

function add_face_data!(anm::Animation, name::Symbol, data::Array)
    add_data!(anm, :face_data, name, fix_data(anm, data, :faces))
end

function remove_vertex_data!(anm::Animation, name::Symbol)
    delete!(anm.vertex_data, name)
end

function remove_face_data!(anm::Animation, name::Symbol)
    delete!(anm.face_data, name)
end

# Returns an Animation of frame meshes loaded from numbered files in specified
# directory with matching name and extension, sorted by the numbering just
# before the extension.
function Animation(directory::String, name::String, ext::String)
    println("Loading...")
    frames = map(load, find_files(directory, name, ext))
    println("done.")

    # collect animation data
    print("Constructing animation data...")
    frame_count  = length(frames)
    vertex_count = length(frames[1].vertices)
    dimensions   = length(frames[1].vertices[1])
    data         = Array{Float64, 3}(dimensions, vertex_count, frame_count)
    for (fi, frame) in enumerate(frames)
        for (vi, vertex) in enumerate(frame.vertices)
            for di in 1:dimensions
                data[di, vi, fi] = vertex[di]
            end end end
    println("done.")

    # collect 3 by face_count matrix containing the triangular face indices,
    # adjusted to native julia (and vtk, obj) one-based indexing.
    print("Collecting face indices...")
    faces      = frames[1].faces
    offset     = typeof(faces[1]).parameters[3]
    face_count = length(faces)
    fv_inds    = LArray(Array{Int64, 2}(3, face_count), [:v_inds, :faces])
    for fi in 1:face_count
        for vi in 1:3
            fv_inds[vi, fi] = faces[fi][vi] - offset
        end end
    println("done.")

    # connectivity computation
    connectivity = Connectivity(fv_inds, vertex_count)

    return Animation(LArray(data, [:coors, :vertices, :frames]),
                     # empty dictionaries for vertex and face data
                     Dict{String, Tuple{Array{Float64, 2}, Int64}}(),
                     Dict{String, Tuple{Array{Float64, 2}, Int64}}(),
                     fv_inds,
                     connectivity)
end

function export_animation_pvd(anm::Animation, dir::String, name::String)
    dir   = ensure_dir(dir)
    pvd   = paraview_collection(name)
    pad   = round(Int64, log10(size(anm, :frames)) + 1)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, fv_inds)
             for fv_inds in cols(arrange(anm.fv_inds, [:v_inds, :faces]))]
    # following zero-based indexing for output filenames
    for fmi in 1:size(anm, :frames)
        vertices = arrange(anm.data[(:frames, fmi)], [:coors, :vertices])
        vtk_filename = joinpath(dir, "$(name)_$(lpad(fmi, pad, 0))")
        vtkfile = vtk_grid(vtk_filename, vertices, cells)
        # vertex data
        for (name, data) in anm.vertex_data
            vdata = arrange(data[(:frames, fmi)], [:data, :vertices])
            vtk_point_data(vtkfile, vdata, string(name))
        end
        # face data
        for (name, data) in anm.face_data
            fc_data = arrange(data[(:frames, fmi)], [:data, :faces])
            vtk_cell_data(vtkfile, fc_data, string(name))
        end
        collection_add_timestep(pvd, vtkfile, convert(Float64, fmi))
    end
    vtk_save(pvd)
end

function export_animation_obj(anm::Animation, dir::String, name::String)
    dir = ensure_dir(dir)
    pad = round(Int64, log10(size(anm, :frames)) + 1)

    # following zero-based indexing for output filenames
    for fmi in 1:size(anm, :frames)
        vertices = arrange(anm.data[(:frames, fmi)], [:coors, :vertices])
        open(joinpath(dir, "$(name)_$(lpad(fmi, pad, 0)).obj"), "w") do f
            # export vertices
            for vcoors in cols(vertices)
                # write emits 'f' instead of 'e' between mantissa and exponent
                @printf(f, "v %e %e %e\n", vcoors[1], vcoors[2], vcoors[3])
            end
            # export faces
            for face in cols(arrange(anm.fv_inds, [:v_inds, :faces]))
                write(f, "f $(face[1]) $(face[2]) $(face[3])\n")
            end
        end
    end
end
