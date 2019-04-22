"Animation is a data structure describing an animation in space and time.

data::LArray{Float, 3}: A 3D LArray with dimension configuration [:coors,
:vertices, :frames], containing x, y, z coordinates at each frame of each vertex
of the animated model.
   
vertex_data::Dict{Symbol, LArray} and face_data::Dict{Symbol, LArray} Two
dictionaries containing vertex and face data, respectively. Data is in the form
of an LArray with labels [:data, <element>] in the case of frame-invariant data
or [:data, <element>, :frames] otherwise (element {:vertices | :faces}).
   
fv_inds::LArray{Int, 2}: An LArray with dimension labels [:v_inds, :faces],
holding the vertex indices of each triangular face of the animated model
(assumed to be static).
   
connectivity::Connectivity: A connectivity encoding the neighbours of each
vertex of the animated model (the set of vertices with which it is connected via
edges)."
mutable struct Animation
    data::LArray{Float, 3}            # [:coors, :vertices, :frames]
    vertex_data::Dict{Symbol, LArray} # [:data, :vertices, :frames]
    face_data::Dict{Symbol, LArray}   # [:data, :faces, :frames]
    fv_inds::LArray{Int, 2}           # [:v_inds, :faces]
    connectivity::Connectivity
end

"Preset configurations for animation data matrices, and relevant extensions for
LArray functions permutedims, arrange and size."
function anm_conf(conf_name::Symbol)
    confs = Dict(:eg           => [[:coors, :frames], [:vertices]],
                 :eg_transpose => [[:vertices], [:coors, :frames]],
                 :frame_column => [[:coors, :vertices], [:frames]])
    return confs[conf_name]
end

function show(io::IO, anm::Animation)
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

permutedims(la::LArray, cname::Symbol) = permutedims(la, anm_conf(cname))

arrange(la::LArray, conf::Symbol) = arrange(la, anm_conf(conf))

size(anm::Animation) = size(anm.data)

function size(anm::Animation, label::Symbol)
    max(size(anm.fv_inds, label), size(anm.data, label))
end

# Connectivity and partitioning convenience functions
neighbours(anm::Animation, vi::Int) = anm.connectivity[vi]

function edge_vertices(part::Partition, anm::Animation)
    edge_vertices(part, anm.connectivity)
end

"Return the partition vertices that are connected to the outside for every
partition in the Partitioning provided"
function edge_vertices(prt::Partitioning, anm::Animation)
    map(part -> edge_vertices(part, anm), prt)
end

function Segmentation(anm::Animation, prt::Partitioning)
    Segmentation(anm.data, (:vertices, unflatten(prt)))
end

# Manipulation of animation (face and vertex) data
"Backend function to add data to the specified field of animation."
function add_data!(anm::Animation, field::Symbol, name::Symbol, data::LArray)
    if haskey(getfield(anm, field), name)
        println("Replacing data with name \"$name\". under $field")
        delete!(getfield(anm, field), name)
    end
    get!(getfield(anm, field), name, data)
end

"Fix the data supplied in an array according to the specified element count and
frame count of animation. 1D arrays are treated as frame-invariant scalar
data. 2D arrays are treated as frame-invariant vector data, except for the case
where the array dimensions match element and frame count of animation, where
they are treated as frame-changing scalar data. 3D arrays are treated as
frame-changing vector data."
function fix_data(anm::Animation, ar::Array, element::Symbol)
    ar_size = size(ar)
    elem_count = size(anm, element)
    frame_count = size(anm, :frames)
    # dimensions indices with sizes elem_count, frame_count and data dimensions
    elem_dim = findfirst(ar_size, elem_count)
    frame_dim = findfirst(ar_size, frame_count)
    data_dim = findfirst(x -> x!=elem_count && x!=frame_count, ar_size)
    # check that at least one dimension size matches element count in animation
    if elem_dim == 0 error("Invalid animation element count") end
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
        if frame_dim == 0 error("Invalid animation frame count") end
        fixed = permutedims(ar, (data_dim, elem_dim, frame_dim))
    end
    if ndims(fixed) == 2 return LArray(fixed, [:data, element])
    else return LArray(fixed, [:data, element, :frames])
    end
end

"Convenience function to fix (see fix_data) and add vertex data to animation."
function add_vertex_data!(anm::Animation, name::Symbol, data::Array)
    add_data!(anm, :vertex_data, name, fix_data(anm, data, :vertices))
end

"Convenience function to fix (see fix_data) and add face data to animation."
function add_face_data!(anm::Animation, name::Symbol, data::Array)
    add_data!(anm, :face_data, name, fix_data(anm, data, :faces))
end

"Remove vertex data from animation."
function remove_vertex_data!(anm::Animation, name::Symbol)
    delete!(anm.vertex_data, name)
end

"Remove face data from animation."
function remove_face_data!(anm::Animation, name::Symbol)
    delete!(anm.face_data, name)
end

# Animation import (constructor) and export
function Animation(filename::String)
    matched = match(Regex("(?=\.)[0-9a-z]+\$"), filename)
    if !isa(matched, Void)
        ext = matched.match
        ext == "obj" && return import_animation_obj([filename])
    else
        println("Unsupported file extension.")
    end
end

function Animation(directory::String, name::String, ext::String)
    fileseries = find_files(directory, name, ext)
    ext == "obj" && return import_animation_obj(fileseries)
    println("Unsupported file extension.")
end

function import_animation_obj(fileseries::Vector{String})
    print("Importing...")
    frame_count  = length(fileseries)

    # determine the number of faces and vertices from the first file
    face_count = 0
    vertex_count = 0
    for line in eachline(fileseries[1])
        startswith(strip(line), "f ") && (face_count += 1)
        startswith(strip(line), "v ") && (vertex_count += 1)
    end
    vertex_data = Array{Float, 3}(3, vertex_count, frame_count)

    # process the files for vertices
    fi = 0                      # frame index
    for filename in fileseries
        li = 0; vi = 0; fi += 1 # line and vertex index
        for line in eachline(filename)
            li += 1
            # split tokens
            line = strip(line)
            (startswith(line, "#") || isempty(line)) && continue
            tokens = split(line)
            # line type dispatch
            line_type = shift!(tokens)
            if line_type == "v"
                vi += 1
                if length(tokens) < 3
                    error("Incomplete line $(fileseries[fi]):$li")
                elseif vi > vertex_count
                    error("Extra vertices $(fileseries[fi]):$li")
                end
                vertex_data[:, vi, fi] = map(t -> parse(Float, t), tokens[1:3])
            end
        end
        vi < vertex_count && error("Missing vertices $(fileseries[fi]):$li")
    end

    # collect 3 by face_count matrix containing the triangular face indices
    face_data = Matrix{Int}(3, face_count)
    li = 0; fi = 0;             # line and face index
    for line in eachline(fileseries[1])
        li += 1
        line = strip(line)
        (startswith(line, "#") || isempty(line)) && continue
        tokens = split(line)
        # line type dispatch
        line_type = shift!(tokens)
        if line_type == "f"
            fi += 1
            if length(tokens) != 3
                error("Invalid face $(fileseries[fi]):$li")
            elseif fi > face_count
                error("Extra faces $(fileseries[fi]):$li")
            end
            face_data[:, fi] =
                map(t -> (sep = search(t, '/');
                          sep == 0 ? round(Int, parse(t)) : round(Int, parse(t[1:sep-1]))),
                    tokens[1:3])
        end
    end
    println("done.")

    # add labels to data
    vertex_data = LArray(vertex_data, [:coors, :vertices, :frames])
    face_data = LArray(face_data, [:v_inds, :faces])

    # connectivity computation
    conn = Connectivity(face_data, vertex_count)

    Animation(vertex_data,
              # empty dictionaries for vertex and face data
              Dict{Symbol, LArray}(),
              Dict{Symbol, LArray}(),
              face_data,
              conn)
end

"Export the animation in PVD format under the specified directory and name."
function export_animation_pvd(anm::Animation, dir::String, name::String)
    dir   = ensure_dir(dir)
    pvd   = paraview_collection(joinpath(dir, name))
    pad   = round(Int, log10(size(anm, :frames)) + 1)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, fv_inds)
             for fv_inds in cols(arrange(anm.fv_inds, [:v_inds, :faces]))]
    # following one-based indexing for output filenames
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
        collection_add_timestep(pvd, vtkfile, convert(Float, fmi))
    end
    vtk_save(pvd)
end

"Export the animation in OBJ format under the specified directory and name."
function export_animation_obj(anm::Animation, dir::String, name::String)
    dir = ensure_dir(dir)
    pad = round(Int, log10(size(anm, :frames)) + 1)

    # following one-based indexing for output filenames
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
