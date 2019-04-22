# This file defines generic means of file processing aiming to provide a unified
# interface for dealing with batch file processing.

"Creates directory under given path if it doesn't exist."
function ensure_dir(path)
    dir = expanduser(path)
    if !isdir(dir) mkpath(dir) end
    return dir
end

"Parse the given simple file (has one element per line)."
function parse_simple_file(filename::String)
    open(expanduser(filename)) do file
        map(parse, readlines(file))
    end
end

"Return filenames under a directory matching extension and part of name."
function find_files(directory::String, name::String, ext::String)
    file_num(fname) = parse(match(Regex("[0-9]+(?=\.$ext\$)"), fname).match)
    directory = expanduser(directory)
    files = filter(file -> ismatch(Regex(".*$name.*[0-9]+\.$ext\$"), file),
                   readdir(directory))
    println("Found $(length(files)) files.")
    return map(fname -> joinpath(directory, fname), sort(files, by=file_num))
end

function process_files(files::Vector{String}, line_processors::Vararg{Function})
    results = map(f -> process_file(f, line_processors...), files)
    data_sizes = map(size, results[1])
    file_count = length(files)
    map((data, data_size) -> reshape(data, data_size..., file_count),
        map(hcat, results...),
        data_sizes)
end

function process_files(directory::String, name::String, ext::String,
                       line_processors::Vararg{Function})
    process_files(find_files(directory, name, ext), line_processors...)
end

"Apply the supplied line processors to each line of the specified file. A line
processor should accept a string (line) argument and return a vector of the
parsed results (maybe empty). Result vectors are horizontally concatenated to 2D
result matrices with rows equal to data entries per line and columns equal to
matching lines."
function process_file(filename::String, line_processors::Vararg{Function})
    lines = open(expanduser(filename)) do file
        map(strip, readlines(file))
    end
    results = map(lp -> map(lp, lines), [line_processors...])
    map(proc_res -> reduce(hcat, remove_if(isempty, proc_res)), results)
end

function vertex_color(line::String)
    l = split(line)
    if l[1] == "v" && length(l) >= 7
        map(t -> parse(Float, t), [l[5], l[6], l[7]])
    else []
    end
end

function vertex_coordinates(line::String)
    l = split(line)
    if l[1] == "v" && length(l) >= 4
        map(t -> parse(Float, t), [l[2], l[3], l[4]])
    else []
    end
end

# Dump matrix to file obeying matrix market format.
function dump_mpi_lasso_matrix(A, filename)
    m, n = size(A)
    open(expanduser(filename), "w") do f
        write(f, "%%MatrixMarket matrix array real general\n")
        write(f, "$m $n\n")
        for element in A[:]
            write(f, "$element\n")
        end
    end
end

# Dump vector to file obeying matrix market format.
function dump_mpi_lasso_vector(A, filename)
    m = length(A)
    open(expanduser(filename), "w") do f
        write(f, "%%MatrixMarket matrix array real general\n")
        write(f, "$m 1\n")
        for element in A[:]
            write(f, "$element\n")
        end
    end
end

