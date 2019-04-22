module GeometryCompression

using Meshes
using WriteVTK

include("utilities.jl")
include("larray.jl")
include("animation.jl")

export Animation, anm_conf, add_data!, get_conf, import_animation,
       export_animation_pvd, export_animation_obj

export LArray, arrange

end
