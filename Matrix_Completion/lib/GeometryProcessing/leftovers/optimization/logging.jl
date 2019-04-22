using CSV

mutable struct ParameterSpace
    keys::Vector{Symbol}
    values::Vector{Any}
    null_key::Symbol
end

function ParameterSpace(parameters::Vector)
    keys = Vector{Symbol}()
    values = Vector{Any}()
    null_key = gensym()
    for p in parameters
        if isa(p, Tuple)
            ensure_type(p[1], Symbol, "Invalid parameter key.")
            push!(keys, p[1])
            push!(values, p[2])
        else
            push!(keys, null_key)
            push!(values, p)
        end
    end
    ParameterSpace(keys, values, null_key)
end

import Base.product
function Base.product(ps::ParameterSpace, flat = true)
    ps_prod = collect(Base.product(ps.values...))
    flat ? vec(ps_prod) : ps_prod
end

function ensure_type(value, DataType, message = "Invalid type")
    isa(value, DataType) ? value : error(message)
end

function format_calls(fn, parameter_specs...)
    format_calls(fn, ParameterSpace([parameter_specs...]))
end

function format_calls(fn, ps::ParameterSpace)
    # specialized function to pair ps.keys with values of a value combination
    function pair_keys(value_comb)
        pair(key, value) = key == ps.null_key ? value : Expr(:kw, key, value)
        map(pair, ps.keys, value_comb)
    end
    map(vcomb -> Expr(:call, fn, pair_keys(vcomb)...), Base.product(ps))
end

