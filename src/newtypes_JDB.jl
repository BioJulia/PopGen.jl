using JuliaDB

abstract type PopObj end

mutable struct PopSample <: PopObj
    name::String
    population::String
    ploidy::Int8
    longitude::Union{Missing, Float32}
    latitude::Union{Missing, Float32}
end

struct PopData <: PopObj
    samples::Vector{PopSample}
    loci::T where T<:IndexedTable
end
