locations(data::PopData)


function locations!(data::PopData, lat::Vector{Union{Missing,T}},
function locations!(data::PopData, lat::Vector{T}, long::Vector{T})
function locations!(data::PopData, lat::Vector{Union{Missing,String}}, long::Vector{Union{Missing,String}})
function locations!(
    data::PopData,
    lat_deg::Vector{Union{Missing,Int}},
    lat_min::Vector{Union{Missing,Float64}},
    long_deg::Vector{Union{Missing,Int}},
    long_min::Vector{Union{Missing,Float64}},
)
function locations!(data::PopData; kwargs...)

function loci(data::PopData)
function loci(data::IndexedTable)
function locus(data::PopData, locus::String)
function meta(data::PopData)
@inline function Base.missing(data::PopData; mode::String = "sample")
@inline function populations(data::PopData; listall::Bool = false)
function populations!(data::PopData, rename::Dict)
function populations!(data::PopData, rename::Vector{String})
function populations!(data::PopData, rename::NamedTuple)
function populations!(data::PopData, rename::Vector{String}, function populations!(data::PopData, oldnames::Vector{String},
function exclude_loci(data::PopData, locus::String)
function exclude_loci(data::PopData, exloci::Vector{String})
function exclude_samples(data::PopData, samp_id::String)
function exclude_samples(data::PopData, samp_ids::Vector{String})
function samples(data::PopData)
function popdata_map_helper(d, x)
function popdata_where_helper(args...)
macro show_only(args...)
