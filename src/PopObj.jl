"""
    PopObj(ind::Array{String,1}, popid::Array{Union{Int64, String},1}, loci::Array{String,1}, ploidy::Int64, genotypes::Dict, NamedTuple{(:x, :y),Tuple{Float64,Float64}})
Type "PopObj", which stores population genetics genotype data
- `ind` ::Array{String,1} of individual names
- `popid` ::Array{Union{Int64,String},1} of population names/numbers
- `loci` ::Array{String,1} of locus names in order of appearance in `genotypes`
- `ploidy` ::Array{Int64,1} integers of ploidy in order of `ind`
- `genotypes` ::Dict of [`ind`] => genotypes::Array{String,1} ordered by `loci`
- `longitude` ::Array{Union{Int64,Float64},1} of longitude values
- `latitude` ::Array{Union{Int64,Float64},1} of latitude values
"""
mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Array{Int64,1}
    genotypes::Dict
    longitude::Array{Union{Float64},1}
    latitude::Array{Union{Float64},1}
end
