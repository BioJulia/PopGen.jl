"""
    PopObj(ind::Array{String,1}, popid::Array{Union{Int64, String},1}, loci::Array{String,1}, ploidy::Int64, genotypes::Dict, NamedTuple{(:x, :y),Tuple{Float64,Float64}})
Type "PopObj", which stores population genetics genotype data
- `ind` ::Array{String,1} of individual names
- `popid` ::Array{Union{Int64,String},1} of population names/numbers
- `loci` ::Array{String,1} of locus names in order of appearance in `genotypes`
- `ploidy` ::Int64 single integer of ploidy
- `genotypes` ::Dict of [`ind`] => genotypes::Array{String,1} ordered by `loci`
- `locdata` ::NamedTuple of decimal degrees longitude (:x) and latitude (:y)
"""
mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genotypes::Dict
    longitude::Array{Union{Int64,Float64},1}
    latitude::Array{Union{Int64,Float64},1}
    #PopObj() = new(ind, popid, loci, ploidy, genotypes, longitude, latitude)
end
