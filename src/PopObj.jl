"""
<<<<<<< HEAD
    PopObj(samples::DataFrame, loci::DataFrame)
The data struct used for the PopGen population genetics ecosystem. You are
STRONGLY discouraged from manually creating dataframes to pass into a PopObj,
and instead should use the provided genepop, csv, or vcf file importers.

- `samples` ::DataFrame individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names/numbers
    - `ploidy` ::Int8 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `loci` ::DataFrame loci and their genotypes
    - columns are named by loci
    - genotypes are Tuples of ::Int16, arraged in order of `.samples.name`
"""
mutable struct PopObj
    samples::DataFrame
    loci::DataFrame
=======
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
>>>>>>> master
end
