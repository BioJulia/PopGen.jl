"""
    PopObj(samples::DataFrame, loci::DataFrame)
- `samples` ::DataFrame individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names/numbers
    - `ploidy` ::Int64 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `loci` ::DataFrame loci and their genotypes
    - genotypes are Tuples of integers, arraged in order if `.sample.name`
"""
mutable struct PopObj
    samples::DataFrame
    loci::DataFrame
end
