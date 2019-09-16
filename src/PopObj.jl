"""
    PopObj(samples::DataFrame, loci::DataFrame)
- `samples` ::DataFrame individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names/numbers
    - `ploidy` ::Int8 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `loci` ::DataFrame loci and their genotypes
    - columns are named by loci
    - genotypes are Tuples of ::Int18, arraged in order of `.sample.name`
"""
mutable struct PopObj
    samples::DataFrame
    loci::DataFrame
end
