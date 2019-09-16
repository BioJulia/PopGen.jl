"""
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
end
