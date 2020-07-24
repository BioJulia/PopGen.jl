## test data available for use in PopGen
export dataset, nancycats, gulfsharks

"""
    dataset(::String)
Load an example dataset from either `"gulfsharks"` (SNP) or `"nancycats"` (microsatellite). Can also use `"sharks"` and `"cats"`
as shorthands. Use `?nancycats` and `?gulfsharks` to learn more about
these datasets.

### Example
```
ncats = dataset("nancycats")
gsharks = dataset("sharks")
```
"""
function dataset(name::String)
    lowercase(name) in ["gulfsharks", "sharks"] && return gulfsharks()
    lowercase(name) in ["nancycats", "cats"] && return nancycats()
end

"""
    nancycats()
Returns a `PopObj` of corresponding "nancycats" dataset as featured in
the R package `adegenet`. This is microsatellite data of 9 loci in 237
individuals across 17 populations.

Example:
```
ncats = nancycats()
```
"""
function nancycats()
    filename = normpath(joinpath(@__DIR__,"../..","data", "datasets.jld2"))
    load(filename, "nancycats")
end

"""
    gulfsharks()
Returns a `PopObj` corresponding the Blacknose shark dataset as used in
Dimens et al. 2019. This is a mid-sized SNV dataset of 2209 haplotypes 
in 212 individuals, across 7 populations.

Example:

```
sharks = gulfsharks()
```
"""
function gulfsharks()
    filename = normpath(joinpath(@__DIR__,"../..","data", "datasets.jld2"))
    load(filename, "gulfsharks")
end
