## test data available for use in PopGen
export dataset, @nancycats, @gulfsharks

"""
    dataset(::String)
Load an example dataset from either `"gulfsharks"` (SNP) or `"nancycats"` (microsatellite). 
Can also use `"sharks"` and `"cats"` as shorthands. Use `?@nancycats` and 
`?@gulfsharks` to learn more about these datasets.

### Example
```
ncats = dataset("nancycats")
gsharks = dataset("sharks")
```
"""
function dataset(name::String)
    if lowercase(name) in ["nancycats", "cats"] 
        filename = normpath(joinpath(@__DIR__,"../..","data", "nancycats.gen"))
    elseif lowercase(name) in ["gulfsharks", "sharks"]
        filename = normpath(joinpath(@__DIR__,"../..","data", "gulfsharks.csv"))
    else
        throw(ArgumentError("Please choose either the \"nancycats\" or \"gulfsharks\" datasets"))
    end
    read_from(filename, silent = true)
end

"""
    @nancycats
Returns `PopData` of corresponding "nancycats" dataset as featured in
the R package `adegenet`. This dataset is composed of 9 microsatellite 
loci in 237 individuals across 17 populations.

Example:
```
ncats = @nancycats
```
"""
macro nancycats()
    return esc(quote
       dataset("nancycats")
    end)
end


"""
    @gulfsharks
Returns `PopData` corresponding the Blacknose shark dataset as used in
Dimens et al. 2019. This is a mid-sized SNP dataset of 2209 haplotypes 
in 212 individuals, across 7 populations.

Example:

```
sharks = @gulfsharks
```
"""
macro gulfsharks()
    return esc(quote
       dataset("gulfsharks")
    end)
end