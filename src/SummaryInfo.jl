"""
    richness(data::PopObj)
Calculates various allelic richness and returns vector of per-locus allelic richness.
To be called internally by functions calculating overall or per-population richness
either rarefied or not.

"""
function richness(data::PopObj)
    rich = map(i -> Iterators.flatten(i |> skipmissing) |> collect |> unique |> length, eachcol(data.loci, false))
    return DataFrame(locus = names(data.loci), richness = rich)
end
