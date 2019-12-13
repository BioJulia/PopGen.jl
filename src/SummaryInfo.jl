"""
    allele_avg(data::PopObj, rounding::Bool = true)
Returns a NamedTuple of the average number of alleles ('avg') and standard
deviation (`stdev`) of a `PopObj`. Use `false` as second argument (no keyword)
to not round results. Default (`true`) rounds to 4 digits.
"""
function allele_avg(data::PopObj, rounding::Bool = true)
    num_alleles = map(i -> Iterators.flatten(i |> skipmissing) |> collect |> unique |> length, eachcol(data.loci, false))
    avg = mean(num_alleles)
    sd = std(num_alleles)
    if rounding == true
        return (avg = round(avg, digits = 4), stdev = round(sd, digits = 4))
    else
        return (avg = avg, stdev = sd)
    end
end


"""
    richness(data::PopObj)
Calculates various allelic richness and returns vector of per-locus allelic richness.
To be called internally by functions calculating overall or per-population richness
either rarefied or not.
"""
function richness(data::PopObj)
    # collapse genotypes into array of alleles |> count number of unique alleles
    rich = map(i -> Iterators.flatten(i |> skipmissing) |> collect |> unique |> length, eachcol(data.loci, false))
    return DataFrame(locus = names(data.loci), richness = rich)
end
