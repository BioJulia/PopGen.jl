export allele_table, allele_avg, richness, summary

#TODO fix this for DataFrames API or delete
"""
    allele_table(data::PopData)
Return a "tidy" IndexedTable of the loci, their alleles, and their alleles' frequencies.
"""
@inline function allele_table(data::PopData)
    frq = @groupby data.loci :locus flatten = true {freq = allele_freq(:genotype)}
end

"""
    allele_avg(data::PopData; rounding::Bool = true)
Returns a NamedTuple of the average number of alleles ('mean') and
standard deviation (`stdev`) of a `PopData`. Use `rounding = false` to
not round results. Default (`true`) rounds to 4 digits.
"""
function allele_avg(data::PopData; rounding::Bool = true, population = false)
    tmp = richness(data)
    if rounding
        (mean = round(mean(tmp.richness), digits = 4), stdev = round(variation(tmp.richness), digits = 4))
    else
        (mean = mean(tmp.richness), stdev = variation(tmp.richness))
    end
end


"""
    richness(data::PopData; population::Bool = false)
Calculates various allelic richness and returns a table of per-locus
allelic richness. Use `population = true` to calculate richness by
locus by population.
"""
function richness(data::PopData; population::Bool = false)
    if !population
        DataFrames.combine(
            groupby(data.loci, :locus),
            :genotype => (geno -> length(unique_alleles(geno))) => :richness
        )
    else
        DataFrames.combine(
            groupby(data.loci, [:locus, :population]),
            :genotype => (geno -> length(unique_alleles(geno))) => :richness
        )
    end
end
