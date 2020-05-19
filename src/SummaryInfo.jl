export allele_table, allele_avg, richness, summary
"""
    allele_table(data::PopData)
Return a "tidy" IndexedTable of the loci, their alleles, and their alleles' frequencies.
"""
@inline function allele_table(data::PopData)
    #@groupby data.loci (:locus, :population) flatten = true {allele = unique_alleles(:genotype)}
    frq = @groupby data.loci :locus flatten = true {freq = allele_freq(:genotype)}
    #transform(tmp, :frequency => frq.columns.freq)
end


"""
    allele_avg(data::PopData, rounding::Bool = true)
Returns a NamedTuple of the average number of alleles ('mean') and
standard deviation (`stdev`) of a `PopData`. Use `rounding = false` to
not round results. Default (`true`) roundsto 4 digits.
"""
function allele_avg(data::PopData; rounding::Bool = true, populations = false)
    tmp = richness(data)
    if !populations
        if rounding
            @with tmp {mean = round(mean(:richness), digits = 4), stdev = round(variation(:richness), digits = 4)}
        else
            @with tmp {mean = mean(:richness), stdev = variation(:richness)}
        end
    else
        if rounding
            @groupby tmp (:locus, :population) {mean = mean(:richness), stdev = round(variation(:richness), digits = 4)}
        else
            @groupby tmp (:locus, :population) {mean = mean(:richness), stdev = variation(:richness)}
        end
    end
end


"""
    richness(data::PopData)
Calculates various allelic richness and returns a table of per-locus
allelic richness. Use `populations = true` to calculate richness by
locus by population.
"""
function richness(data::PopData; populations::Bool = false)
    if !populations
        @groupby data.loci :locus {richness = length(unique_alleles(:genotype))}
    else
        @groupby data.loci (:locus, :population) {richness =  length(unique_alleles(:genotype))}
    end
end
