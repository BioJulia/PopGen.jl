"""
    allele_table(data::PopData)
Return a "tidy" IndexedTable of the loci, their alleles, and their alleles' frequencies.
"""
@inline function allele_table(data::PopData)
    tmp = @groupby data.loci :locus flatten = true {allele = unique_alleles(:genotype)}
    frq = @groupby data.loci :locus flatten = true {freq = allele_freq(:genotype)}
    transform(tmp, :frequency => frq.columns.freq)
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
        @groupby data.loci :locus {richness = length(unique(alleles(:genotype)))}
    else
        @groupby data.loci (:locus, :population) {richness =  length(unique(alleles(:genotype)))}
    end
end


"""
    summary(data::PopData)
Prints a summary of the information contained in a PopData
"""
function Base.summary(data::PopData)
    println("PopData Object")
    if occursin("Int16", string(eltype(data.loci.columns.genotype)))
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    print("  Marker type: "); printstyled(marker, "\n", bold = true)
    ploidy = unique(data.meta.columns.ploidy) |> sort
    if length(ploidy) == 1
        print("  Ploidy: ") ; printstyled(ploidy |> join, "\n", bold = true)
    else
        print("  Ploidy (varies): ")
        print(ploidy[1]), [print(", $i") for i in ploidy[2:end]]
    end
    print("  Number of individuals: ") ; printstyled(length(data.meta.columns.name), "\n", bold = true)
    print("  Number of loci: ") ; printstyled(length(levels(data.loci.columns.locus)), "\n", bold = true)
    print("  Populations: ") ; printstyled(length(unique(data.meta.columns.population)), "\n", bold = true)

    if ismissing.(data.meta.columns.longitude) |> all == true
        print("  Longitude:") ; printstyled(" absent\n", color = :yellow)
    else
        println("  Longitude: present with ", count(i -> i === missing, data.meta.columns.longitude), " missing")
    end
    if ismissing.(data.meta.columns.longitude) |> all == true
        print("  Latitude:") ; printstyled(" absent\n", color = :yellow)
    else
        println("  Latitude: present with ", count(i -> i === missing, data.meta.columns.latitude), " missing")
    end
end
