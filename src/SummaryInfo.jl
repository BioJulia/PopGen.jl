"""
    allele_avg(data::PopData, rounding::Bool = true)
Returns a NamedTuple of the average number of alleles ('avg') and standard
deviation (`stdev`) of a `PopData`. Use `false` as second argument (no keyword)
to not round results. Default (`true`) rounds to 4 digits.
"""
function allele_avg(data::PopData, rounding::Bool = true)
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
    richness(data::PopData)
Calculates various allelic richness and returns vector of per-locus allelic richness.
To be called internally by functions calculating overall or per-population richness
either rarefied or not.
"""
function richness(data::PopData)
    # collapse genotypes into array of alleles |> count number of unique alleles
    rich = map(eachcol(data.loci, false)) do i
        Iterators.flatten(i |> skipmissing) |> collect |> unique |> length
    end
    return DataFrame(locus = names(data.loci), richness = rich)
end


"""
    summary(data::PopData)
Prints a summary of the information contained in a PopData
"""
function Base.summary(data::PopData)
    println("PopData Object")
    if typeof(collect(skipmissing(data.loci.columns.genotype)[1][1])) == Int16
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    print("  Marker type: "); printstyled(marker, "\n", bold = true)
    ploidy = unique(data.samples.columns.ploidy) |> sort
    if length(ploidy) == 1
        print("  Ploidy: ") ; printstyled(ploidy |> join, "\n", bold = true)
    else
        print("  Ploidy (varies): ")
        print(ploidy[1]), [print(", $i") for i in ploidy[2:end]]
    end
    print("  Number of individuals: ") ; printstyled(length(data.samples.columns.name), "\n", bold = true)
    print("  Number of loci: ") ; printstyled(length(levels(data.loci.columns.locus)), "\n", bold = true)
    print("  Populations: ") ; printstyled(length(unique(data.samples.columns.population)), "\n", bold = true)

    if ismissing.(data.samples.columns.longitude) |> all == true
        print("  Longitude:") ; printstyled(" absent\n", color = :yellow)
    else
        println("  Longitude: present with ", count(i -> i === missing, data.samples.columns.longitude), " missing")
    end
    if ismissing.(data.samples.columns.longitude) |> all == true
        print("  Latitude:") ; printstyled(" absent\n", color = :yellow)
    else
        println("  Latitude: present with ", count(i -> i === missing, data.samples.columns.latitude), " missing")
    end
end
