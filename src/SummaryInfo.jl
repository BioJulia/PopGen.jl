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
    rich = map(eachcol(data.loci, false)) do i
        Iterators.flatten(i |> skipmissing) |> collect |> unique |> length
    end
    return DataFrame(locus = names(data.loci), richness = rich)
end


"""
    summary(data::PopObj)
Prints a summary of the information contained in a PopObj
"""
function Base.summary(data::PopObj)

    println(" Object of type PopObj")
    if typeof(skipmissing(data.loci[!, :1])[1][1]) == Int16
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    println(" Marker type: ", marker)
    ploidy = unique(data.samples.ploidy) |> sort
    length(ploidy) == 1 && println(" Ploidy: ", ploidy |> join)
    if length(ploidy) != 1
        print(" Ploidy (varies): ")
        print(ploidy[1]), [print(", $i") for i in ploidy[2:end]]
    end
    println("\n Number of individuals: ", length(data.samples.name))
    println(" Number of loci: ", size(data.loci,2))
    if ismissing.(data.samples.longitude) |> all == true
        printstyled("\n Longitude: absent", color = :yellow)
    else
        println(" Longitude: present with ", count(i -> i === missing, data.samples.longitude), " missing")
    end
    if ismissing.(data.samples.longitude) |> all == true
        printstyled("\n Latitude: absent", color = :yellow)
    else
        println(" Latitude: present with ", count(i -> i === missing, data.samples.latitude), " missing")
    end
    print("\n\n Population names and counts:")
    DataFrames.show(populations(data), summary = false, rowlabel = :num, allrows = true)
end
