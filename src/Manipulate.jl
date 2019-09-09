"""
    sample_names(x::PopObj)
View individual/sample names in a `PopObj`

Equivalent to `PopObj.samples.name`
"""
sample_names(x::PopObj) = x.samples.name

"""
    summary(x::PopObj)
Prints a summary of the information contained in a PopObj
"""
function summary(x::PopObj)
    println("Object of type PopObj:")
    println("\nLongitude:")
    println("$(x.samples.longitude[1:6])" , " \u2026 ", "$(x.samples.longitude[end-5:end])", "\n")
    println("Latitude:")
    println("$(x.samples.latitude[1:6])" , " \u2026 ", "$(x.samples.latitude[end-5:end])", "\n")
    println("Number of individuals: $(length(x.samples.name))")
    println(x.samples.name[1:3] , " \u2026 ", x.samples.name[end-2:end], "\n")
    println("Number of loci: $(size(x.loci,2))")
    println(string.(names(x.loci))[1:3], " \u2026 " , string.(names(x.loci))[end-2:end], "\n" )
    println("Ploidy:")
    println("$(x.samples.ploidy[1:6])", " \u2026 ", "$(x.samples.ploidy[end-7:end])")
    println("Number of populations: $(length(x.samples.population |> unique))","\n")
    println("#Inds | Pop","\n", "--------------")
    popcounts = hcat([sum(x.samples.population .== i) for i in unique(x.samples.population)],unique(x.samples.population))
    for eachpop in 1:length(popcounts)÷2
        println(popcounts[eachpop], " | ", popcounts[eachpop,2])
    end
    println("\nAvailable .samples fields: .name, .population, .ploidy, .longitude, .latitude")
end


"""
    loci(x::PopObj)
View the genotypes of all individuals for specific loci in a `PopObj`.
Default shows all genotypes for all individuals. Use `loci =` to specify a single
locus or array of loci to display.

`loci(wild_rice, "snp_451")`

`loci(wild_rice, ["snp_451", "snp_011", "snp_314"])`
"""
function loci(x::PopObj, loci::Union{String, Array, Nothing}= nothing)
    df = genotypes(x) ;
    if loci != nothing
        if typeof(loci) == String
            loci ∉ x.loci && error("locus $loci not found in PopObj")
            return df[!, [:ind, :population, Symbol(loci)]]
        else
            for locus in loci[2:end]
                locus ∉ x.loci && println("NOTICE: locus \"$locus\" not found in PopObj!")
            end
            println()

            return df[!, append!([:ind, :population], Symbol.(loci))]
        end
    else
        return df
    end
end


"""
    locations(x::PopObj)
View location data (`.longitude` and `.latitude`) in a `PopObj`

Use `locations!` to add spatial data to a `PopObj`
"""
function locations(x::PopObj)
    DataFrame(name = x.samples.name, population = x.samples.population, longitude = x.samples.longitude, latitude = x.samples.latitude)
end

"""
    locations!(x::PopObj; xloc::Array, yloc::Array)
Add location data (longitude `xloc`, latitude `yloc`) to `PopObj`. Takes decimal
degrees or decimal minutes format. **Must** use `-` symbol instead of cardinal directions.
Location data must be in order of `ind`. Replaces existing `PopObj` location data.
- Decimal Degrees : `-11.431`
- Decimal Minutes : `"-11 43.11"` (must use space and double-quotes)

If conversion is not necessary, can directly assign `PopObj.longitude` and `PopObj.latitude`
"""
function locations!(x::PopObj; xloc::Array, yloc::Array)
    # test for decimal degrees vs decimal minutes
    if occursin(" ", string(xloc[1])) == false && occursin(" ", string(yloc[1])) == false
        x.samples.longitude = xloc ;
        x.samples.latitude = yloc ;
    else
        # make sure decimal minutes are Strings
        if typeof(xloc) != Array{String,1}
            xloc = string.(xloc)
        end
        if typeof(yloc) != Array{String,1}
            yloc = string.(yloc)
        end
        # convert xloc to decimal degrees
        xlocConverted = []
        for value in xloc
            tmp = split(value, " ")
            if parse(Float64,tmp[1]) < 0   # if negative, subtract
                decideg = parse(Float64, tmp[1]) - round((parse(Float64,tmp[2])/60), digits = 3)
                push!(xlocConverted, decideg)
            else                           # if positive, add
                decideg = parse(Float64, tmp[1]) + round((parse(Float64,tmp[2])/60), digits = 3)
                push!(xlocConverted, decideg)
            end
        end
        # convert yloc to decimal degrees
        ylocConverted = []
        for value in yloc
            tmp = split(value, " ")
            if parse(Float64,tmp[1]) < 0
                decideg = parse(Float64, tmp[1]) - round((parse(Float64,tmp[2])/60), digits = 3)
                push!(ylocConverted, decideg)
            else
                decideg = parse(Float64, tmp[1]) + round((parse(Float64,tmp[2])/60), digits = 3)
                push!(ylocConverted, decideg)
            end
        end
        x.samples.longitude = xlocConverted
        x.samples.latitude = ylocConverted
        return x.samples
    end
end


"""
    populations(x::PopObj; listall::Bool = false)
View unique population ID's in a `PopObj`.

`listall = true`, displays `ind` and their `population` instead (default = `false`).
"""
function population(x::PopObj; listall::Bool = false)
    if listall == true
        DataFrame(name = x.samples.name, population = x.samples.population)
    else
        println( "   ", " #Inds | Pop " )
        println( "   ", "--------------" )
        popcounts = hcat([sum(x.samples.population .== i) for i in unique(x.samples.population)],unique(x.samples.population))
        for eachpop in 1:length(popcounts)÷2
            println("\t", popcounts[eachpop], "\t", " |", "\t", popcounts[eachpop,2])
        end
    end
end

"""
    populations!(x::PopObj; rename::Dict)
Rename the population ID's of `PopObj.population`.

Uses a `Dict` of `[population] => replacement` to rename populations

Example:

potatopops = Dict(1 => "Idaho", 2 => "Russet")

population!(potatoes, rename = potatopops)
"""
function population!(x::PopObj; rename::Dict)
    for eachkey in keys(rename)
        replace!(x.samples.population, eachkey => rename[eachkey])
    end
    population(x, listall = true)
end


"""
    genotypes(x::PopObj; inds::Union{String, Array, Nothing}= nothing)
Get all the genotypes of specific individuals within a `PopObj`.
- Names must be in quotes

Examples:

genotypes(eggplant, inds = ["ital_001", "ital_101", "spai_031"])

genotypes(eggplant, inds = "ital_001")
"""
function genotypes(x::PopObj; inds::Union{String, Array, Nothing}= nothing)
    df = deepcopy(x.loci)
    insertcols!(df, 1, :ind => x.samples.name) ;
    insertcols!(df, 2, :population => categorical(x.samples.population))
    if inds != nothing
        if typeof(inds) == String
            inds ∉ x.samples.name && error("individual $inds not found in PopObj")
            return df[df.ind .== inds, :]
        else
            tmp = df[df.ind .== inds[1], :]
            for ind in inds[2:end]
                ind ∉ x.samples.name && println("NOTICE: individual \"$ind\" not found in PopObj!")
                tmp = vcat(tmp, df[df.ind .== ind, :])
            end
            println()
            return tmp
        end
    else
        return df
    end
end


#### Find missing ####

"""
    Base.missing(x::PopObj)
Identify and count missing loci in each individual of a `PopObj`. Returns a tuple
of `DataFrames`: loci per individual, number per loci.

Example:

`aardvark = genepop("aardvark.gen", numpop = 5)`  # load file to PopObj

`missing_ind,missing_loc = missing(aardvark)`
"""
function Base.missing(x::PopObj)
    df = deepcopy(x.loci)
    insertcols!(df, 1, :ind => x.samples.name)
    # missing per individual
    nmissing = []
    missing_array = []
    for each in 1:length(df[:,1])
        miss_idx = findall(i -> i == (0,0), df[each,:])
        push!(nmissing, miss_idx |> length)
        push!(missing_array, String.(miss_idx))
    end
    ind_df = DataFrame(ind = x.samples.name,
                       population = x.samples.population,
                       nmissing = nmissing,
                       loci = missing_array
                       )
    # missing per locus
    d = Dict()
    for b in ind_df[!, :4]
        for c in b
            if c ∉ keys(d)
                d[c] = 1
            else
                d[c] += 1
            end
        end
    end

    # add loci without any missing
    countarray = []
    for locus in string.(names(x.loci))
        if locus ∉ keys(d)
            d[locus] = 0
        end
        push!(countarray, d[locus])
    end

    #convert to DF and sort
    loci_df = DataFrame(locus = string.(names(x.loci)), nmissing = countarray)
    return (ind_df, loci_df)
end

##### Removal #####

"""
    remove_loci!(x::PopObj, loci::Union{String, Array{String,1}})
Removes selected loci from a `PopObj`.

Examples:

`remove_loci!(tulips, "north_011")`

`remove_loci!(tulips, ["north_011", "north_003", "south_051"])`
"""
function remove_loci!(x::PopObj, loci::Union{String,Array{String,1}})
    # get individuals indices
    if typeof(loci) == String
        loci ∉ x.loci && error("Locus \"$loci\" not found")
        idx = findfirst(i -> i == loci, x.loci)
    else
        idx = []
        for each in loci
            if each ∉ x.loci
                println("NOTICE: locus \"$each\" not found")
                continue
            end
            push!(idx, findfirst(i -> i == each, x.loci) )
        end
        println()
    end
    deleteat!(x.loci, idx) # remove locus names
    for each in loci
        delete!(x.genotypes, each)  # remove genotypes
    end
    return x
end


"""
    remove_inds!(x::PopObj, inds::Union{Array{String,1}})
Removes selected individuals from a `PopObj`.

Examples:

`remove_inds!(sunflowers, "west_011")`

`remove_inds!(sunflowers, ["west_011", "west_003", "east_051"])`
"""
function remove_inds!(x::PopObj, inds::Union{String, Array{String,1}})
    # get inds indices
    if typeof(inds) == String
        inds ∉ x.samples.name && error("ind \"$inds\" not found")
        idx = findfirst(i -> i == inds, x.samples.name)
    else
        idx = []
        for ind in inds
            if ind ∉ x.samples.name
                println("NOTICE: ind \"$ind\" not found!")
                continue
            end
            push!(idx, findfirst(i -> i == ind, x.samples.name) )
        end
        println()
    end
    deleteat!(x.samples.name, idx)  # remove name(s)
    deleteat!(x.samples.population, idx)    # remove population(s)
    if length(x.samples.longitude) != 0 && length(x.samples.latitude) != 0
        deleteat!(x.samples.longitude, idx)    # remove xloc(s)
        deleteat!(x.samples.longitude, idx)    # remove yloc(s)
    end
    # remove inds from all loci genotypes
    for each in x.loci
        deleteat!(x.genotypes[each], idx)
    end
    return x
end
