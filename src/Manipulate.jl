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
function Base.summary(x::PopObj)
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
    println("$(x.samples.ploidy[1:3])", " \u2026 ", "$(x.samples.ploidy[end-2:end])")
    println("Number of populations: $(length(x.samples.population |> unique))","\n")
    println("#samp_id | Pop","\n", "--------------")
    popcounts = hcat([sum(x.samples.population .== i) for i in unique(x.samples.population)],unique(x.samples.population))
    for eachpop in 1:length(popcounts)÷2
        println(popcounts[eachpop], " | ", popcounts[eachpop,2])
    end
    println("\nAvailable .samples fields: .name, .population, .ploidy, .longitude, .latitude")
end


"""
    isolate_genotypes(x::PopObj; samples::Union{String, Array, Nothing}, loci::Union{String, Array, Nothing})
View the genotypes of specific samples for specific loci in a `PopObj`.
Default shows all genotypes for all individuals.

`isolate_genotypes(nancycats, loci = "fca8")`

`isolate_genotypes(nancycats, samples = "N226", loci = ["fca8", "fca23"])`
"""
function isolate_genotypes(x::PopObj; samples::Union{String, Array, Nothing}= nothing, loci::Union{String, Array, Nothing}= nothing)
    if loci == nothing && samples == nothing
        @warn "please specify either loci= or samples=, otherwise use PopObj.loci"
    end
    df = deepcopy(x.loci)
    insertcols!(df, 1, :name => x.samples.name) ;
    insertcols!(df, 2, :population => categorical(x.samples.population))
    if samples != nothing
        if typeof(samples) == String
            samples ∉ x.samples.name && error("individual $samples not found in PopObj")
            tmp = df[df.name .== samples, :]
        else
            tmp = df[df.name .== samples[1], :]
            for ind in samples[2:end]
                ind ∉ x.samples.name && println("NOTICE: individual \"$ind\" not found in PopObj!")
                tmp = vcat(tmp, df[df.name .== ind, :])
            end
            println()
        end
    else
        tmp = df
    end
    if loci != nothing
        if typeof(loci) == String
            loci ∉ string.(names(x.loci)) && error("locus $loci not found in PopObj")
            return tmp[!, [:name, :population, Symbol(loci)]]
        else
            for locus in loci[2:end]
                locus ∉ string.(names(x.loci)) && println("NOTICE: locus \"$locus\" not found in PopObj!")
            end
            println()
            return tmp[!, append!([:name, :population], Symbol.(loci))]
        end
    else
        return tmp
    end
end


"""
    locations(x::PopObj)
View location data (`.longitude` and `.latitude`) in a `PopObj`

Use `locations!` to add spatial data to a `PopObj`
"""
function locations(x::PopObj)

    DataFrame(name = x.samples.name,
              population = x.samples.population,
              longitude = x.samples.longitude,
              latitude = x.samples.latitude)
end

"""
    locations!(x::PopObj; lat::Array, long::Array)
Add location data (latitude `lat`, longitude `long`) to `PopObj`. Takes decimal
degrees or decimal minutes format. **Must** use `-` symbol instead of cardinal directions.
Location data must be in order of `ind`. Replaces existing `PopObj` location data.
- Decimal Degrees : `-11.431`
- Decimal Minutes : `"-11 43.11"` (must use space and double-quotes)

If conversion is not necessary, can directly assign `PopObj.samples.longitude` and `PopObj.samples.latitude`
"""
function locations!(x::PopObj; lat::Array, long::Array)
    # test for decimal degrees vs decimal minutes
    if occursin(" ", string(lat[1])) == false && occursin(" ", string(long[1])) == false
        x.samples.longitude = lat ;
        x.samples.latitude = long ;
    else
        # make sure decimal minutes are Strings
        if typeof(lat) != Array{String,1}
            lat = string.(lat)
        end
        if typeof(long) != Array{String,1}
            long = string.(long)
        end
        # convert lat to decimal degrees
        latConverted = []
        for value in lat
            tmp = split(value, " ")
            if parse(Float64,tmp[1]) < 0   # if negative, subtract
                decideg = parse(Float64, tmp[1]) - round((parse(Float64,tmp[2])/60), digits = 3)
                push!(latConverted, decideg)
            else                           # if positive, add
                decideg = parse(Float64, tmp[1]) + round((parse(Float64,tmp[2])/60), digits = 3)
                push!(latConverted, decideg)
            end
        end
        # convert long to decimal degrees
        longConverted = []
        for value in long
            tmp = split(value, " ")
            if parse(Float64,tmp[1]) < 0
                decideg = parse(Float64, tmp[1]) - round((parse(Float64,tmp[2])/60), digits = 3)
                push!(longConverted, decideg)
            else
                decideg = parse(Float64, tmp[1]) + round((parse(Float64,tmp[2])/60), digits = 3)
                push!(longConverted, decideg)
            end
        end
        x.samples.longitude = latConverted
        x.samples.latitude = longConverted
        return x.samples

    end
end

"""
    population(x::PopObj; listall::Bool = false)
View unique population ID's in a `PopObj`.

`listall = true`, displays `ind` and their `population` instead (default = `false`).
"""
function population(x::PopObj; listall::Bool = false)
    if listall == true
        DataFrame(name = x.samples.name, population = x.samples.population)
    else
        count = [sum(x.samples.population .== i) for i in unique(x.samples.population)]
        count_conv = Int32.(count)
        popcounts = DataFrame(population = unique(x.samples.population) |> categorical,
                              count = count_conv)
    end
end

"""
    population!(x::PopObj; rename::Dict)
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
    populations(x::PopObj; listall::Bool = false)
View unique population ID's in a `PopObj`.

`listall = true`, displays `ind` and their `population` instead (default = `false`).
"""
function populations(x::PopObj; listall::Bool = false)
    if listall == true
        DataFrame(name = x.samples.name, population = x.samples.population)
    else
        count = [sum(x.samples.population .== i) for i in unique(x.samples.population)]
        count_conv = Int32.(count)
        popcounts = DataFrame(population = unique(x.samples.population) |> categorical,
                              count = count_conv)
    end
end

"""
    populations!(x::PopObj; rename::Dict)
Rename the population ID's of `PopObj.population`.

Uses a `Dict` of `[population] => replacement` to rename populations

Example:

potatopops = Dict(1 => "Idaho", 2 => "Russet")

populations!(potatoes, rename = potatopops)
"""
function populations!(x::PopObj; rename::Dict)
    for eachkey in keys(rename)
        replace!(x.samples.population, eachkey => rename[eachkey])
    end
    population(x, listall = true)
end

#### Find missing ####

"""
    Base.missing(x::PopObj)
Identify and count missing loci in each sample of a `PopObj`. Returns a tuple
of `DataFrames`: loci per sample, number per loci.

Example:

`missing_ind,missing_loc = missing(gulfsharks)`
"""
function Base.missing(x::PopObj)
    df = deepcopy(x.loci)
    insertcols!(df, 1, :ind => x.samples.name)
    # missing per individual
    nmissing = []
    missing_array = []
    for each in 1:length(df[:,1])
        miss_idx = findall(i -> i === missing, df[each,:])
        push!(nmissing, miss_idx |> length)
        push!(missing_array, String.(miss_idx))
    end
    sample_df = DataFrame(name = x.samples.name,
                       population = x.samples.population,
                       nmissing = nmissing,
                       loci = missing_array
                       )
    # missing per locus
    f(x) = map(eachcol(x)) do col
        count(i->i===missing, col)
    end

    loci_df = DataFrame(locus = string.(names(x.loci)), nmissing = f(x.loci))
    return (by_sample = sample_df, by_loci = loci_df)
end

##### Removal #####

"""
    remove_loci!(x::PopObj, loci::Union{String, Array{String,1}})
Removes selected loci from a `PopObj`.

Examples:

`remove_loci!(nancycats, "fca8")`

`remove_loci!(tulips, ["fca8", "fca23"])`
"""
function remove_loci!(x::PopObj, loci::Union{String,Array{String,1}})
    # get individuals indices
    if typeof(loci) == String
        sym_loci = Symbol(loci)
        sym_loci ∉ names(x.loci) && error("Locus \"$loci\" not found")
        return select!(x.loci, Not(sym_loci))
    else
        sym_loci = Symbol.(loci)
        for each in sym_loci
            if each ∉ names(x.loci)
                println("NOTICE: locus \"$each\" not found")
                continue
            end
        end
        return select!(x.loci, Not(sym_loci))
    end
end


"""
    remove_samples!(x::PopObj, samp_id::Union{Array{String,1}})
Removes selected samples from a `PopObj`.

Examples:

`remove_samples!(nancycats, "N100")`

`remove_samples!(nancycats, ["N100", "N102", "N211"])`
"""
function remove_samples!(x::PopObj, samp_id::Union{String, Array{String,1}})
    # get samp_id indices
    if typeof(samp_id) == String
        samp_id ∉ x.samples.name && error("sample \"$samp_id\" not found")
        idx = findfirst(i -> i == samp_id, x.samples.name)
    else
        idx = []
        for ind in samp_id
            if ind ∉ x.samples.name
                println("NOTICE: sample \"$ind\" not found!")
                continue
            end
            push!(idx, findfirst(i -> i == ind, x.samples.name) )
        end
        println()
    end
    deleterows!(x.samples, idx)
    deleterows!(x.loci, idx)
    return x
end
