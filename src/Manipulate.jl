"""
    indnames(x::PopObj)
View individual/sample names in a `PopObj`

Equivalent to `PopObj.inds`
"""
indnames(x::PopObj) = x.ind


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
    if length(x.longitude) == 0 && length(x.latitude) == 0
        @info "location data not provided"
    elseif length(x.longitude) != length(x.latitude)
        @warn "dimensions of longitude and latitude data not equal"
        println("\t\tLengths:")
        println("longitude: $(length(x.longitude)) | latitude: $(length(x.latitude))")
    else
        DataFrame(ind = x.ind, population = x.popid, longitude = x.longitude, latitude = x.latitude)
    end
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
        x.longitude = xloc ;
        x.latitude = yloc ;
        println("  xloc    yloc")
        hcat(x.xloc, x.yloc)
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
        x.longitude = xlocConverted
        x.latitude = ylocConverted
        println("  xloc    yloc")
        hcat(x.xloc, x.yloc)
    end
end


"""
    popid(x::PopObj; listall::Bool = false)
View unique population ID's in a `PopObj`.

`listall = true`, displays `ind` and their `popid` instead (default = `false`).
"""
function popid(x::PopObj; listall::Bool = false)
    if listall == true
        DataFrame(ind = x.ind, population = x.popid)
    else
        println( "   ", " #Inds | Pop " )
        println( "   ", "--------------" )
        popcounts = hcat([sum(x.popid .== i) for i in unique(x.popid)],unique(x.popid))
        for eachpop in 1:length(popcounts)÷2
            println("\t", popcounts[eachpop], "\t", " |", "\t", popcounts[eachpop,2])
        end
    end
end

"""
    popid!(x::PopObj; rename::Dict)
Rename the population ID's of `PopObj.popid`.

Uses a `Dict` of `[popid] => replacement` to rename

Example:

potatopops = Dict(1 => "Idaho", 2 => "Russet")

popid!(potatoes, rename = potatopops)
"""
function popid!(x::PopObj; rename::Dict)
    for eachkey in keys(rename)
        replace!(x.popid, eachkey => rename[eachkey])
    end
    popid(x, listall = true)
end


"""
    genotypes(x::PopObj; inds::Array{String,1})
Get all the genotypes of specific individuals within a `PopObj`.
- Names must be in quotes

Examples:

genotypes(eggplant, inds = ["ital_001", "ital_101", "spai_031"])

genotypes(eggplant, inds = "ital_001")
"""
function genotypes(x::PopObj; inds::Union{String, Array, Nothing}= nothing)
    df = x.genotypes |> DataFrame ;
    insertcols!(df, 1, :ind => x.ind) ;
    insertcols!(df, 2, :population => categorical(x.popid))
    if inds != nothing
        if typeof(inds) == String
            inds ∉ x.ind && error("individual $inds not found in PopObj")
            return df[df.ind .== inds, :]
        else
            tmp = df[df.ind .== inds[1], :]
            for ind in inds[2:end]
                ind ∉ x.ind && println("NOTICE: individual \"$ind\" not found in PopObj!")
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
    Base.missing(x::PopObj; plot::Bool = false)
Identify and count missing loci in each individual of a `PopObj`. Returns a tuple
of `DataFrames`: loci per individual, number per loci.

Example:

`aardvark = genepop("aardvark.gen", numpop = 5)`  # load file to PopObj

`missing_ind,missing_loc = missing(aardvark)`
"""
function Base.missing(x::PopObj)
    df = x.genotypes |> DataFrame
    insertcols!(df, 1, :ind => x.ind)
    # missing per individual
    nmissing = []
    missing_array = []
    for each in 1:length(df[:,1])
        miss_idx = findall(i -> i == (0,0), df[each,:])
        push!(nmissing, miss_idx |> length)
        push!(missing_array, String.(miss_idx))
    end
    ind_df = DataFrame(ind = x.ind,
                       population = string.(x.popid),
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
    for locus in x.loci
        if locus ∉ keys(d)
            d[locus] = 0
        end
        push!(countarray, d[locus])
    end

    #convert to DF and sort
    loci_df = DataFrame(locus = x.loci, nmissing = countarray)
    return (ind_df, loci_df)
end

##### Removal #####

"""
    remove_loci!(x::PopObj; loci::Union{String, Array{String,1}})
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

`remove_ind!(sunflowers, "west_011")`

`remove_ind!(sunflowers, ["west_011", "west_003", "east_051"])`
"""
function remove_inds!(x::PopObj, inds::Union{String, Array{String,1}})
    # get inds indices
    if typeof(inds) == String
        inds ∉ x.ind && error("ind \"$inds\" not found")
        idx = findfirst(i -> i == inds, x.ind)
    else
        idx = []
        for ind in inds
            if ind ∉ x.ind
                println("NOTICE: ind \"$ind\" not found!")
                continue
            end
            push!(idx, findfirst(i -> i == ind, x.ind) )
        end
        println()
    end
    deleteat!(x.ind, idx)  # remove name(s)
    deleteat!(x.popid, idx)    # remove popid(s)
    if length(x.longitude) != 0 && length(x.latitude) != 0
        deleteat!(x.longitude, idx)    # remove xloc(s)
        deleteat!(x.longitude, idx)    # remove yloc(s)
    end
    # remove inds from all loci genotypes
    for each in x.loci
        deleteat!(x.genotypes[each], idx)
    end
    return x
end
