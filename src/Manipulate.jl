#=
These are commands that are for the general manipulation and viewing of the
PopData type. The appear in alphabetical order.
=#

"""
    locations(data::PopData)
View the longitude and latitude data in a `PopData` object. Returns a table
derived from the PopData. Changes made to this table will not alter the source
`PopData` object.

Use `locations!` to add spatial data to a `PopData` object.
"""
function locations(data::PopData)
    select(data.meta, (:longitude, :latitude))
end


"""
    locations!(data::PopData; lat::Vector{Float64}, long::Vector{Float64})
Replaces existing `PopData` location data (latitude `lat`, longitude `long`).
Takes decimal degrees as a `Vector` of any `AbstractFloat`.
## Formatting requirements
- Decimal Degrees format: `-11.431`
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)
### Example
```
ncats = nancycats() ;
x = rand(237) ; y = rand(237)
locations!(ncats, long = x, lat = y)
```
"""
function locations!(data::PopData, lat::Vector{Union{Missing,T}}, long::Vector{Union{Missing,T}}) where T <: AbstractFloat
    length(lat) != length(long) && error("latitude ($(length(lat))) and longitude ($(length(long))) arrays not equal in length")
    length(lat) != length(data.meta.columns.name) && error("lat/long array length ($(length(lat))) and number of samples in PopData $(length(long)) are not equal")

    for i in 1:length(lat)
        data.meta.columns.longitude[i] = lat[i]
        data.meta.columns.latitude[i] = long[i]
    end
end

function locations!(data::PopData, lat::Vector{T}, long::Vector{T}) where T <: AbstractFloat
    # convert to the right type and use locations!()
    lat_adjust = lat |> Vector{Union{Missing, Float32}}
    long_adjust = long |> Vector{Union{Missing, Float32}}
    locations!(data, lat = lat_adjust, long = long_adjust)
end


"""
    locations!(data::PopData; lat::Vector{String}, long::Vector{String})
Replaces existing `PopData` location data (latitude `lat`, longitude `long`). Takes
decimal minutes format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.
## Formatting requirements
- Decimal Minutes: `"-11 43.11"` (must use space and be a `String`)
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)
### Example
```
ncats = nancycats();
x = fill("-11 22.33", 237) ; y = fill("-41 31.52", 237)
locations!(ncats, long = x, lat = y)
```
"""
function locations!(data::PopData, lat::Vector{Union{Missing,String}}, long::Vector{Union{Missing,String}})
    length(lat) != length(long) && error("latitude ($(length(lat))) and longitude $(length(long)) arrays not equal in length")
    length(lat) != length(data.meta.columns.name) && error("lat/long array length ($(length(lat))) and number of samples in PopData $(length(long)) are not equal")
    @info "Converting decimal minutes to decimal degrees"
    # convert coordinates to decimal degrees
    latConverted = Vector{Union{Missing,Float32}}()
    longConverted = Vector{Union{Missing,Float32}}()
    @inbounds for idx in 1:length(lat)
        tmpLat = split(lat[idx], " ")
        tmpLong = split(long[idx], " ")
        lat_deg = parse(Float32,tmpLat[1]) ; lat_min = round(parse(Float32,tmpLat[2])/60, digits = 4)
        long_deg = parse(Float32,tmpLong[1]) ; long_min = round(parse(Float32,tmpLong[2])/60, digits = 4)
        if lat_deg < 0
            # if negative, subtract
            push!(latConverted, lat_deg - lat_min)
        else
            # if positive, add
            push!(latConverted, lat_deg + lat_min)
        end
        if long_deg < 0
            # if negative, subtract
            push!(longConverted, long_deg - long_min)
        else
            # if positive, add
            push!(longConverted, long_deg + long_min)
        end
    end
    for i in 1:length(latConverted)
        data.meta.columns.longitude[i] = latConverted[i]
        data.meta.columns.latitude[i] = longConverted[i]
    end
end


function locations!(
    data::PopData,
    lat_deg::Vector{Union{Missing,Int}},
    lat_min::Vector{Union{Missing,Float64}},
    long_deg::Vector{Union{Missing,Int}},
    long_min::Vector{Union{Missing,Float64}},
)
    if any(length(lat_deg) .!= length.([long_deg, long_min, lat_min]))
        error("input array lengths do not match each other\n  lat_deg: $(length(lat_deg))\n  lat_min: $(length(lat_min))\n  long_deg: $(length(long_deg))\n  long_min: $(length(long_min))")
    elseif length(lat_deg) != length(data.meta.columns.name)
        error("lat/long array length ($(length(lat))) and number of samples in PopData $(length(long)) are not equal")
    end

    @info "Converting decimal minutes to decimal degrees"
    # convert coordinates to decimal degrees
    latConverted = Vector{Union{Missing,Float32}}()
    longConverted = Vector{Union{Missing,Float32}}()
    @inbounds for i = 1:length(lat_deg)
        if lat_deg[i] === missing
            push!(latConverted, missing)
            push!(longConverted, missing)
            continue
        elseif lat_deg[i] < 0
            # if negative, subtract
            push!(latConverted, lat_deg[i] - lat_min[i])
        else
            # if positive, add
            push!(latConverted, lat_deg[i] + lat_min[i])
        end
        if long_deg < 0
            # if negative, subtract
            push!(longConverted, long_deg[i] - long_min[i])
        else
            # if positive, add
            push!(longConverted, long_deg[i] + long_min[i])
        end
    end

    @inbounds for i = 1:length(latConverted)
        data.meta.columns.longitude[i] = latConverted[i]
        data.meta.columns.latitude[i] = longConverted[i]
    end
end


"""
    locations!(data::PopData; kwargs...)
Replaces existing `PopData` location data (latitude, longitude). Requires all four
keyword arguments. Takes decimal minutes format as vectors of degrees (`Int`) and
decimal minutes (`Float`). Recommended
to use `CSV.read` from `CSV.jl` to import your spatial coordinates from a text file.
### Keyword Arguments:
- `lat_deg::Vector{Int}` a vector of postive or negative integers denoting the latitude degrees
    - example: `[11, -12, 15, 11]`
- `lat_min::Vector{Float64}` a vector of positive floating point numbers denoting the latitude decimal minutes
    - example: `[43.12, 41.32, 36.53, 22.41]`
- `long_deg::Vector{Int}` same as `lat_deg` but for longitude degrees
- `long_min::Vector{Float64}` same as `lat_min` but for longitude minutes

## Formatting requirements
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)

### Example
If you have decimal-minutes coordinates for two samples:
- Sample 1  _Long:_ 11 43.12  _Lat:_ 15 36.53
- Sample 2  _Long:_ -12 41.32 _Lat:_ 11 22.41
\n then your inputs would be:
```
lo_deg = [11, -12]
lo_min = [43.12, 41.32]
la_deg  = [15, 11]
la_min  = [36.53, 22.41]
locations!(data, long_deg = lo_deg, long_min = lo_min, lat_deg = la_deg, lat_min = la_min)
```
"""
function locations!(data::PopData; kwargs...)
    kwargs = Dict(kwargs)
    # vectors as strings mode
    if any([haskey(kwargs, :lat), haskey(kwargs, :long)])
        # check for matching lat and long keywords
        if all([haskey(kwargs, :lat), haskey(kwargs, :long)])
            locations!(data, kwargs[:lat], kwargs[:long])
        else
            error("keyword arguments \"lat\" and \"long\" must be supplied together")
        end
    end

    # split-numeric vector mode
    if any([haskey(kwargs, :lat_deg), haskey(kwargs, :lat_min), haskey(kwargs, :long_deg), haskey(kwargs, :long_min)])
        # check for matching lat/long deg/min
        if all([haskey(kwargs, :lat_deg), haskey(kwargs, :lat_min), haskey(kwargs, :long_deg), haskey(kwargs, :long_min)])
            locations!(data, kwargs[:lat_deg], kwargs[:lat_min], kwargs[:long_deg], kwargs[:long_min])
        else
            error("keyword arguments \"lat_deg\", \"lat_min\", \"long_deg\", and \"long_min\" must be supplied together")
        end
    end
end


"""
    loci(data::PopData)
Returns an array of strings of the loci names in a `PopData` object.
"""
function loci(data::PopData)
    levels(data.loci.columns.locus)
end

"""
    loci(data::IndexedTable)
Convenience wrapper to return an array of column names as string in the `loci`
Table of a `PopData` object.
"""
function loci(data::IndexedTable)
    levels(data.columns.locus)
end

"""
    locus(::PopData, ::Union{String, Symbol})
Convenience wrapper to return a vector of all the genotypes of a single locus

### Example
```
locus(gulfsharks(), "contig_475")
```
"""
function locus(data::PopData, locus::String)
    tmp = select(data.loci, (:locus, :genotype)) |>
        @where :locus == locus
    return select(tmp, :genotype)
end

function meta(data::PopData)
    data.meta
end

const metadata = meta

#### Find missing ####
"""
    missing(data::PopData; mode::String = "sample")
Get missing genotype information in a `PopData`. Specify a mode of operation
to return a DataFrame corresponding with that missing information.

#### Modes
- "sample" - returns a count and list of missing loci per individual (default)
- "pop" - returns a count of missing genotypes per population
- "locus" - returns a count of missing genotypes per locus
- "full" - returns a count of missing genotypes per locus per population

### Example:
```
missing(gulfsharks(), mode = "pop")
```
"""
@inline function Base.missing(data::PopData; mode::String = "sample")
    if mode == "sample" || mode == "individual"

        return @groupby data.loci :name {missing = sum(ismissing.(:genotype))}

    elseif mode == "pop" || mode == "population"

        return @groupby data.loci :population {missing = sum(ismissing.(:genotype))}

    elseif mode == "locus" || mode == "loci"

        return @groupby data.loci :locus {missing = sum(ismissing.(:genotype))}

    elseif mode == "detailed" || mode == "full"

        return @groupby data.loci (:locus, :population) {missing = sum(ismissing.(:genotype))}

    else
        @error "Mode \"$mode\" not recognized. Please specify one of: sample, pop, locus, or full"
        missing(data)
    end
end

"""
    populations(data::PopData; listall::Bool = false)
View unique population ID's and their counts in a `PopData`.

- `listall = true` displays all samples and their `population` instead (default = `false`)
"""
@inline function populations(data::PopData; listall::Bool = false)
    if listall == true
        return select(data.meta, (:name, :population))
    else
        if all(ismissing.(data.meta.columns.population)) == true
            @info "no population data present in PopData"
            return populations(data, listall = true)
        end
        return @groupby data.meta :population {count = length(:population)}
    end
end


"""
```
# Replace by matching
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, oldnames::Vector{String}, newnames::Vector{String})

# Generate new population information
populations!(data::PopData, rename::Vector{String}, counts::Vector{T}) where T<:Signed
populations!(data::PopData, rename::NamedTuple)
```
Multiple methods to assign or reassign population names to a `PopData`.

## Rename using a Dictionary
Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`
\n**Example**
```
potatopops = Dict("1" => "Idaho", "2" => "Russet")
populations!(potatoes, potatopops)
```

## Rename using a Vector of Strings
`Vector` of new unique population names in the order that they appear in the PopData.meta
\n**Example**
```
potatopops = ["Idaho", "Russet"]
populations!(potatoes, potatopops)
```
## Rename using two Vectors of Strings
Similar to the `Dict` method, except instead of creating a dictionary of "oldname" => "newname"
you input a Vector{String} of `oldnames` followed by another of `newnames`. Logically, the
new names will replace the old names in the same order as they appear in PopData.meta(e.g. the first newname replaces the first oldname, etc.).

\n**Example**
```
populations!(potatoes, ["russet1", "russet2"], ["north_russet", "south_russet"])
```
## Replace using a NamedTuple
Generate new population names for a `PopData`, overwriting everthing/anything currently
there. Will generate an array of population names from a NamedTuple of
`(names = , counts = )` where `names` is an array of the names of the populations and
`counts` is an array of the number of samples per population.
\n**Example**
\nTo assign names for three populations in a `PopData` named "Starlings" where new
population names are "North", "South", "East" and their sizes are 15, 32, 11:
```
populations!(Starlings, (names = ["North","South", "East"], counts = [15,32,11]))
```

## Replace using a Vector of Strings and Vector of Integers
Just like the NamedTuple method, except without the NamedTuple. Use an Array of Strings
as the second argument to denote population names, and an Array of Integers as the third
argument to denote the number of samples per population.
\n**Example**
```
populations!(Starlings, ["North","South", "East"], [15,32,11])
```
"""
function populations!(data::PopData, rename::Dict)
    msg = ""
    @inbounds for key in keys(rename)
        if key ∉ levels(data.loci.columns.population)
            msg *= "  Population \"$key\" not found in PopData\n"
        else
            replace!(data.meta.columns.population, key => rename[key])
        end
    end
    msg != "" && printstyled("Warnings:", color = :yellow) ; print("\n"*msg)
    recode!(data.loci.columns.population,rename...)
    return
end

function populations!(data::PopData, rename::Vector{String})
    current_popnames = unique(data.meta.columns.population)
    rn_dict = Dict()
    [rn_dict[i] = j for (i,j) in zip(current_popnames, rename)]
    populations!(data, rn_dict)
    return
end


function populations!(data::PopData, rename::NamedTuple)
    tidy_len = length(levels(data.loci.columns.locus))
    popid_array = [fill(i, j) for (i,j) in zip(rename.names, rename.counts)] |>
        Iterators.flatten |> collect
    loci_array = [fill(i, j) for (i,j) in zip(rename.names, rename.counts .* tidy_len)] |>
        Iterators.flatten |> collect
    length(popid_array) != length(data.meta.columns.name) && error("length of names ($(length(popid_array))) does not match sample number ($(length(data.meta.columns.name)))")
    for i in 1:length(data.meta.columns.population)
        data.meta.columns.population[i] = popid_array[i]
    end

    for i in 1:length(data.loci.columns.population)
        data.loci.columns.population[i] = loci_array[i]
    end

    droplevels!(data.loci.columns.population)
    return
end

function populations!(data::PopData, rename::Vector{String}, counts::Vector{T}) where T<:Signed
    populations!(data, (names = rename, counts = counts))
    return
end

function populations!(data::PopData, oldnames::Vector{String}, newnames::Vector{String})
    length(newnames) != length(oldnames) && error("number of current names and new names must match \n\t current: $(length(oldnames))  new: $(length(newnames))")
    rn_dict = Dict()
    [rn_dict[i] = j for (i,j) in zip(oldnames, newnames)]
    populations!(data, rn_dict)
    return
end

const population = populations
const population! = populations!
const popnames! = populations!

##### Exclusion #####
#=
To have a built-in "undo button", exclusion functions return new PopData objects
with the specified loci/samples removed rather than overwriting the original.
=#
#TODO
"""
```
exclude_loci(data::PopData, locus::String)
exclude_loci(data::PopData, loci::Vector{String})
```
Exclude selected loci from a `PopData` object. Returns a new `PopData` object,
leaving the original intact. Synonymous with `omit_loci` and `remove_loci`.

### Examples
```
exclude_loci(nancycats(), "fca8")
exclude_loci!(nancycats(), ["fca8", "fca23"])
```
"""
function exclude_loci(data::PopData, locus::String)
    locus ∉ loci(data) && error("Locus \"$locus\" not found")
    new_table = @where data.loci :locus != locus
    droplevels!(new_table.columns.locus)
    return PopData(data.meta, new_table)
end

function exclude_loci(data::PopData, loci::Vector{String})
    msg = ""
    all_loc = loci(data)
    for each in loci
        if each ∉ all_loci
            msg *= "\n  locus \"$each\" not found"
        end
    end
    new_table =  @where data.loci :locus ∉ loci
    msg != "" && printstyled("Warnings:", color = :yellow) ; print("\n"*msg)
    return PopData(data.meta, new_table)
end

const omit_loci = exclude_loci
const remove_loci = exclude_loci

"""
    exclude_samples!(data::PopData, samp_id::Union{Vector{String}})
Exclude selected samples from a `PopData` object. Returns a new `PopData` object,
leaving the original intact. Synonymous with `omit_samples` and `remove_samples`.

### Examples
```
exclude_samples(nancycats, "N100")
exclude_samples(nancycats, ["N100", "N102", "N211"])
```
"""
function exclude_samples(data::PopData, samp_id::Union{String, Vector{String}})
    # get samp_id indices
    if typeof(samp_id) == String
        samp_id ∉ data.samples.name && error("sample \"$samp_id\" not found")
        idx = findfirst(i -> i == samp_id, data.samples.name)
    else
        idx = Vector{Int}()
        for ind in samp_id
            if ind ∉ data.samples.name
                println("NOTICE: sample \"$ind\" not found!")
                continue
            end
            push!(idx, findfirst(i -> i == ind, data.samples.name))
        end
        println()
    end
    deleterows!(data.samples, idx)
    deleterows!(data.loci, idx)
    return data
end

"""
    samples(data::PopData)
View individual/sample names in a `PopData`
"""
function samples(data::PopData)
    select(data.meta, :name)
end


"""
    view_genotypes(data::PopData; samples::Union{String, Array, Nothing}, loci::Union{String, Array, Nothing})
Returns a dataframe of samples, population, genotypes. View the genotypes of
specific samples for specific loci in a `PopData`. Default shows all genotypes
for all individuals.

### Examples
```
view_genotypes(nancycats, loci = "fca8")
view_genotypes(nancycats, samples = "N226", loci = ["fca8", "fca23"])
```
"""
function view_genotypes(data::PopData; samples::Union{String, Array, Nothing}= nothing, loci::Union{String, Array, Nothing}= nothing)
    if loci == nothing && samples == nothing
        @warn "please specify either loci= or samples=, otherwise use PopData.loci"
    end

    df = deepcopy(data.loci)
    insertcols!(df, 1, :name => data.samples.name)
    insertcols!(df, 2, :population => data.samples.population)
    if samples != nothing
        if typeof(samples) == String
            samples ∉ data.samples.name && error("individual $samples not found in PopData")
            tmp = df[df.name .== samples, :]
        else
            tmp = df[df.name .== samples[1], :]
            for ind in samples[2:end]
                ind ∉ data.samples.name && println("NOTICE: individual \"$ind\" not found in PopData!")
                tmp = vcat(tmp, df[df.name .== ind, :])
            end
            println()
        end
    else
        tmp = df
    end
    if loci != nothing
        if typeof(loci) == String
            loci ∉ string.(names(data.loci)) && error("locus $loci not found in PopData")
            return tmp[!, [:name, :population, Symbol(loci)]]
        else
            for locus in loci[2:end]
                locus ∉ string.(names(data.loci)) && println("NOTICE: locus \"$locus\" not found in PopData!")
            end
            println()
            return tmp[!, append!([:name, :population], Symbol.(loci))]
        end
    else
        return tmp
    end
end
