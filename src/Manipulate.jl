#=
These are commands that are for the general manipulation and viewing of the
PopData type. The appear in alphabetical order.
=#
export locations, locations!, loci, locus, get_genotypes, get_sample_genotypes, missing, populations, population, populations!, population!, reindex, exclude_loci, remove_loci, omit_loci, exclude_samples, remove_samples, omit_samples, samples

"""
    locations(data::PopData)
View the longitude and latitude data in a `PopData` object. Returns a table
derived from the PopData. Changes made to this table will not alter the source
`PopData` object.

Use `locations!` to add spatial data to a `PopData` object.
"""
function locations(data::PopData)
    @view data.meta[!, [:longitude, :latitude]]
end


"""
    locations!(data::PopData; long::Vector{Float64}, lat::Vector{Float64})
Replaces existing `PopData` location data (longitude `long`, latitude `lat`).
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
function locations!(data::PopData, long::Vector{Union{Missing,T}}, lat::Vector{Union{Missing,T}}) where T <: AbstractFloat
    long_len = length(long)
    lat_len = length(lat)
    long_len != lat_len && error("latitude ($lat_len) and longitude ($long_len) arrays not equal in length")
    long_len != length(data.meta.name) && error("lat/long array length ($long_len) and number of samples in PopData ($long_len) are not equal")

    data.meta.longitude .= long
    data.meta.latitude .= lat
    return
end

function locations!(data::PopData, long::Vector{T}, lat::Vector{T}) where T <: AbstractFloat
    # convert to the right type and use locations!()
    lat_adjust = lat |> Vector{Union{Missing, Float32}}
    long_adjust = long |> Vector{Union{Missing, Float32}}
    locations!(data, lat = lat_adjust, long = long_adjust)
end


"""
    locations!(data::PopData; long::Vector{String}, lat::Vector{String})
Replaces existing `PopData` location data (longitude `long`, latitude `lat`). Takes
decimal minutes format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.
## Formatting requirements
- Decimal Minutes: `"-11 43.11"` (must use space and be a `String`)
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
### Example
```
ncats = nancycats();
x = fill("-11 22.33", 237) ; y = fill("-41 31.52", 237)
locations!(ncats, long = x, lat = y)
```
"""
function locations!(data::PopData, long::Vector{String}, lat::Vector{String})
    long_len = length(long)
    lat_len = length(lat)
    lat_len != long_len && error("latitude ($lat_len) and longitude ($long_len) arrays not equal in length")
    lat_len != length(data.meta.name) && error("lat/long array length ($lat_len) and number of samples in PopData ($long_len) are not equal")
    @info "Converting decimal minutes to decimal degrees"
    # convert coordinates to decimal degrees
    data.meta.longitude .= convert_coord.(long)
    data.meta.latitude .= convert_coord.(lat)
    return
end

function locations!(
    data::PopData,
    lat_deg::Vector{Union{Missing,Int}},
    lat_min::Vector{Union{Missing,Float64}},
    long_deg::Vector{Union{Missing,Int}},
    long_min::Vector{Union{Missing,Float64}}
)
    long_string = string.(lat_deg, " ", lat_min)
    lat_string = string.(long_deg, " ", long_min)
    locations!(data, long_string, lat_string)
end

function locations!(
    data::PopData,
    lat_deg::Vector{Int},
    lat_min::Vector{Float64},
    long_deg::Vector{Int},
    long_min::Vector{Float64}
)
    long_string = string.(lat_deg, " ", lat_min)
    lat_string = string.(long_deg, " ", long_min)
    locations!(data, long_string, lat_string)
end


"""
    locations!(data::PopData; kwargs...)
Replaces existing `PopData` location data (longitude, latitude). Requires all four
keyword arguments. Takes decimal minutes format as vectors of degrees (`Int`) and
decimal minutes (`Float`). Recommended to use `CSV.read` from `CSV.jl` to import
your spatial coordinates from a text file.
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
    levels(data.loci.locus)
end

# TODO all the manipulation functions below this
"""
    get_genotypes(data::PopObj; samples::Union{String, Vector{String}}, loci::Union{String, Vector{String}})
Return the genotype(s) of one or more `samples` for one or more
specific `loci` (both as keywords) in a `PopData` object.
### Examples
```
cats = nancycats();
get_genotype(cats, samples = "N115" , loci = "fca8")
get_genotypes(cats, samples = ["N1", "N2"] , loci = "fca8")
get_genotype(cats, samples = "N115" , loci = ["fca8", "fca37"])
get_genotype(cats, samples = ["N1", "N2"] , loci = ["fca8", "fca37"])
```

"""
function get_genotypes(data::PopData; sample::Union{String, Vector{String}}, locus::Union{String, Vector{String}})
    if typeof(sample) == String
        sample = [sample]
    end
    if typeof(locus) == String
        locus = [locus]
    end
    @where data.loci :name in sample && :locus in locus
end

"""
    get_sample_genotypes(data::PopData, sample::String)
Return all the genotypes of a specific sample in a `PopData` object.
This is an extension for the internal function `get_genotypes`.
```
cats = nancycats()
get_sample_genotypes(cats, "N115")
```
"""
function get_sample_genotypes(data::PopObj, sample::String)
    @where data.loci :name == sample
end


"""
    locus(data::PopData, locus::Union{String, Symbol})
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

        return @groupby data.loci :name {missing = count(ismissing,:genotype)}

    elseif mode == "pop" || mode == "population"

        return @groupby data.loci :population {missing = count(ismissing, :genotype)}

    elseif mode == "locus" || mode == "loci"

        return @groupby data.loci :locus {missing = count(ismissing, :genotype)}

    elseif mode == "detailed" || mode == "full"

        return @groupby data.loci (:locus, :population) {missing = count(ismissing, :genotype)}

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
    show(populations(data))
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
    show(populations(data))
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

"""
    reindex(data::PopData, col::Union{String, Symbol})
Re-index and sort the `loci` table of a `PopData` object by column
`col`. Returns a new `PopData` object.

### Example
sharks = gulfsharks()
reindex(sharks, :population)
"""
function JuliaDB.reindex(data::PopData, col::Union{String, Symbol})
    if typeof(col) == String
        col = Symbol(col)
    end
    return PopData(data.meta, reindex(data.loci, col))
end

##### Exclusion #####
#=
To have a built-in "undo button", exclusion functions return new PopData objects
with the specified loci/samples removed rather than overwriting the original.
=#
"""
```
exclude_loci(data::PopData, locus::String)
exclude_loci(data::PopData, loci::Vector{String})
```
Exclude selected loci from a `PopData` object. Returns a new `PopData` object,
leaving the original intact. Synonymous with `omit_loci` and `remove_loci`.

### Examples
```
new_cats = exclude_loci(nancycats(), "fca8")
very_new_cats = exclude_loci(nancycats(), ["fca8", "fca23"])
```
"""
function exclude_loci(data::PopData, locus::String)
    locus ∉ loci(data) && error("Locus \"$locus\" not found")
    new_table = @where data.loci :locus != locus
    return PopData(data.meta, new_table)
end

function exclude_loci(data::PopData, exloci::Vector{String})
    msg = ""
    all_loci = loci(data)
    for each in exloci
        if each ∉ all_loci
            msg *= "\n  locus \"$each\" not found"
        end
    end
    new_table =  @where data.loci :locus ∉ exloci
    msg != "" && printstyled("Warnings:", color = :yellow) ; println(msg)
    return PopData(data.meta, new_table)
end

const omit_loci = exclude_loci
const remove_loci = exclude_loci

"""
    exclude_samples(data::PopData, samp_id::Union{Vector{String}})
Exclude selected samples from a `PopData` object. Returns a new `PopData` object,
leaving the original intact. Synonymous with `omit_samples` and `remove_samples`.

### Examples
```
exclude_samples(nancycats, "N100")
exclude_samples(nancycats, ["N100", "N102", "N211"])
```
"""
function exclude_samples(data::PopData, samp_id::String)
    samp_id ∉ samples(data) && error("Sample \"$samp_id\" not found")
    new_loc_table = @where data.loci :name != samp_id
    new_meta_table = @where data.meta :name != samp_id
    return PopData(new_meta_table, new_loc_table)
end


function exclude_samples(data::PopData, samp_ids::Vector{String})
    msg = ""
    all_samp = samples(data)
    for each in samp_ids
        if each ∉ all_samp
            msg *= "\n  Sample \"$each\" not found"
        end
    end
    new_loc_table = @where data.loci :name ∉ samp_ids
    new_meta_table = @where data.meta :name ∉ samp_ids
    msg != "" && printstyled("Warnings:", color = :yellow) ; println(msg)
    return PopData(new_meta_table, new_loc_table)
end

const omit_samples = exclude_samples
const remove_samples = exclude_samples

"""
    samples(data::PopData)
View individual/sample names in a `PopData`
"""
function samples(data::PopData)
    select(data.meta, :name)
end
