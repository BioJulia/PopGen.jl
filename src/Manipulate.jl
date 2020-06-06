#=
These are commands that are for the general manipulation and viewing of the
PopData type. The appear in alphabetical order.
=#
export locations, locations!, loci, locus, get_genotypes, get_sample_genotypes, missing, populations, population, populations!, population!, exclude_loci, remove_loci, omit_loci, exclude_samples, remove_samples, omit_samples, samples

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
Takes **decimal degrees** as a `Vector` of any `AbstractFloat`.
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
    locations!(data, long = long_adjust, lat = lat_adjust)
end


"""
    locations!(data::PopData; long::Vector{String}, lat::Vector{String})
Replaces existing `PopData` location data (longitude `long`, latitude `lat`). Takes
**decimal minutes** or **degrees minutes seconds** format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.
## Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
### NOTE
If you read in the coordinate data as 4 vectors (longitude degrees, longitude minutes, latitude degrees, latitude minutes),
then the easiest course of action would be to merge them into two vectors of strings
(one for longitude, one for latitude):
```
long_string = string.(lat_deg, " ", lat_min)
lat_string = string.(long_deg, " ", long_min)
```
and use these as inputs into `locations!`

### Example
```
ncats = nancycats();
x = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)
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

function locations!(data::PopData; kwargs...)
    kwargs = Dict(kwargs)
    # check for matching lat and long keywords
    if all([haskey(kwargs, :lat), haskey(kwargs, :long)])
        locations!(data, kwargs[:long], kwargs[:lat])
    else
        error("keyword arguments \"lat\" and \"long\" must be supplied together")
    end
end


"""
    loci(data::PopData)
Returns an array of strings of the loci names in a `PopData` object.
"""
function loci(data::PopData)
    levels(data.loci.locus)
end

"""
    get_genotypes(data::PopObj; samples::Union{String, Vector{String}}, loci::Union{String, Vector{String}})
Return the genotype(s) of one or more `samples` for one or more
specific `loci` (both as keywords) in a `PopData` object.
### Examples
```
cats = nancycats();
get_genotype(cats, sample = "N115" , locus = "fca8")
get_genotypes(cats, sample = ["N1", "N2"] , locus = "fca8")
get_genotype(cats, sample = "N115" , locus = ["fca8", "fca37"])
get_genotype(cats, sample = ["N1", "N2"] , locus = ["fca8", "fca37"])
```
"""
function get_genotypes(data::PopData, sample::Union{String, Vector{String}}, locus::Union{String, Vector{String}})
    if typeof(sample) == String
        sample = [sample]
    end
    if typeof(locus) == String
        locus = [locus]
    end
    @where(data.loci, :name .∈ Ref(sample), :locus .∈ Ref(locus))
end


"""
    get_sample_genotypes(data::PopData, sample::String)
Return all the genotypes of a specific sample in a `PopData` object.
```
cats = nancycats()
get_sample_genotypes(cats, "N115")
```
"""
function get_sample_genotypes(data::PopObj, sample::String)
    @where(data.loci, :name .== sample)[! , :genotype]
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
    @where(data.loci, :locus .== locus)[! , :genotype]
end


#### Find missing ####
"""
    missing(data::PopData; by::String = "sample")
Get missing genotype information in a `PopData`. Specify a mode of operation
to return a DataFrame corresponding with that missing information.

#### Modes
- "sample" - returns a count and list of missing loci per individual (default)
- "pop" - returns a count of missing genotypes per population
- "locus" - returns a count of missing genotypes per locus
- "full" - returns a count of missing genotypes per locus per population

### Example:
```
missing(gulfsharks(), by = "pop")
```
"""
@inline function Base.missing(data::PopData; by::String = "sample")
    if by ∈ ["sample", "individual"]
        return @by(data.loci, :name, missing = count(ismissing, :genotype))

    elseif by ∈ ["pop", "population"]
        return @by(data.loci, :population, missing = count(ismissing, :genotype))

    elseif by ∈ ["locus", "loci"]
        return @by(data.loci, :locus, missing = count(ismissing, :genotype))

    elseif by ∈ ["detailed", "full"]
        return @by(data.loci, [:locus, :population], missing = count(ismissing, :genotype))
    else
        @error "Mode \"$by\" not recognized. Please specify one of: sample, pop, locus, or full"
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
        return @select(data.meta, :name, :population)
    else
        if all(ismissing.(data.meta.population)) == true
            @info "no population data present in PopData"
            return populations(data, listall = true)
        end
        return @by(data.meta, :population, count = length(:population))
    end
end

#TODO all the populations! functions
"""
```
# Replace by matching
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```
Multiple methods to rename or reassign population names to a `PopData`.

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
## Reassign using samples and new population assignments
Completely reassign populations for each individual. Takes two vectors of strings
as input: one of the sample names, and the other with their new corresponding
population name.

\n**Example**
```
populations!(potatoes, ["potato_1", "potato_2"], ["north_russet", "south_russet"])
```

"""
function populations!(data::PopData, rename::Dict)
    msg = ""
    @inbounds for key in keys(rename)
        if key ∉ levels(data.loci.population)
            msg *= "  Population \"$key\" not found in PopData\n"
        else
            replace!(data.meta.population, key => rename[key])
        end
    end
    msg != "" && printstyled("Warnings:", color = :yellow) ; print("\n"*msg)
    recode!(data.loci.population,rename...)
    return
end

function populations!(data::PopData, rename::Vector{String})
    current_popnames = unique(data.meta.population)
    rn_dict = Dict{String, String}()
    [rn_dict[string(i)] = j for (i,j) in zip(current_popnames, rename)]
    populations!(data, rn_dict)
    return
end

function populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
    for (i,j) in zip(samples, populations)
        data.meta[data.meta.name .== i, :population] .= j
        data.loci[data.loci.name .== i, :population] .= j
    end
    droplevels!(data.loci.population)
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
    new_table = @where(data.loci, :locus .!= locus)
    droplevels!(new_table.locus)
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
    new_table =  @where data.loci :locus .∉ Ref(exloci)
    msg != "" && printstyled("Warnings:", color = :yellow) ; println(msg)
    droplevels!(new_table.locus)
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
    new_loc_table = @where(data.loci, :name .!= samp_id)
    new_meta_table = @where(data.meta, :name .!= samp_id)
    droplevels!(new_loc_table.locus)
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
    new_loc_table = @where(data.loci, :name .∉  Ref(samp_ids))
    new_meta_table = @where(data.meta, :name .∉  Ref(samp_ids))
    msg != "" && printstyled("Warnings:", color = :yellow) ; println(msg)
    droplevels!(new_loc_table.locus)
    return PopData(new_meta_table, new_loc_table)
end

const omit_samples = exclude_samples
const remove_samples = exclude_samples

"""
    samples(data::PopData)
View individual/sample names in a `PopData`
"""
function samples(data::PopData)
    data.meta.name
end
