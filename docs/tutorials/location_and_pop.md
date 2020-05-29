---
id: populations
title: Location and population data
sidebar_label: Location and population data
---

Population genetics involves a focus on... populations (gasp!). The commands below show you how to view and modify both population information (names), and location information (geographic coordinates). 

## Location Data

### View location data

```julia
locations(data::PopData)
```

View location data (`.longitude` and `.latitude`) in a `PopData`,  returning a table the longitude and latitude information in `meta`. 

:::: tabs card stretch
::: tab locations

```julia
julia> locations(sharks)
```

:::
::: tab output

```
212×2 SubDataFrame
│ Row │ longitude │ latitude │
│     │ Float64   │ Float64  │
├─────┼───────────┼──────────┤
│ 1   │ -80.5993  │ 28.3062  │
│ 2   │ -80.5995  │ 28.3079  │
│ 3   │ -80.5996  │ 28.3023  │
│ 4   │ -80.4225  │ 28.6123  │
│ 5   │ -80.3578  │ 27.8666  │
│ 6   │ -80.3579  │ 27.8666  │
⋮
│ 206 │ -86.5374  │ 30.064   │
│ 207 │ -86.5376  │ 30.0696  │
│ 208 │ -86.0905  │ 29.9065  │
│ 209 │ -87.3661  │ 30.0532  │
│ 210 │ -87.3662  │ 30.0522  │
│ 211 │ -85.7143  │ 29.8234  │
│ 212 │ -85.7143  │ 29.8234  │
```

:::
::::

### Add location data

Location data can be added using one of the methods of `locations!`. As indicated by the bang `!`, your `PopData` will be edited in place, and there will be no return output. If your data is in anything other than Decimal-Degrees format, this function will convert your long/lat into Decimal Degrees. To import those data into Julia, you'll likely want to use the wonderful `CSV.jl` package first. 

#### Formatting requirements

- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)

There are three main ways of adding location data:
:::: tabs card stretch
::: tab Already in decimal degrees

```julia
locations!(data::PopObj, long::Vector{T}, lat::Vector{T}) where T<:AbstractFloat
```

This method is pretty straightforward, and it tolerates vectors with `missing` data.

```julia
# generate some fake location data
julia> long = rand(212) .* 10 ; lat = rand(212) .* -10

julia> locations!(sharks, long, lat)
```

:::
::: tab Decimal minutes as strings
It would likely be most convenient if you imported your coordinate data as vectors of strings, which would look something like this:

```julia
long = ["-43 54 11", "22 23 11N"]
lat = ["11 44 31", "-25 41 94"]
```

For this, the method is

```julia
locations!(data::PopData; long::Vector{String}, lat::Vector{String})
```

which uses the `long` and `lat` keywords.

:::caution Missing values
This method tolerates `missing` values, but you will need to `replace!` instances of `missing` with the string `"missing"`.
:::
::::


## Population Names

### View population names

```julia
populations(data::PopData; listall::Bool = false)
```

Just as you can view population names with `PopData.meta.population`, you can also view them with the `populations` command, which by default shows you a summary of the number of individuals in each population.  

:::: tabs card stretch
::: tab populations

``` julia
julia> populations(sharks)
```

:::
::: tab output

```
7×2 DataFrame
│ Row │ population     │ count │
│     │ String         │ Int64 │
├─────┼────────────────┼───────┤
│ 1   │ Cape Canaveral │ 21    │
│ 2   │ Georgia        │ 30    │
│ 3   │ South Carolina │ 28    │
│ 4   │ Florida Keys   │ 65    │
│ 5   │ Mideast Gulf   │ 28    │
│ 6   │ Northeast Gulf │ 20    │
│ 7   │ Southeast Gulf │ 20    │
```

:::
::::

You can use the keyword `listall=true` to display each individual and their associated population as a table. 
:::: tabs card stretch
::: tab listall = true

``` julia
julia> populations(sharks, listall=true)
```

:::
::: tab output

```
212×2 DataFrame
│ Row │ name    │ population     │
│     │ String  │ String         │
├─────┼─────────┼────────────────┤
│ 1   │ cc_001  │ Cape Canaveral │
│ 2   │ cc_002  │ Cape Canaveral │
│ 3   │ cc_003  │ Cape Canaveral │
│ 4   │ cc_005  │ Cape Canaveral │
│ 5   │ cc_007  │ Cape Canaveral │
│ 6   │ cc_008  │ Cape Canaveral │
⋮
│ 206 │ seg_025 │ Southeast Gulf │
│ 207 │ seg_026 │ Southeast Gulf │
│ 208 │ seg_027 │ Southeast Gulf │
│ 209 │ seg_028 │ Southeast Gulf │
│ 210 │ seg_029 │ Southeast Gulf │
│ 211 │ seg_030 │ Southeast Gulf │
│ 212 │ seg_031 │ Southeast Gulf │
```

:::
::::

:::note alias
You can use the command `population` for the same functionality. We made the commands `population` and `populations` synonymous so you wouldn't have to memorize if the name was singular or plural-- it just works! This also applies to `populations!` and `population!`
:::

### Rename populations

There are a handful of methods to alter `PopData` population names depending on what you find most convenient. Each of these methods start with `populations!()` and vary in their inputs. It's for that reason this function has an obnoxiously long docstring. For simplicity, the methods will be separated into categories. However, all the methods for `populations!` are unified in that they edit `PopData` in place, and print (rather than return) a table of the new population names and counts courtesy of `populations()`.

#### Replace by matching

These methods require that some kind of population information is already present, in the sense that the samples in `PopData` aren't all in one population.
:::: tabs card stretch
::: tab using a Dictionary

```julia
populations!(data::PopData, rename::Dict)
```

Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`.

``` julia
# create a dictionary of name conversions
julia> new_popnames = Dict(
    		"Cape Canaveral" => "Atlantic",
			"Georgia" => "Atlantic",
			"South Carolina" => "Atlantic",
    		"Florida Keys" => "Gulf",
    		"Mideast Gulf" => "Gulf",
    		"Northeast Gulf" => "Gulf",
    		"Southeast Gulf" => "Gulf"
		);	

julia> populations!(sharks, new_popnames)

2×2 DataFrame
│ Row │ population │ count │
│     │ String     │ Int64 │
├─────┼────────────┼───────┤
│ 1   │ Atlantic   │ 79    │
│ 2   │ Gulf       │ 133   │
```

:::
::: tab Using a Vector of names

```julia
populations!(data::PopData, rename::Vector{String})
```

`Vector` of new unique population names in the order that they appear in the `PopData.meta`.

```julia
julia> new_popnames = ["Atlantic", "Atlantic", "Atlantic", "Gulf", "Gulf", "Gulf", "Gulf"] ;

julia> populations!(sharks, new_popnames)

2×2 DataFrame
│ Row │ population │ count │
│     │ String     │ Int64 │
├─────┼────────────┼───────┤
│ 1   │ Atlantic   │ 79    │
│ 2   │ Gulf       │ 133   │
```

:::

**********

#### Reassign population information

You may want outright overwrite all current population information. This is particularly useful when importing from VCF format when population information is not provided. This method will completely replace the population names of `PopData` regardless of what they currently are. 

This method takes a vector of sample names and a vector of the new population names of the samples in the order that they appear in the name-vector.

```julia
# creating a vector of sample names
julia> ch_names = samples(sharks)[1:5]
5-element Array{String,1}:
 "cc_001"
 "cc_002"
 "cc_003"
 "cc_005"
 "cc_007"
```

and we then also create the vector of these samples' new population names:

```julia
julia> popnames = ["North Cape", "North Cape", "North Cape", "South Cape", "South Cape"] ;
```

Now we can combine them with `populations!` to restore the population names to how they were originally
::::: tabs card stretch
::: tab Replace using a NamedTuple

```julia
julia> populations!(x, ch_names, popnames)
212×2 DataFrame
│ Row │ name    │ population │
│     │ String  │ String     │
├─────┼─────────┼────────────┤
│ 1   │ cc_001  │ North Cape │
│ 2   │ cc_002  │ North Cape │
│ 3   │ cc_003  │ North Cape │
│ 4   │ cc_005  │ South Cape │
│ 5   │ cc_007  │ South Cape │
│ 6   │ cc_008  │ Atlantic   │
⋮
│ 206 │ seg_025 │ Gulf       │
│ 207 │ seg_026 │ Gulf       │
│ 208 │ seg_027 │ Gulf       │
│ 209 │ seg_028 │ Gulf       │
│ 210 │ seg_029 │ Gulf       │
│ 211 │ seg_030 │ Gulf       │
│ 212 │ seg_031 │ Gulf       │
```

:::

::::