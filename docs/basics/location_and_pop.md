---
id: populations
title: Location and population data
sidebar_label: Location and population data
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

Population genetics involves a focus on... populations (gasp!). The commands below show you how to view and modify both population information (names), and location information (geographic coordinates). 

## Location Data

### View location data

```julia
locations(data::PopData)
```

View location data (`.longitude` and `.latitude`) in a `PopData`,  returning a table the longitude and latitude information in `meta`. 

```julia
julia> locations(sharks)
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

### Add geographical coordinates
```julia
locations!(data::PopData; long::Vector{T}, lat::Vector{T}) where T<:AbstractFloat
locations!(data::PopData; long::Vector{String}, lat::Vector{String})
```

Location data can be added using one of the methods of `locations!`. As indicated by the bang `!`, your `PopData` will be edited in place, and there will be no return output. If your data is in anything other than Decimal-Degrees format, this function will convert your long/lat into Decimal Degrees. To import those data into Julia, you'll likely want to use the wonderful `CSV.jl` package first. The functions accept keywords `long` and `lat`, or can be used without them so long as the vectors are input in that order. 

<Tabs
  block={true}
  defaultValue="dm"
  values={[
    { label: 'decimal minutes', value: 'dm', },
    { label: 'other formats', value: 'other', },
  ]
}>
<TabItem value="dm">

This method is pretty straightforward and tolerates vectors with `missing` data.
#### Formatting requirements
- Coordinates must be decimal-minutes as either `Float32` or `Float64` (e.g. `-21.321`)

```julia
# generate some fake location data
julia> long = rand(212) .* 10 ; lat = rand(212) .* -10

julia> locations!(sharks, long, lat)
```

</TabItem>
<TabItem value="other">

#### Formatting requirements

- Coordinates as `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)

If not already in decimal-minutes format, it would likely be most convenient if you imported your coordinate data as vectors of strings, which would look something like this:

```julia
long = ["-43 54 11", "22 23 11N"]
lat = ["11 44 31", "-25 41 13"]
```

:::caution Missing values
This method tolerates `missing` values, but you will need to `replace!` instances of `missing` with the string `"missing"`.
:::

</TabItem>
</Tabs>

## Population Names

### View population names

```julia
populations(data::PopData; listall::Bool = false)
```

Just as you can view population names with `PopData.meta.population`, you can also view them with the `populations` command, which by default shows you a summary of the number of individuals in each population.  

<Tabs
  block={true}
  defaultValue="pop"
  values={[
    { label: 'populations', value: 'pop', },
    { label: 'listall = true', value: 'all', },
  ]
}>
<TabItem value="pop">

``` julia
julia> populations(sharks)
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

</TabItem>
<TabItem value="all">

You can use the keyword `listall = true` to display each individual and their associated population as a table. 

``` julia
julia> populations(sharks, listall=true)
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

</TabItem>
</Tabs>

:::note alias
You can use the command `population` for the same functionality. We made the commands `population` and `populations` synonymous so you wouldn't have to memorize if the name was singular or plural-- it just works! This also applies to `populations!` and `population!`
:::

### Rename populations
```julia
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```

There are a handful of methods to alter `PopData` population names depending on what you find most convenient. Each of these methods start with `populations!()` and vary in their inputs. It's for that reason this function has an uncharacteristically long docstring. However, all the methods for `populations!` are unified in that they edit `PopData` in place, and print (rather than return) a table of the new population names and counts courtesy of `populations()`.

<Tabs
  block={true}
  defaultValue="dict"
  values={[
    { label: 'using a Dictionary', value: 'dict', },
    { label: 'using a Vector of names', value: 'vec', },
	{ label: 'reassign by sample', value: 'samp', },
  ]
}>
<TabItem value="dict">

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

</TabItem>
<TabItem value="vec">

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

</TabItem>
<TabItem value="samp">

```julia
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```

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

Now we can combine them with `populations!` to rename the first 5 Cape Canaveral samples.

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

</TabItem>
</Tabs>