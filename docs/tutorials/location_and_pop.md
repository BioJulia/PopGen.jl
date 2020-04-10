# Location and population data

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
Table with 212 rows, 2 columns:
longitude  latitude
───────────────────
28.3062    -80.5993
28.3079    -80.5995
28.3023    -80.5996
28.6123    -80.4225
27.8666    -80.3578
27.8666    -80.3579
27.8682    -80.3482
27.8711    -80.3482
⋮
30.0696    -86.5376
29.9065    -86.0905
30.0532    -87.3661
30.0522    -87.3662
29.8234    -85.7143
29.8234    -85.7143
```
:::
::::

### Add location data
Location data can be added using one of the methods of `locations!`. As indicated by the bang `!`, your `PopData` will be edited in place, and there will be no return output. If your data is in Decimal Minutes format, this function will convert your long/lat into Decimal Degrees. To import those data into Julia, you'll likely want to use the wonderful `CSV.jl` package first. 

- Decimal Degrees : `-11.431`
- Decimal Minutes : `-11 43.11` (notice the space)

::: warning Must use minus sign
Your data **must** use the minus sign `-` (if appropriate) instead of cardinal directions. `11 43.11W` is **not** valid.
:::

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
It would likely be most convenient if you imported your decimal minutes data as vectors of strings, which would look something like this:
```julia
lat = ["11 44.31", "-25 41.94"]
long = ["-43 54.11", "22 23.11"]
```

For this, the method is

```julia
locations!(data::PopData; lat::Vector{String}, long::Vector{String})
```

which uses the `lat` and `long` keywords.

::: warning Missing values
This method tolerates `missing` values, but you will need to `replace!` the string `"missing"` with values of `missing`.
:::

::: tab Decimal minutes as separate vectors
Alternatively, can input four vectors of numbers with the associated keyword arguments:

| Input                       | Type    | Keyword Argument |
| --------------------------- | ------- | ---------------- |
| Vector of longitude degrees | Integer | `long_deg`       |
| Vector of longitude minutes | Float   | `long_min`       |
| Vector of latitude degrees  | Integer | `lat_deg`        |
| Vector of latitude minutes  | Float   | `lat_min`        |

This method is easier or more tedious depending on what you consider a more practical approach. For example, if you have decimal-minutes coordinates for two samples:

|          | Longitude | Latitude |
| -------- | --------- | -------- |
| Sample 1 | 11 43.12  | 15 36.53 |
| Sample 2 | -12 41.32 | 11 22.41 |

then your inputs would be:

```
lo_deg = [11, -12]
lo_min = [43.12, 41.32]
la_deg  = [15, 11]
la_min  = [36.53, 22.41]
```

and you would then use `locations!` like this:

```julia
locations!(data, long_deg = lo_deg, long_min = lo_min, lat_deg = la_deg, lat_min = la_min)
```
:::
::::


## Population Names

### View population names

```julia
populations(data::PopData; listall::Bool = false)
```

Just as you can view population names with `PopData.meta.columns.population`, you can also view them with the `populations` command, which by default shows you a summary of the number of individuals in each population.  

:::: tabs card stretch
::: tab populations
``` julia
julia> populations(sharks)
```
:::
::: tab output
```
Table with 7 rows, 2 columns:
population        count
───────────────────────
"Cape Canaveral"  21
"Florida Keys"    65
"Georgia"         30
"Mideast Gulf"    28
"Northeast Gulf"  20
"South Carolina"  28
"Southeast Gulf"  20
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
Table with 212 rows, 2 columns:
name       population
───────────────────────────
"cc_001"   "Cape Canaveral"
"cc_002"   "Cape Canaveral"
"cc_003"   "Cape Canaveral"
"cc_005"   "Cape Canaveral"
"cc_007"   "Cape Canaveral"
"cc_008"   "Cape Canaveral"
"cc_009"   "Cape Canaveral"
"cc_010"   "Cape Canaveral"
⋮
"seg_026"  "Southeast Gulf"
"seg_027"  "Southeast Gulf"
"seg_028"  "Southeast Gulf"
"seg_029"  "Southeast Gulf"
"seg_030"  "Southeast Gulf"
"seg_031"  "Southeast Gulf"
```
:::
::::

::: details alias functions
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

Table with 2 rows, 2 columns:
population  count
─────────────────
"Atlantic"  79
"Gulf"      133
```
:::
::: tab Using a Vector of names
```julia
populations!(data::PopData, rename::Vector{String})
```

`Vector` of new unique population names in the order that they appear in the `PopData.meta`.

```julia
julia> new_popnames = ["Atlantic", "Atlantic", "Atlantic", "Gulf", "Gulf", "Gulf", "Gulf"]

julia> populations!(sharks, new_popnames)

Table with 2 rows, 2 columns:
population  count
─────────────────
"Atlantic"  79
"Gulf"      133
```
:::
::: tab Using a Vector of oldnames and new names
```julia
populations!(data::PopData, oldnames::Vector{String, newnames::Vector{String}})
```

Similar to the `Dict` method, except instead of creating a dictionary of "oldname" => "newname", you input a Vector{String} of `oldnames` followed by another of `newnames`. Logically, the new names will replace the old names in the same order as they appear in `PopData.meta` (e.g. the first newname replaces the first oldname, etc.).

```julia
julia> old_pop = ["Cape Canaveral", "Florida Keys", "Georgia", "Mideast Gulf", "Northeast Gulf", "South Carolina", "Southeast Gulf"] ;

julia> new_pop = ["Atlantic", "Gulf", "Atlantic", "Gulf", "Gulf", "Atlantic", "Gulf"] ;

julia> populations(sharks, old_pop, new_pop)
Table with 2 rows, 2 columns:
population  count
─────────────────
"Atlantic"  79
"Gulf"      133
```
:::

**********

#### Generate new population information

You may want outright overwrite all current population information. This is particularly useful when importing from VCF format when population information is not provided. This method will completely replace the population names of `PopData` regardless of what they currently are. 

::: warning Double-check your population counts
If you're playing along and getting errors that the lengths don't match, then get make sure you're using the right population counts. You can get those numbers with `populations(sharks)`.
:::

```julia
counts = [21, 65, 30, 28, 20, 28, 20]
```

and we then also create the vector of the names in the order in which they appear:

```julia
popnames = ["Cape Canaveral", "Florida Keys", "Georgia", "Mideast Gulf", "Northeast Gulf", "South Carolina", "Southeast Gulf"]
```

Now we can combine them with `populations!` to restore the population names to how they were originally
::::: tabs card stretch
::: tab Replace using a NamedTuple
```julia
julia> populations!(sharks, (counts = counts, names = popnames))
Table with 7 rows, 2 columns:
population        count
───────────────────────
"Cape Canaveral"  21
"Florida Keys"    65
"Georgia"         30
"Mideast Gulf"    28
"Northeast Gulf"  20
"South Carolina"  28
"Southeast Gulf"  20
```
:::
::: tab Replace Using Vectors of names & counts
This is just about the same as using the `NamedTuple`, but perhaps some users will prefer this format.
```julia
julia> populations!(sharks, popnames, counts)
Table with 7 rows, 2 columns:
population        count
───────────────────────
"Cape Canaveral"  21
"Florida Keys"    65
"Georgia"         30
"Mideast Gulf"    28
"Northeast Gulf"  20
"South Carolina"  28
"Southeast Gulf"  20
```
:::
::::
