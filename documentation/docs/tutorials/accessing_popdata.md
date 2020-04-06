# Directly accessing elements

A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code "tabs", where the left-most tab is the input named after what it's accomplishing, and the right tab is the output of running the command. This guide is to show you how to directly access  `PopData` elements, but there are shortcut commands to view just about every element of the data within. 

::: danger don't manually edit or sort
There are specific relationships between the record entries in `PopData` objects, so **do not use** `sort`, `sort!`, or manually arrange/add/delete anything in PopData. There are included functions to remove samples or loci, rename things, add location data, etc. 
:::
## Loading in the data

Let's keep things simple by loading in the nancycats data and calling it `ncats`.


:::: tabs type:board-card stretch

::: tab load nancycats
``` julia
julia> ncats = nancycats() ; summary(ncats)
```
:::

::: tab output
```
PopData Object
  Marker type: Microsatellite
  Ploidy: 2
  Number of individuals: 237
  Number of loci: 9
  Populations: 17
  Longitude: absent
  Latitude: absent
```
:::
::::


Now that we have nancycats loaded in, we can use standard Julia accessor conventions to view the elements within our PopData. The IndexedTable format requires a little extra work, so we must use the convention `PopData.meta.columns.colname` to directly access the columns we want.

## The metadata table

### .meta

To view the entire `meta` table.

:::: tabs type:board-card stretch

::: tab PopData meta field
```julia
julia> ncats.meta
```
:::

::: tab output
```julia
Table with 237 rows, 5 columns:
name    population  ploidy  latitude  longitude
───────────────────────────────────────────────
"N1"    "1"         2       missing   missing
"N2"    "1"         2       missing   missing
"N3"    "1"         2       missing   missing
"N4"    "1"         2       missing   missing
"N5"    "1"         2       missing   missing
"N6"    "1"         2       missing   missing
"N7"    "1"         2       missing   missing
"N8"    "1"         2       missing   missing
⋮
"N231"  "17"        2       missing   missing
"N232"  "17"        2       missing   missing
"N233"  "17"        2       missing   missing
"N234"  "17"        2       missing   missing
"N235"  "17"        2       missing   missing
"N236"  "17"        2       missing   missing
"N237"  "17"        2       missing   missing
```
:::

::::

### .name

This will access the names of the samples.


:::: tabs type:board-card stretch

::: tab meta .name field
``` julia
julia> ncats.meta.columns.name
```
:::

::: tab output
```
237-element Array{String,1}:
 "N1"  
 "N2"  
 "N3"  
 "N4"  
 "N5"  
 "N6"  
 "N7"  
 "N8"  
 ⋮     
 "N230"
 "N231"
 "N232"
 "N233"
 "N234"
 "N235"
 "N236"
 "N237"
```
:::
::::

### .population
This will access the names of the populations associated with each sample, in the same order as the  samples.


:::: tabs type:board-card stretch
::: tab meta .population field
``` julia
julia> ncats.meta.columns.population
```
:::

::: tab output
```
237-element Array{String,1}:
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 ⋮   
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
```
:::
::::

These ID's aren't super informative. Later, we'll change them using the `popid!` command.

###  .ploidy

This shows you the ploidy of the data per individual

:::: tabs type:board-card stretch
::: tab meta .ploidy
``` julia
julia> ncats.meta.columns.ploidy
```
:::
::: tab output
```
237-element Array{Int8,1}:
 2
 2
 2
 2
 2
 2
 2
 2
 ⋮
 2
 2
 2
 2
 2
 2
 2
 2
```
:::
::::

### .latitude

This accesses the latitude information of the PopObj. If there is none, like in the nancycats data, it returns a vector of `missing`.

:::: tabs type:board-card stretch
::: tab meta .latitude field
```julia
julia> ncats.meta.columns.latitude
```
:::
::: tab output
```
237-element Array{Union{Missing, Float32},1}:
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 ⋮      
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
```
:::
::::

### .longitude

This accesses the longitude information of the PopObj. Like before, if there is none, like in the nancycats data, it returns an array of `missing`.

:::: tabs type:board-card stretch
::: tab meta .longitude field
```julia
julia> ncats.meta.columns.longitude
```
:::
::: tab output
```
237-element Array{Union{Missing, Float32},1}:
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 ⋮      
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
```
:::
::::

:::: tabs type:border-card stretch
::: tab actually seeing some location info 
The nancycats data has some weird coordinate system for information, so those data were omitted. If you want a proof of concept for `.longitude` and `.latitude`, load in `gulfsharks` and try it out. We'll use `hcat` (horizontal concatination) to horizontally bind the individual names, their latitude, and longitude. Later, you'll see that the `locations` command does this.
:::    
::: tab gulfsharks location data
``` julia
julia> sharks = gulfsharks() ;    # semicolon just supresses printing output

julia> hcat(sharks.meta.columns.name, sharks.meta.columns.latitude, sharks.meta.columns.longitude)

212×3 Array{Any,2}:
"cc_001"   28.3062  -80.5993
"cc_002"   28.3079  -80.5995
"cc_003"   28.3023  -80.5996
"cc_005"   28.6123  -80.4225
"cc_007"   27.8666  -80.3578
"cc_008"   27.8666  -80.3579
"cc_009"   27.8682  -80.3482
"cc_010"   27.8711  -80.3482
"cc_012"   28.4815  -80.4303
"cc_013"   28.2421  -80.2494
⋮                           
"seg_021"  29.9466  -86.0399
"seg_023"  29.9969  -85.6494
"seg_024"  29.6966  -87.4403
"seg_025"  30.064   -86.5374
"seg_026"  30.0696  -86.5376
"seg_027"  29.9065  -86.0905
"seg_028"  30.0532  -87.3661
"seg_029"  30.0522  -87.3662
"seg_030"  29.8234  -85.7143
"seg_031"  29.8234  -85.7143
```
:::
::::
--------------------

## The genotype table

### .loci

This will show you the entire `loci` table.

:::: tabs type:board-card stretch
::: tab PopData loci field 
```julia
julia> ncats.loci
```
:::

::: tab output
```
Table with 2133 rows, 4 columns:
name    population  locus    genotype
─────────────────────────────────────
"N1"    "1"         "fca8"   missing
"N1"    "1"         "fca23"  (4, 9)
"N1"    "1"         "fca43"  (4, 4)
"N1"    "1"         "fca45"  (1, 3)
"N1"    "1"         "fca77"  (9, 9)
"N1"    "1"         "fca78"  (3, 6)
"N1"    "1"         "fca90"  (9, 9)
"N1"    "1"         "fca96"  (8, 8)
⋮
"N237"  "17"        "fca45"  missing
"N237"  "17"        "fca77"  (6, 9)
"N237"  "17"        "fca78"  (3, 7)
"N237"  "17"        "fca90"  (8, 8)
"N237"  "17"        "fca96"  missing
"N237"  "17"        "fca37"  (10, 10)
```
:::
::::

### locus names

This will access the names of the loci as they appear in the data. Since everything but the genotypes in `.loci` are coded as Categorical, we need to use `levels()` from `CategoricalArrays.jl` or `unique()` from Base to pull out the unique loci. 

:::: tabs type:board-card stretch
::: tab loci .locus
```julia
julia> levels(ncats.loci.columns.locus)
# or #
julia> unique(ncats.loci.columns.locus)
```
:::
::: tab output
```
9-element Array{String,1}:
 "fca8" 
 "fca23"
 "fca43"
 "fca45"
 "fca77"
 "fca78"
 "fca90"
 "fca96"
 "fca37"
```
:::
::::


### view genotypes

Because the genotype data is in "tidy" format, accessing genotypes in a meaningful way is not immediately obvious. We can of course follow the same convention of `data.loci.genotype` as we have above, but a list of all the genotypes across all individuals, loci, and populations isn't terribly useful. Instead, we can use the JuliaDB or JuliaDBMeta interfaces to retrieve this information. Here is an example using `@where` from JuliaDBMeta:

:::: tabs type:board-card stretch
::: tab single locus
```julia
julia> @where ncats.loci :locus == "fca8"
Table with 237 rows, 4 columns:
name    population  locus   genotype
────────────────────────────────────
"N1"    "1"         "fca8"  missing
"N2"    "1"         "fca8"  missing
"N3"    "1"         "fca8"  (9, 13)
"N4"    "1"         "fca8"  (8, 9)
"N5"    "1"         "fca8"  (8, 9)
"N6"    "1"         "fca8"  (9, 13)
"N7"    "1"         "fca8"  (9, 9)
"N8"    "1"         "fca8"  (9, 13)
⋮
"N232"  "17"        "fca8"  (8, 12)
"N233"  "17"        "fca8"  (8, 12)
"N234"  "17"        "fca8"  (8, 13)
"N235"  "17"        "fca8"  (9, 12)
"N236"  "17"        "fca8"  (10, 13)
"N237"  "17"        "fca8"  (9, 12)
```
:::
::: tab multiple loci
```julia
julia> @where ncats.loci :locus in ["fca8", "fca23"]
Table with 474 rows, 4 columns:
name    population  locus    genotype
─────────────────────────────────────
"N1"    "1"         "fca8"   missing
"N1"    "1"         "fca23"  (4, 9)
"N2"    "1"         "fca8"   missing
"N2"    "1"         "fca23"  (9, 9)
"N3"    "1"         "fca8"   (9, 13)
"N3"    "1"         "fca23"  (4, 9)
"N4"    "1"         "fca8"   (8, 9)
"N4"    "1"         "fca23"  (5, 5)
⋮
"N235"  "17"        "fca8"   (9, 12)
"N235"  "17"        "fca23"  (4, 8)
"N236"  "17"        "fca8"   (10, 13)
"N236"  "17"        "fca23"  (2, 4)
"N237"  "17"        "fca8"   (9, 12)
"N237"  "17"        "fca23"  (2, 9)
```
:::
::::

Now that you're somewhat familiar with the parts of `PopData`, [have a look at the commands](view_and_sort.md) to view and manipulate `PopData` objects.