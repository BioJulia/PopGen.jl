# Directly accessing elements

A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code "tabs", where the left-most tab is the input named after what it's accomplishing, and the right tab is the output of running the command. This guide is to show you how to directly access  `PopData` elements, but there are shortcut commands to view just about every element of the data within. 

::: danger don't manually edit or sort
There are specific relationships between the record entries in `PopData` objects, so **do not use** `sort`, `sort!`, or manually arrange/add/delete anything in PopData. There are included functions to remove samples or loci, rename things, add location data, etc. 
:::
## Loading in the data

Let's keep things simple by loading in the nancycats data and calling it `ncats`.


:::: tabs card stretch

::: tab load nancycats
``` julia
julia> ncats = nancycats()
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


Now that we have nancycats loaded in, we can use standard Julia accessor conventions to view the elements within our PopData. The IndexedTable format requires a little extra work, so we must use the convention `PopData.meta.colname` to directly access the columns we want.

## The metadata table

### .meta

To view the entire `meta` table.

:::: tabs card stretch

::: tab PopData meta field
```julia
julia> ncats.meta
```
:::

::: tab output
```julia
237×5 DataFrame
│ Row │ name   │ population │ ploidy │ latitude │ longitude │
│     │ String │ String     │ Int8   │ Float32? │ Float32?  │
├─────┼────────┼────────────┼────────┼──────────┼───────────┤
│ 1   │ N215   │ 1          │ 2      │ missing  │ missing   │
│ 2   │ N216   │ 1          │ 2      │ missing  │ missing   │
│ 3   │ N217   │ 1          │ 2      │ missing  │ missing   │
│ 4   │ N218   │ 1          │ 2      │ missing  │ missing   │
│ 5   │ N219   │ 1          │ 2      │ missing  │ missing   │
│ 6   │ N220   │ 1          │ 2      │ missing  │ missing   │
⋮
│ 231 │ N294   │ 17         │ 2      │ missing  │ missing   │
│ 232 │ N295   │ 17         │ 2      │ missing  │ missing   │
│ 233 │ N296   │ 17         │ 2      │ missing  │ missing   │
│ 234 │ N297   │ 17         │ 2      │ missing  │ missing   │
│ 235 │ N281   │ 17         │ 2      │ missing  │ missing   │
│ 236 │ N289   │ 17         │ 2      │ missing  │ missing   │
│ 237 │ N290   │ 17         │ 2      │ missing  │ missing   │
```
:::

::::

### .name

This will access the names of the samples.


:::: tabs card stretch

::: tab meta .name field
``` julia
julia> ncats.meta.name
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


:::: tabs card stretch
::: tab meta .population field
``` julia
julia> ncats.meta.population
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

:::: tabs card stretch
::: tab meta .ploidy
``` julia
julia> ncats.meta.ploidy
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

:::: tabs card stretch
::: tab meta .latitude field
```julia
julia> ncats.meta.latitude
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

:::: tabs card stretch
::: tab meta .longitude field
```julia
julia> ncats.meta.longitude
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

:::: tabs card stretch
::: tab actually seeing some location info 
The nancycats data has some weird coordinate system for information, so those data were omitted. If you want a proof of concept for `.longitude` and `.latitude`, load in `gulfsharks` and try it out. We'll use `DataFrames.select` to isolate just the information we want. Later, you'll see that the `locations` command does this.
:::    
::: tab gulfsharks location data

``` julia
julia> sharks = gulfsharks() ;    # semicolon just supresses printing output

julia> select(sharks.meta, :name, :longitude, :latitude)
212×3 DataFrame
│ Row │ name    │ longitude │ latitude │
│     │ String  │ Float64   │ Float64  │
├─────┼─────────┼───────────┼──────────┤
│ 1   │ cc_001  │ -80.5993  │ 28.3062  │
│ 2   │ cc_002  │ -80.5995  │ 28.3079  │
│ 3   │ cc_003  │ -80.5996  │ 28.3023  │
│ 4   │ cc_005  │ -80.4225  │ 28.6123  │
│ 5   │ cc_007  │ -80.3578  │ 27.8666  │
│ 6   │ cc_008  │ -80.3579  │ 27.8666  │
⋮
│ 206 │ seg_025 │ -86.5374  │ 30.064   │
│ 207 │ seg_026 │ -86.5376  │ 30.0696  │
│ 208 │ seg_027 │ -86.0905  │ 29.9065  │
│ 209 │ seg_028 │ -87.3661  │ 30.0532  │
│ 210 │ seg_029 │ -87.3662  │ 30.0522  │
│ 211 │ seg_030 │ -85.7143  │ 29.8234  │
│ 212 │ seg_031 │ -85.7143  │ 29.8234  │
```
:::
::::
--------------------

## The genotype table

### .loci

This will show you the entire `loci` table.

:::: tabs card stretch
::: tab PopData loci field 
```julia
julia> ncats.loci
```
:::

::: tab output
```
2133×4 DataFrame
│ Row  │ name │ population │ locus │ genotype   │
│      │ Cat… │ Cat…       │ Cat…  │ Tuple…?    │
├──────┼──────┼────────────┼───────┼────────────┤
│ 1    │ N215 │ 1          │ fca8  │ missing    │
│ 2    │ N216 │ 1          │ fca8  │ missing    │
│ 3    │ N217 │ 1          │ fca8  │ (135, 143) │
│ 4    │ N218 │ 1          │ fca8  │ (133, 135) │
│ 5    │ N219 │ 1          │ fca8  │ (133, 135) │
│ 6    │ N220 │ 1          │ fca8  │ (135, 143) │
⋮
│ 2127 │ N294 │ 17         │ fca37 │ (208, 208) │
│ 2128 │ N295 │ 17         │ fca37 │ (208, 208) │
│ 2129 │ N296 │ 17         │ fca37 │ (208, 220) │
│ 2130 │ N297 │ 17         │ fca37 │ (208, 208) │
│ 2131 │ N281 │ 17         │ fca37 │ (208, 208) │
│ 2132 │ N289 │ 17         │ fca37 │ (208, 208) │
│ 2133 │ N290 │ 17         │ fca37 │ (208, 208) │
```
:::
::::

### locus names

This will access the names of the loci as they appear in the data. Since everything but the genotypes in `.loci` are coded as Categorical, we need to use `levels()` from `CategoricalArrays.jl` or `unique()` from Base to pull out the unique loci. 

:::: tabs card stretch
::: tab loci .locus
```julia
julia> levels(ncats.loci.locus)
# or #
julia> unique(ncats.loci.locus)
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

Because the genotype data is in "tidy" format, accessing genotypes in a meaningful way is not immediately obvious. We can of course follow the same convention of `data.loci.genotype` as we have above, but a list of all the genotypes across all individuals, loci, and populations isn't terribly useful. Instead, we can use the DataFrames or DataFramesMeta interfaces to retrieve this information. Here is an example using `@where` from DataFramesMeta:

:::: tabs card stretch
::: tab single locus
```julia
julia> @where(ncats.loci, :locus .== "fca8")
2133×4 DataFrame
│ Row  │ name │ population │ locus │ genotype   │
│      │ Cat… │ Cat…       │ Cat…  │ Tuple…?    │
├──────┼──────┼────────────┼───────┼────────────┤
│ 1    │ N215 │ 1          │ fca8  │ missing    │
│ 2    │ N216 │ 1          │ fca8  │ missing    │
│ 3    │ N217 │ 1          │ fca8  │ (135, 143) │
│ 4    │ N218 │ 1          │ fca8  │ (133, 135) │
│ 5    │ N219 │ 1          │ fca8  │ (133, 135) │
│ 6    │ N220 │ 1          │ fca8  │ (135, 143) │
⋮
│ 2127 │ N294 │ 17         │ fca8  │ (208, 208) │
│ 2128 │ N295 │ 17         │ fca8  │ (208, 208) │
│ 2129 │ N296 │ 17         │ fca8  │ (208, 220) │
│ 2130 │ N297 │ 17         │ fca8  │ (208, 208) │
│ 2131 │ N281 │ 17         │ fca8  │ (208, 208) │
│ 2132 │ N289 │ 17         │ fca8  │ (208, 208) │
│ 2133 │ N290 │ 17         │ fca8  │ (208, 208) │
```
:::
::: tab multiple loci
```julia
julia> cats.loci[ncats.loci.locus .∈ Ref(["fca8", "fca23"]), :]
474×4 DataFrame
│ Row │ name │ population │ locus │ genotype   │
│     │ Cat… │ Cat…       │ Cat…  │ Tuple…?    │
├─────┼──────┼────────────┼───────┼────────────┤
│ 1   │ N215 │ 1          │ fca8  │ missing    │
│ 2   │ N216 │ 1          │ fca8  │ missing    │
│ 3   │ N217 │ 1          │ fca8  │ (135, 143) │
│ 4   │ N218 │ 1          │ fca8  │ (133, 135) │
│ 5   │ N219 │ 1          │ fca8  │ (133, 135) │
│ 6   │ N220 │ 1          │ fca8  │ (135, 143) │
⋮
│ 468 │ N294 │ 17         │ fca23 │ (136, 146) │
│ 469 │ N295 │ 17         │ fca23 │ (130, 136) │
│ 470 │ N296 │ 17         │ fca23 │ (136, 146) │
│ 471 │ N297 │ 17         │ fca23 │ (130, 130) │
│ 472 │ N281 │ 17         │ fca23 │ (136, 144) │
│ 473 │ N289 │ 17         │ fca23 │ (130, 136) │
│ 474 │ N290 │ 17         │ fca23 │ (130, 146) │
```
:::
::::

Now that you're somewhat familiar with the parts of `PopData`, [have a look at the commands](view_and_sort.md) to view and manipulate `PopData` objects.