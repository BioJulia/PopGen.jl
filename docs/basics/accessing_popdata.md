---
id: accessing
title: Accessing elements
sidebar_label: Accessing elements
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code "tabs", where the left-most tab is the input named after what it's accomplishing, and the right tab is the output of running the command. This guide is to show you how to directly access  `PopData` elements, but there are shortcut commands to view just about every element of the data within. 

:::danger don't manually edit or sort
There are specific relationships between the record entries in `PopData` objects, so **do not use** `sort`, `sort!`, or manually arrange/add/delete anything in PopData. There are included functions to remove samples or loci, rename things, add location data, etc. 
:::

## Loading in the data

Let's keep things simple by loading in the nancycats data and calling it `ncats`.


``` julia
julia> ncats = @nancycats
PopData Object
  Marker type: Microsatellite
  Ploidy: 2
  Number of individuals: 237
  Number of loci: 9
  Populations: 17
  Longitude: absent
  Latitude: absent
```

Now that we have nancycats loaded in, we can use standard Julia accessor conventions to view the elements within our PopData. The DataFrames uses the convention `PopData.meta.colname` to directly access the columns we want.

## The metadata table

### .meta

To view the entire `meta` table.

```julia
julia> ncats.meta
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

### .name

This will access the names of the samples.

``` julia
julia> ncats.meta.name
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

### .population

This will access the names of the populations associated with each sample, in the same order as the  samples.

``` julia
julia> ncats.meta.population
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

These ID's aren't super informative. Later, we'll change them using the `populations!` command.

### .ploidy

This shows you the ploidy of the data per individual

``` julia
julia> ncats.meta.ploidy
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

### .latitude

This accesses the latitude information of the PopObj. If there is none, like in the nancycats data, it returns a vector of `missing`.

```julia
julia> ncats.meta.latitude
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

### .longitude

This accesses the longitude information of the PopObj. Like before, if there is none, like in the nancycats data, it returns an array of `missing`.


```julia
julia> ncats.meta.longitude
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



:::note actually seeing some location info 
The nancycats data has some weird coordinate system for information, so those data were omitted. If you want a proof of concept for `.longitude` and `.latitude`, load in `gulfsharks` and try it out. We'll use `DataFrames.select` to isolate just the information we want. Later, you'll see that the `locations` command does this.
``` julia
julia> sharks = @gulfsharks ;    # semicolon just supresses printing output

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

--------------------

## The genotype table

### .loci

This will show you the entire `loci` table.

```julia
julia> ncats.loci
2133×4 DataFrame
│ Row  │ name │ population │ locus │ genotype   │
│      │ Str… │ String     │ Str…  │ Tuple…?    │
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

### locus names

This will access the names of the loci as they appear in the data. We need to use  `unique()` from Base to pull out the unique loci. 

```julia
julia> unique(ncats.loci.locus)
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

### locus population

This will access the population of the individual for a genotype of a locus. We need to use `unique()` from Base to pull out the unique populations. 

```julia
julia> unique(ncats.loci.locus)
17-element Array{String,1}:
 "1" 
 "2"
 "3"
 "4"
 "5"
 ⋮
 "14"
 "15"
 "16"
 "17"
```

### view genotypes

Because the genotype data is in "tidy" format, accessing genotypes in a meaningful way is not immediately obvious. We can of course follow the same convention of `data.loci.genotype` as we have above, but a list of all the genotypes across all individuals, loci, and populations isn't terribly useful. Instead, we can use the DataFrames or DataFramesMeta interfaces to retrieve this information. Here are examples using `@where` from DataFramesMeta:
<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'single locus', value: 's', },
    { label: 'multiple loci', value: 'm', },
    { label: 'population', value: 'p', },
  ]
}>
<TabItem value="s">

```julia
julia> @where(ncats.loci, :locus .== "fca8")
2133×4 DataFrame
│ Row  │ name │ population │ locus │ genotype   │
│      │ Str… │ Str…       │ Str…  │ Tuple…?    │
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

</TabItem>
<TabItem value="m">

```julia
julia> @where(ncats.loci, :locus .∈ Ref(["fca8", "fca23"])
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

</TabItem>
<TabItem value="p">

```julia
julia> @where(x.loci, :population .== "9")
81×4 DataFrames.DataFrame
│ Row │ name │ population │ locus │ genotype   │
│     │ Cat… │ Cat…       │ Cat…  │ Tuple…?    │
├─────┼──────┼────────────┼───────┼────────────┤
│ 1   │ N104 │ 9          │ fca8  │ (121, 135) │
│ 2   │ N105 │ 9          │ fca8  │ (137, 137) │
│ 3   │ N106 │ 9          │ fca8  │ (135, 135) │
│ 4   │ N107 │ 9          │ fca8  │ (121, 135) │
│ 5   │ N108 │ 9          │ fca8  │ (135, 135) │
│ 6   │ N109 │ 9          │ fca8  │ (137, 137) │
│ 7   │ N111 │ 9          │ fca8  │ (135, 135) │
⋮
│ 74  │ N105 │ 9          │ fca37 │ (182, 182) │
│ 75  │ N106 │ 9          │ fca37 │ (182, 208) │
│ 76  │ N107 │ 9          │ fca37 │ (182, 208) │
│ 77  │ N108 │ 9          │ fca37 │ (208, 208) │
│ 78  │ N109 │ 9          │ fca37 │ (182, 208) │
│ 79  │ N111 │ 9          │ fca37 │ (208, 214) │
│ 80  │ N112 │ 9          │ fca37 │ (182, 206) │
│ 81  │ N113 │ 9          │ fca37 │ (208, 214) │
```

</TabItem>
</Tabs>

Now that you're somewhat familiar with the parts of `PopData`, [have a look at the commands](view_data.md) to view and manipulate `PopData` objects.