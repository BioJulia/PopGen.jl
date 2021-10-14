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
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
```

Now that we have nancycats loaded in, we can use standard Julia accessor conventions to view the elements within our PopData. The DataFrames uses the convention `PopData.sampleinfo.colname` to directly access the columns we want.

## The metadata (data about the data)
Some critical information about the data is front-loaded into a PopData object to eliminate constantly getting these values in calculations.
To view this information, use `popdata.metadata` or `popdata.info`. Each
of these fields can be directly accessed to retrieve those values (e.g. `popdata.metadata.samples`, `popdata.metadata.loci`, etc.)
```
julia> ncats.metadata
 ploidy:        2
 # loci:        9
 # samples:     237
 # populations: 17
 biallelic:     false
 ```

Included in `metadata` are two DataFrames, one for sample information, and another for locus information.
<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'sample information', value: 's', },
    { label: 'locus information', value: 'l', },
  ]
}>
<TabItem value="s">

### .sampleinfo

To view the sample information, you can use `popdata.metadata.sampleinfo`, or the shortcut `PopData.sampleinfo`

```julia
julia> ncats.sampleinfo
237×3 DataFrame
 Row │ name      population  ploidy 
     │ String7…  String      Int8   
─────┼──────────────────────────────
   1 │ N217      1                2
   2 │ N218      1                2
   3 │ N219      1                2
   4 │ N220      1                2
   5 │ N221      1                2
   6 │ N222      1                2
  ⋮  │    ⋮          ⋮         ⋮
 232 │ N197      14               2
 233 │ N198      14               2
 234 │ N199      14               2
 235 │ N200      14               2
 236 │ N201      14               2
 237 │ N206      14               2
                    222 rows omitted

```

Using the standard DataFrames `getindex` methods, we can access these columns like so:

``` julia
julia> ncats.sampleinfo.name
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

</TabItem>
<TabItem value="l">

### .locusinfo

To view the locus information, you can use `popdata.metadata.locusinfo`, of the shortcut `popdata.locusinfo`. Locus information is not mandatory,
but present if needed for future analyses.

```julia
julia> ncats.locusinfo
9×4 DataFrame
 Row │ chromosome  locus   cm       bp   
     │ Int8        String  Float64  Int64 
─────┼───────────────────────────────────
   1 │          0  fca8       0       0
   2 │          0  fca23      0       0
   3 │          0  fca43      0       0
   4 │          0  fca45      0       0
   5 │          0  fca77      0       0
   6 │          0  fca78      0       0
   7 │          0  fca90      0       0
   8 │          0  fca96      0       0
   9 │          0  fca37      0       0
```

</TabItem>
</Tabs>
--------------------

## The genotype table

### .genodata

This will show you the entire `genodata` table.

```julia
julia> ncats.genodata
2133×4 DataFrame
  Row │ name    population  locus   genotype
      │ String  String      String  Tuple…?
──────┼────────────────────────────────────────
    1 │ N215    1           fca8    missing
    2 │ N216    1           fca8    missing
    3 │ N217    1           fca8    (135, 143)
    4 │ N218    1           fca8    (133, 135)
    5 │ N219    1           fca8    (133, 135)
    6 │ N220    1           fca8    (135, 143)
  ⋮   │   ⋮         ⋮         ⋮         ⋮
 2128 │ N295    17          fca37   (208, 208)
 2129 │ N296    17          fca37   (208, 220)
 2130 │ N297    17          fca37   (208, 208)
 2131 │ N281    17          fca37   (208, 208)
 2132 │ N289    17          fca37   (208, 208)
 2133 │ N290    17          fca37   (208, 208)
                              2121 rows omitted
```

Because the genotype data is in "tidy" format, accessing genotypes in a meaningful way is fairly straightforward if you have any experience
with dataframe manipulation. For a deeper look into indexing `PopData`,
read [Advanced PopData Indexing](indexing) 

Now that you're somewhat familiar with the parts of `PopData`, have a look at [the commands](view_data.md) to view and manipulate `PopData` objects.
