---
id: viewdata
title: Viewing data
sidebar_label: Viewing data
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code "tabs".

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

Now that we have nancycats loaded in, we can use standard Julia accessor conventions to view the elements within our PopData. The DataFrames uses the convention `dataframe.colname` to directly access the columns we want.

## The metadata (data about the data)
Some critical information about the data is front-loaded into a PopData object to eliminate constantly getting these values in calculations.
To view this information, use `metadata()`.
```
julia> metadata(ncats)
 ploidy:      2
 loci:        9
 samples:     237
 populations: 17
 biallelic:   false
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

### sampleinfo

To view the sample information, you can use `sampleinfo()`

```julia
julia> sampleinfo(ncats)
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
julia> sinfo = sampleinfo(ncats) ;
julia> sinfo.name
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

### locusinfo

To view the locus information, you can use `locusinfo()`. Locus information is not mandatory,
but present if needed for future analyses.

```julia
julia> locusinfo(ncats)
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

-----

## The genotype table

### genodata

You can view the genotype information with `genodata()`.

```julia
julia> genodata(ncats)
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

Because the genotype data is in long format (aka "tidy"), accessing genotypes in a meaningful way is fairly 
straightforward if you have any experience with dataframe manipulation. For a deeper look into indexing `PopData`,
read [Advanced PopData Indexing](advancedindexing) 



The functions here help you inspect your `PopData` and pull information from it easily.

## View specific information

### sample names

```julia
samplenames(data::PopData)
```
View individual/sample names in a `PopData`. 
``` julia
julia> samplenames(sharks)
212-element Array{String,1}:
 "cc_001" 
 "cc_002" 
 "cc_003" 
 "cc_005" 
 "cc_007" 
 ⋮        
 "seg_027"
 "seg_028"
 "seg_029"
 "seg_030"
 "seg_031"
```

### locus names
```julia
loci(data::PopData)
```
Returns a vector of strings of the loci names in a `PopData`
```julia
julia> loci(sharks)
2213-element Array{String,1}:
 "contig_35208"
 "contig_23109"
 "contig_4493" 
 "contig_10742"
 "contig_14898"
 ⋮             
 "contig_43517"
 "contig_27356"
 "contig_475"  
 "contig_19384"
 "contig_22368"
 "contig_2784" 
```

## View genotypes
### all genotypes in one locus or sample
```julia
genotypes(data::PopData, samplelocus::String)
```
Returns a vector (view) of genotypes for a locus, or sample, depending on which the function finds in your data. Don't worry too much about
the wild type signature of the return vector. 

``` julia
julia> genotypes(sharks, "contig_2784")
212-element view(::PooledArrays.PooledVector{Union{Missing, Tuple{Int8, Int8}}, UInt8, Vector{UInt8}}, [468097, 468098, 468099, 468100, 468101, 468102, 468103, 468104, 468105, 468106  …  468299, 468300, 468301, 468302, 468303, 468304, 468305, 468306, 468307, 468308]) with eltype Union{Missing, Tuple{Int8, Int8}}:
 (1, 1)
 (1, 1)
 (1, 1)
 ⋮
 (1, 1)
 (1, 1)
 (1, 1)


julia> genotypes(sharks, "cc_001")
2209-element view(::PooledArrays.PooledVector{Union{Missing, Tuple{Int8, Int8}}, UInt8, Vector{UInt8}}, [1, 213, 425, 637, 849, 1061, 1273, 1485, 1697, 1909  …  466189, 466401, 466613, 466825, 467037, 467249, 467461, 467673, 467885, 468097]) with eltype Union{Missing, Tuple{Int8, Int8}}:
 (1, 2)
 (1, 1)
 (1, 2)
 ⋮
 (2, 2)
 (1, 1)
 (1, 1)
```

### one sample, one locus
```julia
genotype(data::PopData, sample::String => locus::String)
```
Returns the genotype of the `sample` at the `locus`. Uses `Pair` notation.
```julia
julia> genotype(sharks, "cc_001" => "contig_2784")
(1, 1)
```

### many samples, one locus
```julia
genotype(data::PopData, samples::Vector{String} => loci::String)
```
Returns a subdataframe of the genotypes of the `samples` at the `locus`. Uses `Pair` notation.
```julia
julia> genotypes(sharks, samplenames(sharks)[1:3] => "contig_2784")
3×4 SubDataFrame
 Row │ name     population     locus        genotype 
     │ String7  String         String       Tuple…?  
─────┼───────────────────────────────────────────────
   1 │ cc_001   CapeCanaveral  contig_2784  (1, 1)
   2 │ cc_002   CapeCanaveral  contig_2784  (1, 1)
   3 │ cc_003   CapeCanaveral  contig_2784  (1, 1)
```

### one sample, many loci
```julia
genotype(data::PopData, sample::String => loci::Vector{String})
```
Returns a subdataframe of the genotypes of the `sample` at the `loci`. Uses `Pair` notation.
```julia
julia> genotypes(sharks, "cc_001" => loci(sharks)[1:3])
3×4 SubDataFrame
 Row │ name     population     locus         genotype 
     │ String7  String         String        Tuple…?  
─────┼────────────────────────────────────────────────
   1 │ cc_001   CapeCanaveral  contig_35208  (1, 2)
   2 │ cc_001   CapeCanaveral  contig_23109  (1, 1)
   3 │ cc_001   CapeCanaveral  contig_4493   (1, 2)
```

### many samples, many loci
```julia
genotype(data::PopData, samples::Vector{String} => loci::Vector{String})
```
Returns a subdataframe of the genotypes of the `samples` at the `loci`. Uses `Pair` notation.

```julia
julia> genotypes(sharks, samplenames(sharks)[1:3] => loci(sharks)[1:3])
9×4 SubDataFrame
 Row │ name     population     locus         genotype 
     │ String7  String         String        Tuple…?  
─────┼────────────────────────────────────────────────
   1 │ cc_001   CapeCanaveral  contig_35208  (1, 2)
   2 │ cc_002   CapeCanaveral  contig_35208  (1, 2)
   3 │ cc_003   CapeCanaveral  contig_35208  (1, 1)
   4 │ cc_001   CapeCanaveral  contig_23109  (1, 1)
   5 │ cc_002   CapeCanaveral  contig_23109  (1, 2)
   6 │ cc_003   CapeCanaveral  contig_23109  missing  
   7 │ cc_001   CapeCanaveral  contig_4493   (1, 2)
   8 │ cc_002   CapeCanaveral  contig_4493   (1, 1)
   9 │ cc_003   CapeCanaveral  contig_4493   (1, 1)
```
