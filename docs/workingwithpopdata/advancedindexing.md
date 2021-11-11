---
id: advancedindexing
title: Advanced PopData Indexing
sidebar_label: Advanced Indexing
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

In version `0.7.0`, we introduce a powerful new way to index PopData...
by directly piggy-backing off of the incredible work done in `DataFrames.jl`. 
Now, you can index and subset PopData [almost] as though you were directly
subsetting the `genodata` dataframe, and it will return a new subsetted
PopData object (or other stuff). We'll go through some examples using the `nancycats` data.
The conceptual syntax (the arrows are for demonstration) looks like:
```julia
# return new PopData
popdata[column -> condition]

# return a new genodata table
popdata[column -> condition, :]

# return a specific column
popdata[column -> condition, :column]

# return specific columns
popdata[column -> condition, [:col1, :col2]]
```

```julia
#julia> using PopGen
julia> ncats = @nancycats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
```

### Basic conditional indexing 
Basic conditional indexing is a fancy way of saying "pulling out specific
information". Let's say we wanted to omit locus `fca8`.

```julia
julia> ncats[genodata(ncats).locus .!= "fca8"]
PopData{Diploid, 8 Microsatellite loci}
  Samples: 237
  Populations: 17
```

Or, maybe we only want loci `fca8` and `fca23`. We use the `∈` (`\in<TAB>`) operator and wrap the loci in `Ref()` to keep the set from being broadcasted.

```julia
julia> ncats[genodata(ncats).locus .∈  Ref(["fca8", "fca23"])]
PopData{Diploid, 2 Microsatellite loci}
  Samples: 237
  Populations: 17
```

Perhaps we want only populations 1 through 5. Again, we bind the set in `Ref()` to prevent broadcasting over its elements. We also need to change the integers to strings because population names are always strings.
```julia
julia> ncats[genodata(ncats).population .∈  Ref(string.(1:5))]
PopData{Diploid, 9 Microsatellite loci}
  Samples: 82
  Populations: 5
```

Maybe we just wanted to know the names of the samples in population `5`. Although for something like this you can just as well index the `sampleinfo` dataframe.
```julia
julia> ncats[genodata(ncats).population .== "5", :name] |> unique
15-element Vector{InlineStrings.String7}:
 "N55"
 "N56"
 "N57"
 "N58"
 "N59"
 "N60"
 "N61"
 "N62"
 "N63"
 "N64"
 "N65"
 "N66"
 "N67"
 "N68"
 "N69"
```

### Advanced conditional indexing
Just like in `DataFrames.jl`, we can chain conditions with a broadcasted 
"and" operator (`.&`) and really pull out information of interest. This also works for a broadcasted
"or" operator (`.|`). Something to keep in mind is that each statement needs to be wrapped in
parentheses like:
```julia
popdata[(statement1) .& (statement2)]
```

Let's find all the samples in population `2` that are heterozygous for allele `133` in locus `fca8` and return just a dataframe.
```julia
julia> gd = genodata(ncats) ;
julia> ncats[(gd.locus .== "fca8") .& (gd.population .== "2") .& (ishet.(gd.genotype, 133)), :]

6×4 DataFrame
 Row │ name      population  locus   genotype   
     │ String7…  String      String  Tuple…?    
─────┼──────────────────────────────────────────
   1 │ N141      2           fca8    (129, 133)
   2 │ N142      2           fca8    (129, 133)
   3 │ N146      2           fca8    (129, 133)
   4 │ N151      2           fca8    (129, 133)
   5 │ N154      2           fca8    (133, 135)
   6 │ N155      2           fca8    (131, 133)
```

How about which samples are missing data for locus `fca8`?
```julia
julia> gd = genodata(ncats) ;
julia> ncats[(gd.locus .== "fca8") .& (ismissing.(gd.genotype)), :name]
20-element PooledArrays.PooledVector{InlineStrings.String7, UInt8, Vector{UInt8}}:
 "N215"
 "N216"
 "N188"
 "N189"
 "N190"
 "N191"
 "N192"
 ⋮
 "N197"
 "N198"
 "N199"
 "N200"
 "N201"
 "N206"
```

This should get you started on thinking of ways to explore your data :smile:. 
