---
id: populationdata
title: Population data
sidebar_label: Population data
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

Needless to say, population information is crucial for population genetics, so there are several handy tools for dealing with that information.
If you need to see the population for every sample, then use `sampleinfo(popdata)` to retrieve the dataframe containing sample information.

### View unique population names

```julia
populations(data::PopData; counts::Bool = false)
```
If `counts = false`, returns a Vector of the unique populations present in the `PopData`. If `counts = true`, returns a
table of sample counts per population.

<Tabs
  block={true}
  defaultValue="unq"
  values={[
    { label: 'unique populations', value: 'unq', },
    { label: 'counts per population', value: 'cou', },
  ]
}>
<TabItem value="unq">

Return a vector of the unique populations. 

``` julia
julia> populations(sharks)
7-element Array{String,1}:
 "CapeCanaveral"
 "Georgia"
 "SouthCarolina"
 "FloridaKeys"
 "MideastGulf"
 "NortheastGulf"
 "SoutheastGulf"
```

</TabItem>
<TabItem value="cou">

Retrun a table of the populations and their counts

``` julia
julia> populations(sharks, counts = true)
7×2 DataFrame
 Row │ population      count 
     │ String          Int64 
─────┼───────────────────────
   1 │ Cape Canaveral     21
   2 │ Georgia            30
   3 │ South Carolina     28
   4 │ Florida Keys       65
   5 │ Mideast Gulf       28
   6 │ Northeast Gulf     20
   7 │ Southeast Gulf     20
```

</TabItem>
</Tabs>



### Rename populations
```julia
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```

There are a handful of methods to alter `PopData` population names depending on what you find most convenient. Each of these methods start with `populations!()` and vary in their inputs. It's for that reason this function has an uncharacteristically long docstring. However, all the methods for `populations!` are unified in that they edit `PopData` in place.

<Tabs
  block={true}
  defaultValue="dict"
  values={[
    { label: 'with a Dictionary', value: 'dict', },
    { label: 'with a Vector of names', value: 'vec', },
	{ label: 'reassign by sample', value: 'samp', },
  ]
}>
<TabItem value="dict">

```julia
populations!(data::PopData, rename::Dict)
```

:::tip
Recommended for renaming existing populations	
:::

Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`.

``` julia
# create a dictionary of name conversions
julia> new_popnames = 
Dict(
  "CapeCanaveral" => "Atlantic",
  "Georgia" => "Atlantic",
  "SouthCarolina" => "Atlantic",
  "FloridaKeys" => "Gulf",
  "MideastGulf" => "Gulf",
  "NortheastGulf" => "Gulf",
  "SoutheastGulf" => "Gulf"
);	

julia> populations!(sharks, new_popnames)
julia> populations(sharks, counts = true)
2×2 DataFrame
 Row │ population  count 
     │ String      Int64 
─────┼───────────────────
   1 │ Atlantic       79
   2 │ Gulf          133
```

</TabItem>
<TabItem value="vec">

:::caution not recommended
These methods _are_ available, but the `Dict` method is recommended instead of (1) and the reassign-by-sample method is recommended
instead of (2)
:::

```julia
populations!(data::PopData, rename::Vector{String})
```

1. rename the unique populations
    - **condition**: `length(rename) == length(unique(populations))`
    - `rename` is a vector of new unique population names in the order that they appear in `sampleinfo(popdata)`.
2. rename the population association for every sample
    - **condition**: `length(rename) == length(samplenames(data))`
    - `rename` is a vector of new populations names for the samples in the order that they appear in `sampleinfo(popdata)`

```julia
julia> new_popnames = ["Atlantic", "Atlantic", "Atlantic", "Gulf", "Gulf", "Gulf", "Gulf"] ;

julia> populations!(sharks, new_popnames)
julia> populations(sharks, counts = true)
2×2 DataFrame
 Row │ population  count 
     │ String      Int64 
─────┼───────────────────
   1 │ Atlantic       79
   2 │ Gulf          133
```

</TabItem>
<TabItem value="samp">

:::tip
Recommended for assigning population ID's for specific samples.	
:::

```julia
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```

You may want outright overwrite all current population information. This is particularly useful when importing from VCF format when population information is not provided. This method will completely replace the population names of `PopData` regardless of what they currently are. 

This method takes a vector of sample names and a vector of the new population names of the samples in the order that they appear in the name-vector.

```julia
# creating a vector of sample names
julia> ch_names = samplenames(sharks)[1:5]
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
julia> populations!(sharks, ch_names, popnames)
julia> sampleinfo(sharks)[1:6,:]
6×5 DataFrame
 Row │ name     population     ploidy  longitude  latitude 
     │ String7  String         Int8    Float64    Float64  
─────┼─────────────────────────────────────────────────────
   1 │ cc_001   North Cape          2    28.3062  -80.5993
   2 │ cc_002   North Cape          2    28.3079  -80.5995
   3 │ cc_003   North Cape          2    28.3023  -80.5996
   4 │ cc_005   South Cape          2    28.6123  -80.4225
   5 │ cc_007   South Cape          2    27.8666  -80.3578
   6 │ cc_008   CapeCanaveral       2    27.8666  -80.3579
```

</TabItem>
</Tabs>
