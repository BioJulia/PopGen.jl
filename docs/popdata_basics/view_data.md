---
id: viewdata
title: Viewing Data
sidebar_label: Viewing data
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

The functions here help you inspect your `PopData` and pull information from it easily.

## Individuals / Samples

### view individuals' names

```julia
samples(data::PopData)
```

View individual/sample names in a `PopData`. 

``` julia
julia> samples(sharks)
212-element Array{String,1}:
 "cc_001" 
 "cc_002" 
 "cc_003" 
 "cc_005" 
 "cc_007" 
 "cc_008" 
 "cc_009" 
 "cc_010" 
 "cc_012" 
 "cc_013" 
 ⋮        
 "seg_023"
 "seg_024"
 "seg_025"
 "seg_026"
 "seg_027"
 "seg_028"
 "seg_029"
 "seg_030"
 "seg_031"
```

## Display Specific Loci and/or Samples

### Get loci names

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
 "contig_8483" 
 "contig_8065" 
 "contig_14708"
 "contig_2307" 
 ⋮             
 "contig_24711"
 "contig_18959"
 "contig_43517"
 "contig_27356"
 "contig_475"  
 "contig_19384"
 "contig_22368"
 "contig_2784" 
```

### View genotypes at a locus

```julia
locus(data::PopData, locus::String)
```

Returns a Vector of genotypes for a locus

``` julia
julia> locus(sharks, "contig_2784")
212-element Array{Union{Missing, Tuple{Int8,Int8}},1}:
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 ⋮
 (1, 1)
 (1, 1)
 (1, 1)
 (1, 1)
 missing
 (1, 1)
 (1, 1)
 (1, 1)
```