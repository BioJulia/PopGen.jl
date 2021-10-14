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

## View genotypes by sample (or anything)

This can be done fairly easily using DataFramesMeta macro `@where`

<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'single sample', value: 's', },
    { label: 'multiple samples', value: 'm', },
    { label: 'name and locus', value: 'nl', },
    { label: 'advanced usage', value: 'advanced', },
  ]
}>
<TabItem value="s">

```julia
julia> @where(sharks.genodata, :name .== "cc_001")
2209×4 DataFrame
  Row │ name    population      locus         genotype
      │ String  String          String        Tuple…?
──────┼────────────────────────────────────────────────
    1 │ cc_001  Cape Canaveral  contig_35208  (1, 2)
    2 │ cc_001  Cape Canaveral  contig_23109  (1, 1)
    3 │ cc_001  Cape Canaveral  contig_4493   (1, 2)
    4 │ cc_001  Cape Canaveral  contig_10742  (1, 1)
    5 │ cc_001  Cape Canaveral  contig_14898  (1, 2)
    6 │ cc_001  Cape Canaveral  contig_8483   (1, 1)
  ⋮   │   ⋮           ⋮              ⋮           ⋮
 2205 │ cc_001  Cape Canaveral  contig_27356  (1, 1)
 2206 │ cc_001  Cape Canaveral  contig_475    (1, 2)
 2207 │ cc_001  Cape Canaveral  contig_19384  (2, 2)
 2208 │ cc_001  Cape Canaveral  contig_22368  (1, 1)
 2209 │ cc_001  Cape Canaveral  contig_2784   (1, 1)
                                      2198 rows omitted
```

</TabItem>
<TabItem value="m">

```julia
julia> @where(sharks.genodata, :name .∈ Ref(["cc_001", "cc_002"]))
4418×4 DataFrame
  Row │ name    population      locus         genotype
      │ String  String          String        Tuple…?
──────┼────────────────────────────────────────────────
    1 │ cc_001  Cape Canaveral  contig_35208  (1, 2)
    2 │ cc_002  Cape Canaveral  contig_35208  (1, 2)
    3 │ cc_001  Cape Canaveral  contig_23109  (1, 1)
    4 │ cc_002  Cape Canaveral  contig_23109  (1, 2)
    5 │ cc_001  Cape Canaveral  contig_4493   (1, 2)
    6 │ cc_002  Cape Canaveral  contig_4493   (1, 1)
  ⋮   │   ⋮           ⋮              ⋮           ⋮
 4414 │ cc_002  Cape Canaveral  contig_19384  (2, 2)
 4415 │ cc_001  Cape Canaveral  contig_22368  (1, 1)
 4416 │ cc_002  Cape Canaveral  contig_22368  (1, 1)
 4417 │ cc_001  Cape Canaveral  contig_2784   (1, 1)
 4418 │ cc_002  Cape Canaveral  contig_2784   (1, 1)
                                      4407 rows omitted
```

</TabItem>
<TabItem value="nl">

It also means that you can combine different queries with commas. Here is an example of an approach combining a name and locus criteria:

```julia
julia> @where(sharks.genodata, :name .∈ Ref(["cc_001", "cc_002"]), :locus .== "contig_2784")
2×4 DataFrame
│ Row │ name   │ population     │ locus       │ genotype │
│     │ Str…   │ String         │ String      │ Tuple…?  │
├─────┼────────┼────────────────┼─────────────┼──────────┤
│ 1   │ cc_001 │ Cape Canaveral │ contig_2784 │ (1, 1)   │
│ 2   │ cc_002 │ Cape Canaveral │ contig_2784 │ (1, 1)   │
```

</TabItem>
<TabItem value="advanced">
 
Here's an advanced example for writing a query that only returns heterozygous genotypes for locus `contig_1784`

```julia
julia> @where(sharks.genodata, ishet.(:genotype) .== true, :locus .== "contig_2784")
6×4 DataFrame
 Row │ name     population      locus        genotype
     │ String   String          String       Tuple…?
─────┼────────────────────────────────────────────────
   1 │ key_003  Florida Keys    contig_2784  (1, 2)
   2 │ key_055  Florida Keys    contig_2784  (1, 2)
   3 │ key_072  Florida Keys    contig_2784  (1, 2)
   4 │ meg_001  Mideast Gulf    contig_2784  (1, 2)
   5 │ meg_011  Mideast Gulf    contig_2784  (1, 2)
   6 │ seg_016  Southeast Gulf  contig_2784  (1, 2)
```

</TabItem>
</Tabs>
