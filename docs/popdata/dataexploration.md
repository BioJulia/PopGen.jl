---
id: dataexploration
title: Data Exploration
sidebar_label: Data exploration
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

## Allele frequency table
```julia
allelefreqtable(data::PopData; by::Union{String, Symbol} = "global")
```
Return a table of the observed `global` (default) or `population` allele frequencies in a PopData object. Use this if you want to see what the frequencies are for every allele at every locus.

<Tabs
  block={true}
  defaultValue="g"
  values={[
    { label: 'global', value: 'g', },
    { label: 'population', value: 'p', },
  ]
}>
<TabItem value="g">

```julia
julia> cats = @nancycats ;

julia> allelefreqtable(cats)
108×4 DataFrame
 Row │ locus   allele  count  frequency  
     │ String  Int16?  Int64  Float64    
─────┼───────────────────────────────────
   1 │ fca8       135    105  0.241935
   2 │ fca8       143     44  0.101382
   3 │ fca8       133     33  0.0760369
   4 │ fca8       137     83  0.191244
  ⋮  │   ⋮       ⋮       ⋮        ⋮
 105 │ fca37      226      2  0.00421941
 106 │ fca37      216      7  0.0147679
 107 │ fca37      224      2  0.00421941
 108 │ fca37      204      6  0.0126582
                         100 rows omitted
```

</TabItem>
<TabItem value="p">

```julia
julia> cats = @nancycats ;

julia> allelefreqtable(cats, by = "population")
839×5 DataFrame
 Row │ locus   population  allele  count  frequency 
     │ String  String      Int16?  Int64  Float64   
─────┼──────────────────────────────────────────────
   1 │ fca8    1              135      9  0.5625
   2 │ fca8    1              143      4  0.25
   3 │ fca8    1              133      2  0.125
   4 │ fca8    1              137      1  0.0625
  ⋮  │   ⋮         ⋮         ⋮       ⋮        ⋮
 836 │ fca37   16             210      5  0.208333
 837 │ fca37   17             208     22  0.846154
 838 │ fca37   17             182      3  0.115385
 839 │ fca37   17             220      1  0.0384615
                                    831 rows omitted
```

</TabItem>
</Tabs>

## Genotype frequency table
```julia
genofreqtable(data::PopData; by::Union{String, Symbol} = "global")
```
Return a table of the observed `global` (default) or `population` genotype frequencies in a PopData object. Use this if you want to see what the frequencies are for every genotype at every locus.

<Tabs
  block={true}
  defaultValue="g"
  values={[
    { label: 'global', value: 'g', },
    { label: 'population', value: 'p', },
  ]
}>
<TabItem value="g">

```julia
julia> cats = @nancycats ;

julia> genofreqtable(cats)
341×4 DataFrame
 Row │ locus   genotype    count  frequency  
     │ String  Tuple…      Int64  Float64    
─────┼───────────────────────────────────────
   1 │ fca8    (135, 143)     16  0.0737327
   2 │ fca8    (133, 135)      9  0.0414747
   3 │ fca8    (135, 135)     23  0.105991
   4 │ fca8    (137, 143)      8  0.0368664
  ⋮  │   ⋮         ⋮         ⋮        ⋮
 338 │ fca37   (206, 220)      1  0.00421941
 339 │ fca37   (208, 218)      1  0.00421941
 340 │ fca37   (184, 184)      3  0.0126582
 341 │ fca37   (208, 210)      3  0.0126582
                             333 rows omitted
```

</TabItem>
<TabItem value="p">

```julia
julia> cats = @nancycats ;

julia> genofreqtable(cats, by = "population")
1094×5 DataFrame
  Row │ locus   population  genotype    count  frequency         
      │ String  String      Tuple…      Int64  Float64           
──────┼──────────────────────────────────────────────────        
    1 │ fca8    1           (135, 143)      3  0.375
    2 │ fca8    1           (133, 135)      2  0.25
    3 │ fca8    1           (135, 135)      2  0.25
    4 │ fca8    1           (137, 143)      1  0.125
  ⋮   │   ⋮         ⋮           ⋮         ⋮        ⋮
 1091 │ fca37   17          (208, 208)     10  0.769231
 1092 │ fca37   17          (182, 182)      1  0.0769231
 1093 │ fca37   17          (182, 208)      1  0.0769231
 1094 │ fca37   17          (208, 220)      1  0.0769231
                                        1086 rows omitted 
```

</TabItem>
</Tabs>

## Missing Data

```julia
missingdata(data::PopData; by::Union{String, Symbol} = "sample")
```

Get missing genotype information in a `PopData` object. Specify a mode of operation using the `by =` keyword to return a table corresponding with that missing information type.

|     by                |  what it does                                                  |
| :-------------------: |  ------------------------------------------------------------  |
| `"sample"`            |  returns a count of missing loci per individual (default)      |
|  `"population"`       | returns a count of missing genotypes per population            |
| `"locus"`             |  returns a count of missing genotypes per locus                |
|  `"locusxpopulation"` |  returns a count of missing genotypes per locus per population |

<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'sample', value: 's', },
    { label: 'population', value: 'p', },
    { label: 'locus', value: 'l', },
    { label: 'locusxpopulation', value: 'f', },
  ]
}>
<TabItem value="s">

```
julia> sharks = @gulfsharks ;

julia> missingdata(sharks)
212×2 DataFrame
 Row │ name     missing
     │ String   Int64
─────┼──────────────────
   1 │ cc_001       124
   2 │ cc_002        94
   3 │ cc_003       100
   4 │ cc_005         0
   5 │ cc_007         2
   6 │ cc_008         1
   7 │ cc_009         2
  ⋮  │    ⋮        ⋮
 206 │ seg_025        0
 207 │ seg_026        0
 208 │ seg_027        2
 209 │ seg_028       25
 210 │ seg_029        0
 211 │ seg_030        1
 212 │ seg_031        1
        198 rows omitted
```

</TabItem>
<TabItem value="p">

```
julia> sharks = @gulfsharks ;

julia> missingdata(sharks, by = "population")
7×2 DataFrame
 Row │ population      missing
     │ String          Int64
─────┼─────────────────────────
   1 │ Cape Canaveral      666
   2 │ Georgia             423
   3 │ South Carolina      233
   4 │ Florida Keys       1241
   5 │ Mideast Gulf         99
   6 │ Northeast Gulf      472
   7 │ Southeast Gulf     1504
```

</TabItem>
<TabItem value="l">s

```
julia> sharks = @gulfsharks ;

julia> missingdata(sharks, by = "locus")
2209×2 DataFrame
  Row │ locus         missing
      │ String        Int64
──────┼───────────────────────
    1 │ contig_35208        0
    2 │ contig_23109        6
    3 │ contig_4493         3
    4 │ contig_10742        2
    5 │ contig_14898        0
    6 │ contig_8483         0
    7 │ contig_8065         0
  ⋮   │      ⋮           ⋮
 2203 │ contig_18959        0
 2204 │ contig_43517        6
 2205 │ contig_27356        2
 2206 │ contig_475          0
 2207 │ contig_19384        5
 2208 │ contig_22368        3
 2209 │ contig_2784         7
             2195 rows omitted
```

</TabItem>
<TabItem value="f">

```
julia> sharks = @gulfsharks ;

julia> missingdata(sharks, by = "locusxpopulation")
15463×3 DataFrame
   Row │ locus         population      missing
       │ String        String          Int64
───────┼───────────────────────────────────────
     1 │ contig_35208  Cape Canaveral        0
     2 │ contig_35208  Georgia               0
     3 │ contig_35208  South Carolina        0
     4 │ contig_35208  Florida Keys          0
     5 │ contig_35208  Mideast Gulf          0
     6 │ contig_35208  Northeast Gulf        0
     7 │ contig_35208  Southeast Gulf        0
   ⋮   │      ⋮              ⋮            ⋮
 15457 │ contig_2784   Cape Canaveral        0
 15458 │ contig_2784   Georgia               2
 15459 │ contig_2784   South Carolina        1
 15460 │ contig_2784   Florida Keys          2
 15461 │ contig_2784   Mideast Gulf          1
 15462 │ contig_2784   Northeast Gulf        0
 15463 │ contig_2784   Southeast Gulf        1
                             15449 rows omitted
```

</TabItem>
</Tabs>


## Pairwise Identical Genotypes
While not a substitute for a [kinship analysis](docs/analyses/kinship), it may be useful to know or verify how similar your data are in a very literal sense:
how many identical genotypes do two individuals have across all loci? To do this, we use `pairwiseidentical()` to perform an all x all comparison of identical genotypes. This can be done for all individuals in a `PopData` object, or restricted to a specific set of individuals:

<Tabs
  block={true}
  defaultValue="a"
  values={[
    { label: 'all samples', value: 'a', },
    { label: 'some samples', value: 's', },
  ]
}>
<TabItem value="a">

```julia
julia> cats = @nancycats;

julia> pairwiseidentical(cats)
27966×4 DataFrame
   Row │ sample_1  sample_2  identical  n     
       │ String    String    Float64    Int64 
───────┼──────────────────────────────────────
     1 │ N215      N216           0.5       8
     2 │ N215      N217           0.25      8
     3 │ N215      N218           0.38      8
     4 │ N215      N219           0.38      8
   ⋮   │    ⋮         ⋮          ⋮        ⋮
 27963 │ N297      N290           0.29      7
 27964 │ N281      N289           0.25      8
 27965 │ N281      N290           0.43      7
 27966 │ N289      N290           0.14      7
                            27958 rows omitted
```

</TabItem>
<TabItem value="s">

```julia
julia> cats = @nancycats;

julia> interesting_cats = samplenames(cats)[1:5]
5-element Array{String,1}:
 "N215"
 "N216"
 "N217"
 "N218"
 "N219"

julia> pairwiseidentical(cats, interesting_cats)
10×4 DataFrame
 Row │ sample_1  sample_2  identical  n     
     │ String    String    Float64    Int64 
─────┼──────────────────────────────────────
   1 │ N215      N216           0.5       8 
   2 │ N215      N217           0.25      8 
   3 │ N215      N218           0.38      8 
   4 │ N215      N219           0.38      8 
   5 │ N216      N217           0.12      8 
   6 │ N216      N218           0.25      8 
   7 │ N216      N219           0.38      8 
   8 │ N217      N218           0.0       9 
   9 │ N217      N219           0.11      9 
  10 │ N218      N219           0.33      9 
```

</TabItem>
</Tabs>

## Allelic Richness
If you were curious about allelic richness (number of alleles per locus), then you can use `richness()` to find that out. Use `by = "population"` to return a table by locus by population.

<Tabs
  block={true}
  defaultValue="l"
  values={[
    { label: 'by locus', value: 'l', },
    { label: 'by locusxpopulation', value: 'p', },
  ]
}>
<TabItem value="l">

```julia
julia> cats = @nancycats;

julia> richness(cats)
9×2 DataFrame
 Row │ locus   richness 
     │ String  Int64    
─────┼──────────────────
   1 │ fca8          16
   2 │ fca23         11
   3 │ fca43         10
   4 │ fca45          9
   5 │ fca77         12
   6 │ fca78          8
   7 │ fca90         12
   8 │ fca96         12
   9 │ fca37         18
```

</TabItem>
<TabItem value="p">

```julia
julia> richness(cats, by = "population")
153×3 DataFrame
 Row │ locus   population  richness 
     │ String  String      Int64    
─────┼──────────────────────────────
   1 │ fca8    1                  4
   2 │ fca8    2                  6
   3 │ fca8    3                  7
   4 │ fca8    4                 10
  ⋮  │   ⋮         ⋮          ⋮
 150 │ fca37   14                 3
 151 │ fca37   15                 4
 152 │ fca37   16                 3
 153 │ fca37   17                 3
                    145 rows omitted
```

</TabItem>
</Tabs>

## Average Number of Alleles
Similar to richness, if you wanted to know the average number of alleles per locus, use `alleleavg()`. Use `rounding = false` if you don't want the answer rounded to 4 decimal places.
```julia
julia> alleleavg(@nancycats)
(mean = 12.0, stdev = 0.2668)

julia> alleleavg(@nancycats, rounding = false)
(mean = 12.0, stdev = 0.2667968432263687)
```

## Summary Statistics
Population genetics famously includes all manner of coefficients with which to summarize data. Use `summary()` to view FST, DST, HT, etc. (like `Hierfstat::basic.stats`). 

<Tabs
  block={true}
  defaultValue="g"
  values={[
    { label: 'global', value: 'g', },
    { label: 'by locus', value: 'l', },
  ]
}>
<TabItem value="g">

```julia
julia> summary(@nancycats)
1×10 DataFrame
 Row │ Het_obs  HS       HT       DST      HT′      DST′     FST      FST′     FIS      DEST
     │ Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────
   1 │  0.6299   0.7083   0.7717   0.0634   0.7757   0.0674   0.0821   0.0869   0.1108    0.231
```

</TabItem>
<TabItem value="l">

```julia
julia> summary(@nancycats, by = "locus")
9×11 DataFrame
 Row │ locus   Het_obs  HS       HT       DST      HT′      DST′     FST      FST′     FIS      DEST
     │ String  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ fca8     0.667    0.779    0.8619   0.0829   0.8671   0.0881   0.0962   0.1016   0.1438   0.3987
   2 │ fca23    0.6838   0.7439   0.7994   0.0555   0.8029   0.0589   0.0694   0.0734   0.0809   0.2302
   3 │ fca43    0.6814   0.7442   0.7937   0.0495   0.7968   0.0526   0.0623   0.066    0.0844   0.2054
   4 │ fca45    0.71     0.7085   0.7642   0.0557   0.7679   0.0594   0.0729   0.0774  -0.0021   0.2039
   5 │ fca77    0.6295   0.7828   0.8659   0.0831   0.8711   0.0883   0.096    0.1014   0.1958   0.4067
   6 │ fca78    0.5773   0.6339   0.6773   0.0434   0.6801   0.0462   0.0641   0.0679   0.0893   0.1261
   7 │ fca90    0.6454   0.7408   0.8144   0.0736   0.819    0.0782   0.0904   0.0955   0.1287   0.3017
   8 │ fca96    0.6259   0.6747   0.7657   0.091    0.7714   0.0967   0.1189   0.1254   0.0723   0.2973
   9 │ fca37    0.4485   0.5671   0.6027   0.0356   0.6049   0.0379   0.0591   0.0626   0.2091   0.0874
```

</TabItem>
</Tabs>

:::tip prime symbol
The column names above use the unicode prime symbol `′` to better reflect the actual coefficient ("FST prime" etc.). To print that character, press `\prime<TAB>`, which reads "backslash, the word 'prime', and the TAB button".
:::
