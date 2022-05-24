---
id: dataexploration
title: DataExploration.jl
sidebar_label: DataExploration.jl
---

## PopGen.jl/src/DataExplortation.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 🔵 allelefreqtable
```julia
allelefreqtable(data::PopData; by::String = "global")
```
Return a table of the observed `global` (default) or `population` allele frequencies in a PopData object.

**Example**
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

----

### 🔵 genofreqtable
```julia
genofreqtable(data::PopData; by::String = "global")
```
Return a table of the observed `global` (default) or `population` genotype frequencies in a PopData object.

**Example**
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

----

### 🔵 missingdata
```julia
missingdata(data::PopData; by::String = "sample")
```
Get missing genotype information in a `PopData`. Specify a mode of operation to return a DataFrame corresponding with that missing information.

**Modes**
- `"sample"` - returns a count and list of missing loci per individual (default)
- `"pop"` - returns a count of missing genotypes per population
- `"locus"` - returns a count of missing genotypes per locus
- `"full"` - returns a count of missing genotypes per locus per population

**Example**
```julia
missingdata(@gulfsharks, by = "pop")
```

-----

### 🔵 pairwiseidentical
```julia
pairwiseidentical(data::PopData)
```
Return a pairwise matrix of the percent of identical genotypes at each locus between all pairs of individuals.

**Example**
```
julia> cats = @nancycats ;
julia> pairwiseidentical(cats)
237×237 Named Matrix{Float64}
A ╲ B │     N215      N216  …      N289      N290
──────┼──────────────────────────────────────────
N215  │      1.0       0.5  …  0.142857  0.166667
N216  │      0.5       1.0     0.142857  0.166667
N217  │     0.25     0.125        0.125  0.142857
N218  │    0.375      0.25         0.25  0.142857
N219  │    0.375     0.375         0.25  0.142857
⋮              ⋮         ⋮  ⋱         ⋮         ⋮
N296  │      0.5  0.333333          0.0       0.0
N297  │ 0.166667  0.166667     0.428571  0.285714
N281  │ 0.142857  0.142857         0.25  0.428571
N289  │ 0.142857  0.142857          1.0  0.142857
N290  │ 0.166667  0.166667  …  0.142857       1.0
```

```julia
pairwiseidentical(data::PopData, sample_names::Vector{String})
```
Return a pairwise matrix of the percent of identical genotypes at 
each nonmissing locus between all pairs of provided `sample_names`.

**Example**
```julia
julia> cats = @nancycats ;
julia> interesting_cats = samplenames(cats)[1:5]
5-element Array{String,1}:
 "N215"
 "N216"
 "N217"
 "N218"
 "N219"
 
julia> pairwiseidentical(cats, interesting_cats)
5×5 Named Matrix{Float64}
A ╲ B │     N217      N218      N219      N220      N221
──────┼─────────────────────────────────────────────────
N217  │      1.0       0.0  0.111111  0.222222  0.111111
N218  │      0.0       1.0  0.333333  0.111111  0.444444
N219  │ 0.111111  0.333333       1.0  0.111111  0.333333
N220  │ 0.222222  0.111111  0.111111       1.0  0.222222
N221  │ 0.111111  0.444444  0.333333  0.222222       1.0
```
