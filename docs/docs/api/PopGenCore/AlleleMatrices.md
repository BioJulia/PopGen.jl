---
id: allelematrices
title: AlleleMatrices.jl
sidebar_label: AlleleMatrices.jl
---

## PopGenCore.jl/src/AlleleMatrices.jl
| 📦  not exported | 🟪  exported by PopGenCore.jl | 🔵  exported by PopGen.jl |
|:---:|:---:|:---:|

### 🟪 matrix
```julia
matrix(data::PopData, matrixtype::String = "frequency"; missings = "mean", scale = false, center = false)
```
Return a matrix of allele counts or frequencies per genotype where rows are samples
and columns are the occurence count or frequency of an allele for that locus in that sample.
Loci and alleles are sorted alphanumerically. The `scale` and `center` keywords are only relevant for allele frequencies, not counts.

**Positional Arguments**

- `data`: a PopData object
- `matrixtype`: a `String` or `Symbol` of `count` or `frequency` (default: `frequency`)

**Keyword Arguments**

- `missings`: a `String` denoting how to handle missing values when outputting `frequency` (default: `mean`)
    - `"missing"`: fallback method to keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')

**Example**
```
julia> cats = @nancycats ;
julia> cnts = matrix(cats, "count") ;  cnts[1:5,1:6]
5×6 Matrix{Union{Missing, Int8}}:
  missing   missing   missing   missing   missing   missing
  missing   missing   missing   missing   missing   missing
 0         0         0         0         0         0
 0         0         0         0         0         0
 0         0         0         0         0         0

julia> frq = matrix(cats, "frequency") ;  frq[1:5,1:6]
5×6 Matrix{Union{Missing, Float32}}:
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0

julia> frq = matrix(cats, :frequency, missings = "zero") ;  frq[1:5,1:6]
 5×6 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 
 julia> frq = matrix(cats, missings = "mean", scale = true, center = true) ;  frq[1:5,1:6]
 5×6 Matrix{Float32}:
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.07
```

----

### 📦 _setcounts

----

### 📦 _matrix
```julia
_matrix(data::PopData, ::Val{:count})
```
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
Missing values are preserved as `-1`.

----
```julia
_matrix(data::PopData, ::Val{:frequencyzero})
```
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by zeros.

---
```julia
_matrix(data::PopData, ::Val{:frequencymean})
```
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by the global mean allele frequency.

---
```julia
_matrix(data::PopData, ::Val{:frequencymissing})
```
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are kept as `missing`.

----

### 🟪 featurematrix
```julia
featurematrix(data::PopData, matrixtype::Union{String, Symbol} = "genotype")
```
**Positional Arguments**

    - `data`: a PopData object
    - `matrixtype`: a `String` or `Symbol` of `genotype`, or `allele` (default: `genotype`)

**genotype feature matrix**

Return a matrix of dummy-encoded genotypes (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`. For biallelic loci, `0` encodes homozygous for allele 1, `1` encodes for a heterozygote,
and `2` encodes for homozygous allele 2.

**allele feature matrix**

Return a matrix of dummy-encoded alleles (0,1), where rows correspond with samples and columns correspond to alleles within loci, such
that there are as many columns per locus as alleles for that locus. Missing alleles (from missing genotypes) are encoded as `-1`.

**Example**
```
julia> cats = @nancycats ;
julia> featurematrix(cats)
237×9 Matrix{Int8}:
 -1   0   0   0   0  0   0   0   0
 -1   1   1   1   0  0   1   0   0
  0   0   2   2   1  1   2   0   1
  1   2   3   3   2  0   0   1   0
  1   3   4   4   3  0   3   0   0
  ⋮                  ⋮
 49   0   1  -1  36  0   0  -1  13
 48   6   8  -1  25  1   2  -1   0
 29   9  23  -1   7  3  26  14   0
  3   5   8  -1   2  1   3  14   0
 29  10  16  -1   4  3   2  -1   0

julia> featurematrix(cats, "allele")
237×108 Matrix{Int8}:
 -1  -1  -1  -1  -1  …  0  0  0  0  0  0
 -1  -1  -1  -1  -1     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  ⋮                  ⋱           ⋮
  0   0   0   0   0     0  0  0  1  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0  …  0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0 
```

---
```julia
_featurematrix(data::PopData, ::Val{:genotype})
```
Return a matrix of dummy-encoded genotypes (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`. For biallelic loci, `0` encodes homozygous for allele 1, `1` encodes for a heterozygote,
and `2` encodes for homozygous allele 2.

---
```julia
_featurematrix(data::PopData, ::Val{:allele})
```
Return a matrix of dummy-encoded alleles (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`.