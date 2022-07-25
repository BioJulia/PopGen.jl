---
id: pairwisekinship
title: PairwiseKinship.jl
sidebar_label: PairwiseKinship.jl
---
## PopGen.jl/src/Kinship/PairwiseKinship.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 🔵 kinship
```julia
kinship(data::PopData; method::Function, iterations::Int = 0, interval::Vector{Float64} = [0.025, 0.975])
kinship(data::PopData, samplenames::AbstractVector{T}; method::Function, iterations::Int = 0, interval::Vector{Float64} = [0.025, 0.975]) where T<:AbstractString
```
Calculate pairwise relatedness estimates for all or specific individuals in a `PopData` object using 
the specified `method` (see below). Returns a `NamedMatrix` if not performing bootstrapping, otherwise returns a `DataFrame` (since bootstrapping provides more output information). To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `[0.0275, 0.975]` (i.e. 95%), however that can be changed by supplying a `Vector{Float64}` of `[low, high]` to the keyword `interval`. **Note:** samples must be diploid.

**Arguments**
- `data` : A PopData object
- `samplenames`: Vector of sample names (optional)

**Keyword Arguments**
- `method::Function` : A method function (see below)
- `iterations::Int64` : The number of iterations to perform bootstrapping (default: `0`, will not perform bootstrapping)
- `interval::Vector{Float64}` : A Vector of [low, high] indicating the confidence intervals you would like for bootstrapping (default: `[0.275, 0.975]`, i.e. 95%)

**Methods**
| Method | Type | Method Call |
|:----|:-----|:-----|
| [Blouin 1996](https://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.1996.00094.x) | moments-based | `Blouin` |
| [Li & Horvitz 1953](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1716461/)| moments-based | `LiHorvitz` |
| [Loiselle 1995](https://bsapubs.onlinelibrary.wiley.com/doi/abs/10.1002/j.1537-2197.1995.tb12679.x) | moments-based | `Loiselle` |
| [Lynch 1988](https://pubmed.ncbi.nlm.nih.gov/3193879/) | moments-based | `Lynch` |
| [Lynch/Li 1993](https://pubmed.ncbi.nlm.nih.gov/8514326/) | moments-based | `LynchLi` |
| [Lynch & Ritland 1999](https://www.genetics.org/content/152/4/1753.short) | moments-based | `LynchRitland` |
| [Moran 1950](https://www.jstor.org/stable/2332142?origin=crossref&seq=1#metadata_info_tab_contents) | moments-based | `Moran` |
| [Queller & Goodnight 1989](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1989.tb04226.x) | moments-based | `QuellerGoodnight` |
| [Ritland 1996](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1996.tb02347.x) | moments-based | `Ritland` |


```julia
julia> cats = @nancycats ; 

julia> kin = kinship(cats, method = Moran)
237×237 Named Matrix{Float64}
A ╲ B │         N215          N216  …          N289          N290
──────┼──────────────────────────────────────────────────────────
N215  │ 8.13724e-316       1.62338  …       1.04589       1.15351
N216  │      1.62338       0.29485         0.957724        1.1637
N217  │     0.673577      0.587163         0.547427      0.709887
N218  │     0.896935       0.72942         0.919448      0.791255
⋮                  ⋮             ⋮  ⋱             ⋮             ⋮
N297  │     0.757915      0.858834          1.15432        1.2677
N281  │     0.686057      0.604236         0.942749       1.08762
N289  │      1.04589      0.957724              0.0         1.104
N290  │      1.15351        1.1637  …         1.104           0.0

julia> kinship(cats, method = Moran, iterations = 100)
27966×7 DataFrame
   Row │ sample1  sample2  Moran     bootmean  std       CI_lower      CI_upper 
       │ String   String   Float64   Float64   Float64   Float64       Float64  
───────┼────────────────────────────────────────────────────────────────────────
     1 │ N215     N216     1.62338   0.376626  0.27286    0.00274863   0.916719
     2 │ N215     N217     0.673577  0.202888  0.20094    0.00105976   0.59871
     3 │ N215     N218     0.896935  0.206272  0.232048   7.58373e-5   0.786113
     4 │ N215     N219     0.988931  0.236503  0.221345  -0.00053018   0.718204
   ⋮   │    ⋮        ⋮        ⋮         ⋮         ⋮           ⋮           ⋮
 27964 │ N281     N289     0.942749  0.220475  0.200358   0.001656     0.799307
 27965 │ N281     N290     1.08762   0.285053  0.289967   0.000299019  1.09343
 27966 │ N289     N290     1.104     0.277445  0.235519   0.00186445   0.858206
                                                              27959 rows omitted
```
----

### 📦 _kinship_noboot_nofreq
```julia
_kinship_noboot_nofreq(data::PopData, method::Function)
```
----

### 📦 _kinship_noboot_freq
```julia
_kinship_noboot_freq(data::PopData, method::Function)
```
----

### 📦 _kinship_boot_nofreq
```julia
_kinship_boot_nofreq(data::PopData, method::Function, iterations::Int, interval::Vector{Float64} = [0.025, 0.975])
```

### 📦 _kinship_boot_freq
```julia
_kinship_boot_freq(data::PopData, method::Function, iterations::Int, interval::Vector{Float64} = [0.025, 0.975])
```

### 🔵 kinshiptotable
```julia
kinshiptotable(kinshipresults::T, methd::Symbol) where T<:NamedMatrix
```
Converts the `NamedMatrix` result from the non-bootstrapped `kinship()` results into a `DataFrame`.
The second positonal argument (`methd`) is the name of the value column (default: `kinship`). For
better analysis workflow, it would be useful to specify the method for this column, to
keep track of which estimator was used (e.g., `Blouin`, `LynchLi`, etc.)
**Example**
```julia
julia> cats = @nancycats ; kin = kinship(cats, method = Moran) ;
julia> kinshiptotable(kin, :Moran)
22366×3 DataFrame
   Row │ sample1  sample2  Moran      
       │ String   String   Float64      
───────┼────────────────────────────────
     1 │ cc_001   cc_002    0.00688008
     2 │ cc_001   cc_003   -0.0286812
     3 │ cc_001   cc_005   -0.000749142
     4 │ cc_001   cc_007    0.0516361
     5 │ cc_001   cc_008    0.0261128
     6 │ cc_001   cc_009   -0.00187027
     7 │ cc_001   cc_010    0.0182852
   ⋮   │    ⋮        ⋮          ⋮
 22361 │ seg_028  seg_029  -0.0472928
 22362 │ seg_028  seg_030  -0.0172853
 22363 │ seg_028  seg_031  -0.00240921
 22364 │ seg_029  seg_030  -0.0278483
 22365 │ seg_029  seg_031   0.0297876
 22366 │ seg_030  seg_031  -0.0371295
                      22353 rows omitted
```