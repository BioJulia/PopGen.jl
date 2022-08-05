"""
    missingdata(data::PopData; by::Union{String, Symbol} = "sample")
Get missing genotype information in a `PopData`. Specify a mode of operation
to return a DataFrame corresponding with that missing information.

#### Modes
- "sample" - returns a count and percent and list of missing loci per individual (default)
- "population" - returns a count and percent of missing genotypes per population
- "locus" - returns a count and percent of missing genotypes per locus
- "locusxpopulation" - returns a count and percent of missing genotypes per locus per population

### Example:
```
missingdata(@gulfsharks, by = "pop")
```
"""
function missingdata(data::PopData; by::Union{String, Symbol} = "sample")
    if string(by) ∈ ["sample", "population", "locus", "locusxpopulation"]
        _missingdata(data, Val(Symbol(by)))
    else
        throw(ArgumentError("Mode \"$by\" not recognized. Please specify one of: sample, population, locus, or full"))
    end
end

function _missingdata(data::PopData, ::Val{:sample})
    df = DataFrames.combine(DataFrames.groupby(data.genodata, :name), :genotype => (i -> count(ismissing, i)) => :missing)
    df[:, :percent] = round.(df.missing ./ data.metadata.loci, digits = 3)
    df
end

function _missingdata(data::PopData, ::Val{:population})
    df = DataFrames.combine(DataFrames.groupby(data.genodata, :population), :genotype => (i -> count(ismissing, i)) => :missing, :genotype => length => :n)
    df[:, :percent] = round.(df.missing ./ df.n, digits = 3)
    select!(df, 1, 2, 4)
    df
end

function _missingdata(data::PopData, ::Val{:locus})
    df = DataFrames.combine(DataFrames.groupby(data.genodata, :locus), :genotype => (i -> count(ismissing, i)) => :missing)
    df[:, :percent] = round.(df.missing ./ data.metadata.samples, digits = 3)
    df
end

function _missingdata(data::PopData, ::Val{:locusxpopulation})
    df = DataFrames.combine(DataFrames.groupby(data.genodata, [:locus, :population]), :genotype => (i -> count(ismissing, i)) => :missing, :genotype => length => :n)
    df[:, :percent] = round.(df.missing ./ df.n, digits = 3)
    select!(df, 1, 2,3,5)
    df
end


function _pwiseidenticalhelper(x::T,y::U)::Float64 where T<:AbstractArray where U<:AbstractArray
    μsum = 0
    l = 0
    for i in 1:length(x)
        eq = x[i] == y[i]
        μsum += eq === missing ? 0 : eq
        l += eq === missing ? 0 : 1
    end
    return μsum / l
end

"""
    pairwiseidentical(data::PopData)
Return a pairwise matrix of the percent of identical genotypes at each locus between all pairs of individuals.

### Example:
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

"""
function pairwiseidentical(data::PopData)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    result = NamedArray(zeros(Float64, n, n))
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    @inbounds for i in 1:n-1
        @inbounds v1 = view(locmtx,i,:)
        @inbounds for j in i+1:n
            @inbounds v2 = view(locmtx,j,:)
            res = _pwiseidenticalhelper(v1, v2)
            @inbounds result[i,j] = res
            @inbounds result[j,i] = res
        end
    end
    # fill in diagonal
    for i in 1:n
        @inbounds result[i,i] = 1.0
    end
    return result
end

"""
    pairwiseidentical(data::PopData, sample_names::Vector{String})
Return a pairwise matrix of the percent of identical genotypes at 
each nonmissing locus between all pairs of provided `sample_names`.

### Example:
```
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
"""
function pairwiseidentical(data::PopData, sample_names::Vector{T}) where T<:AbstractString
    all_samples = samplenames(data)
    missingsamples = setdiff(sample_names, all_samples)
    if !isempty(missingsamples)
        throw(ArgumentError("Samples not found in the PopData:\n  " * join(missingsamples, "\n  ")))
    end
    sampidx = indexin(sample_names, all_samples)
    locmtx = locimatrix(data)
    n = length(sampidx)
    result = NamedArray(zeros(Float64, n, n))
    s_names = T == String ? sample_names : string.(sample_names)
    setnames!(result, s_names,1)
    setnames!(result, s_names,2)
    @inbounds for i in 1:n-1
        @inbounds v1 = view(locmtx,sampidx[i],:)
        @inbounds for j in i+1:n
            @inbounds v2 = view(locmtx,sampidx[j],:)
            res = _pwiseidenticalhelper(v1, v2)
            @inbounds result[i,j] = res
            @inbounds result[j,i] = res
        end
    end
    # fill in diagonal
    for i in 1:n
        @inbounds result[i,i] = 1.0
    end
    return result
end

"""
    genofreqtable(data::PopData; by::Union{Symbol,String} = "global")
Return a table of the observed `global` (default) or `population` genotype frequencies in a PopData object.

### Example:
```
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
"""
function genofreqtable(data::PopData; by::Union{Symbol,String} = "global")
    strby = lowercase(string(by))
    if strby == "global"
        grp = groupby(dropmissing(data.genodata, :genotype), [:locus, :genotype])
        counts = DataFrames.combine(grp, nrow => :count)
        counts = DataFrames.combine(
            groupby(counts, :locus), :genotype, :count,
            :count => (x -> x / sum(x)) => :frequency
        )
    elseif strby == "population"
        grp = groupby(dropmissing(data.genodata, :genotype), [:locus, :population, :genotype])
        counts = DataFrames.combine(grp, nrow => :count)
        counts = DataFrames.combine(
            groupby(counts, [:locus, :population]), :genotype, :count,
            :count => (x -> x / sum(x)) => :frequency
        )
    else
        throw(ArgumentError("Please use by = \"global\" (default) or \"population\""))
    end
end


"""
    allelefreqtable(data::PopData; by::Union{Symbol,String} = "global")
Return a table of the observed `global` (default) or `population` allele frequencies in a PopData object.

### Example:
```
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
"""
function allelefreqtable(data::PopData; by::Union{Symbol,String} = "global")
    strby = lowercase(string(by))
    if strby == "global"
        grp = groupby(dropmissing(data.genodata, :genotype), :locus)
        _alleles = DataFrames.combine(grp, :genotype => alleles => :allele, ungroup = false)
        counts = DataFrames.combine(
            _alleles,
            :allele,
            nrow => :n,
            ungroup = false
        )
        counts = DataFrames.combine(
            groupby(DataFrame(counts), [:locus, :allele]),
            nrow => :count,
            [:n, :allele] => ((n,al) -> length(al)/first(n)) => :frequency
        )
    elseif strby == "population"
        grp = groupby(dropmissing(data.genodata, :genotype), [:locus, :population])
        _alleles = DataFrames.combine(grp, :genotype => alleles => :allele, ungroup = false)
        counts = DataFrames.combine(
            _alleles,
            :allele,
            nrow => :n,
            ungroup = false
        )
        counts = DataFrames.combine(
            groupby(DataFrame(counts), [:locus, :population, :allele]),
            nrow => :count,
            [:n, :allele] => ((n,al) -> length(al)/first(n)) => :frequency
        )
    else
        throw(ArgumentError("Please use by = \"global\" (default) or \"population\""))
    end
end
