
export pairwise_identical, missing_data, geno_freqtable, allele_freqtable

#TODO update docs (API/dataexploration) and tests
"""
    missing_data(data::PopData; by::String = "sample")
Get missing genotype information in a `PopData`. Specify a mode of operation
to return a DataFrame corresponding with that missing information.

#### Modes
- "sample" - returns a count and list of missing loci per individual (default)
- "pop" - returns a count of missing genotypes per population
- "locus" - returns a count of missing genotypes per locus
- "full" - returns a count of missing genotypes per locus per population

### Example:
```
missing_data(@gulfsharks, by = "pop")
```
"""
@inline function missing_data(data::PopData; by::String = "sample")
    if by ∈ ["sample", "individual"]
        DataFrames.combine(
            DataFrames.groupby(data.loci, :name),
            :genotype => (i -> count(ismissing, i)) => :missing
        )
    elseif by ∈ ["pop", "population"]
                DataFrames.combine(
            DataFrames.groupby(data.loci, :population),
            :genotype => (i -> count(ismissing, i)) => :missing
        )
    elseif by ∈ ["locus", "loci"]
        DataFrames.combine(
            DataFrames.groupby(data.loci, :locus),
            :genotype => (i -> count(ismissing, i)) => :missing
        )

    elseif by ∈ ["detailed", "full"]
        DataFrames.combine(
            DataFrames.groupby(data.loci, [:locus, :population]),
            :genotype => (i -> count(ismissing, i)) => :missing
        )
    else
        @error "Mode \"$by\" not recognized. Please specify one of: sample, pop, locus, or full"
        missing_data(data)
    end
end

"""
    pairwise_identical(data::PopData)
Return a table of the percent of identical genotypes at each locus between all pairs of individuals.

### Example:
```
julia> cats = @nancycats ;

julia> pairwise_identical(cats)
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

"""
function pairwise_identical(data::PopData)
    sample_names = collect(samples(data))
    pairwise_identical(data, sample_names)
end

"""
    pairwise_identical(data::PopData, sample_names::Vector{String})
Return a table of the percent of identical genotypes at each locus
between all pairs of provided `sample_names`.

### Example:
```
julia> cats = @nancycats ;

julia> interesting_cats = samples(cats)[1:5]
5-element Array{String,1}:
 "N215"
 "N216"
 "N217"
 "N218"
 "N219"

julia> pairwise_identical(cats, interesting_cats)
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
"""
function pairwise_identical(data::PopData, sample_names::Vector{String})
    errs = ""
    all_samples = samples(data)
    if sample_names != all_samples
        [errs *= "\n  $i" for i in sample_names if i ∉ all_samples]
        errs != "" && error("Samples not found in the PopData: " * errs)
    end
    sample_pairs = pairwise_pairs(sample_names)
    n = length(sample_pairs)
    perc_ident_vec = Vector{Float64}(undef, n)
    n_vec = Vector{Int}(undef, n)
    popdata_idx = groupby(data.loci, :name)
    idx = 0
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    @inbounds @sync for i in 1:length(sample_pairs)
        Base.Threads.@spawn begin
            @inbounds geno_1 = popdata_idx[(sample_pairs[i][1],)].genotype
            @inbounds geno_2 = popdata_idx[(sample_pairs[i][2],)].genotype
            len_1 = nonmissing(geno_1)
            len_2 = nonmissing(geno_2)
            shared_geno = minimum([len_1, len_2])
                shared_geno = minimum([len_1, len_2])
                shared = sum(skipmissing(geno_1 .== geno_2))
                @inbounds perc_ident_vec[i] = round(shared/shared_geno, digits = 2)
                @inbounds n_vec[i] = shared_geno
            next!(p)
        end
    end
    DataFrame(:sample_1 => map(i -> i[1], sample_pairs), :sample_2 => map(i -> i[2], sample_pairs), :identical => perc_ident_vec, :n => n_vec)
end

"""
    geno_freqtable(data::PopData; by::String = "global")
Return a table of the observed `global` (default) or `population` genotype frequencies in a PopData object.

### Example:
```
julia> cats = @nancycats ;

julia> geno_freqtable(cats)

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

julia> geno_freqtable(cats, by = "population")
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
function geno_freqtable(data::PopData; by::String = "global")
    if by == "global"
        grp = groupby(dropmissing(data.loci, :genotype), [:locus, :genotype])
        counts = DataFrames.combine(grp, nrow => :count)
        counts = DataFrames.combine(
            groupby(counts, :locus), :genotype, :count,
            :count => (x -> x / sum(x)) => :frequency
        )
    elseif by in ["local", "population"]
        grp = groupby(dropmissing(data.loci, :genotype), [:locus, :population, :genotype])
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
    allele_freqtable(data::PopData; by::String = "global")
Return a table of the observed `global` (default) or `population` allele frequencies in a PopData object.

### Example:
```
julia> cats = @nancycats ;

julia> allele_freqtable(cats)
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

julia> allele_freqtable(cats, by = "population")
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
function allele_freqtable(data::PopData; by::String = "global")
    if by == "global"
        grp = groupby(dropmissing(data.loci, :genotype), :locus)
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
    elseif by in ["local", "population"]
        grp = groupby(dropmissing(data.loci, :genotype), [:locus, :population])
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
