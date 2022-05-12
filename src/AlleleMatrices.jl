"""
    _allelematrix(data::PopData; by::String = "count", missings = "mean", scale = false, center = false)
Return a matrix of allele counts or frequencies per genotype where rows are samples
and columns are the occurence count or frequency of an allele for that locus in that sample.
Loci and alleles are sorted alphanumerically. Setting `scale` or `center` as `true` will
compute allele frequencies regardless of the `by` keyword.

### Keyword Arguments
- `by`: a `String` of `count` or `frequency` (default: `frequency`)
- `missings`: a `String` denoting how to handle missing values when outputting `frequency` (default: `mean`)
    - `"missing"`: fallback method to keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')

**Example**

```
julia> cats = @nancycats ;

julia> frq = _allelematrix(cats) ;  frq[1:5,1:6]
5×6 Matrix{Union{Missing, Int8}}:
  missing   missing   missing   missing   missing   missing
  missing   missing   missing   missing   missing   missing
 0         0         0         0         0         0
 0         0         0         0         0         0
 0         0         0         0         0         0

julia> frq = _allelematrix(cats, by = "frequency") ;  frq[1:5,1:6]
5×6 Matrix{Union{Missing, Float32}}:
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0

 julia> frq = _allelematrix(cats, by = "frequency", missings = "zero") ;  frq[1:5,1:6]
 5×6 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0

 julia> frq = _allelematrix(cats, missings = "mean", scale = true, center = true) ;  frq[1:5,1:6]
 5×6 Matrix{Float32}:
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
```
 """
function _allelematrix(data::PopData; by::String = "frequency", missings::String = "mean", scale::Bool = false, center::Bool = false)
    if occursin(lowercase(by), "counts") && !any([scale, center])
        return _countmatrix(data)
    elseif occursin(lowercase(by), "frequency") 
        freqs = 
            if lowercase(missings) == "mean"
                _freqmatrix_mean(data)
            elseif lowercase(missings) == "zero"
                _freqmatrix_zero(data)
            elseif lowercase(missings) == "missing"
                _freqmatrix_missing(data)
            else
                throw(ArgumentError("use one of \"zero\" or \"mean\" for handling missing values"))
            end
        if any([scale, center])
            return _freqmatrix_scale(freqs, scale, center)
        else
            return freqs
        end
    else
        throw(ArgumentError("Choose from either \"count\" or \"frequency\" methods"))
    end
end

#=
function _countset(query::Genotype,reference::Vector{T}) where T<:Union{Int8, Int16}
    @inbounds [count(==(x), query, init = T(0)) for x in reference]
end

function _countset(query::Missing, reference)
    fill(eltype(reference)(-1), length(reference))
end
=#

function _setcounts(q, r)
    l = 0
    @inbounds for i in r
        l += length(i)
    end
    cnt = Vector{eltype(eltype(r))}(undef, l)
    idx = 0
    @inbounds for (i,j) in enumerate(r)
        @inbounds geno = q[i]
        @inbounds @simd for h in j 
            idx += 1
            @inbounds cnt[idx] = geno === missing ? -1 : count(==(h), geno)
        end
    end
    return cnt
end


"""
    _countmatrix(data::PopData)
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
Missing values are preserved as `-1`.
"""
function _countmatrix(data::PopData)
    gmtx = locimatrix(data)
    allalleles = Tuple(uniquealleles(i) for i in eachcol(gmtx))
    #return allalleles
    mapreduce(hcat, eachrow(gmtx)) do smple
        _setcounts(smple, allalleles)
        #[j for i in _countset.(smple, allalleles) for j in i]
    end |> permutedims
end


"""
    _freqmatrix_zero(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by zeros.
"""
function _freqmatrix_zero(data::PopData)
    # divide each row (sample) by the ploidy of that sample
    out = _countmatrix(data)
    replace!(out, -1 => 0)
    out ./ data.sampleinfo.ploidy
end

"""
    _freqmatrix_mean(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by the global mean allele frequency.
"""
function _freqmatrix_mean(data::PopData)
    counts = @inbounds _countmatrix(data) ./ data.sampleinfo.ploidy
    map(eachcol(counts)) do alcol
        colmean = mean([x for x in alcol if x >= 0])
        replace!(x -> x < 0 ? colmean : x, alcol)
    end
    return counts
end

"""
    _freqmatrix_missing(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are kept as `missing`.
"""
function _freqmatrix_missing(data::PopData)
    out = _countmatrix(data)
    replace!(out, -1 => missing)
    out ./ data.sampleinfo.ploidy
end


"""
    _freqmatrix_scale(freqs::Matrix{Float32}, scale::Bool = true, center::Bool = true)
Returns a Z-score scaled matrix of allele frequencies where rows are samples 
and columns are the frequency of an allele for that locus in that sample.
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')
"""
function _freqmatrix_scale(freqs::Matrix{Float64}, scale::Bool = true, center::Bool = true)
    mtx = standardize(ZScoreTransform, freqs, dims = 1, scale = scale, center = center)
    # replace almost-zero values caused by missing values with 0.0
    replace!(x ->  0 < x < (10^-9) ? 0.0 : x, mtx)
    return mtx
end
