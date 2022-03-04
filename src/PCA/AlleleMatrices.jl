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

"""
    _countmatrix(data::PopData)
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
This is differentiated from `allele_countmatrix` by preserving missing values as
`missing`
"""
function _countmatrix(data::PopData)
    geno_mtx = locimatrix(data)
    n_samples, n_loci = size(data)
    keytype = eltype(data.genodata.genotype |> eltype |> nonmissingtype)
    # create a vector of empty dicts
    @inbounds all_alleles = [Dict{keytype, Vector{Union{Missing,Int8}}}() for i in 1:n_loci]
    # fill the dicts with keys for every allele for that locus, and the values are a vector of zeros
    @inbounds for (i,j) in enumerate(eachcol(geno_mtx))
       [get!(all_alleles[i], k, zeros(Int8, n_samples)) for k in uniquealleles(j)]
    end
    # for every sample row
    @inbounds for (idx,sample_row) in enumerate(eachrow(geno_mtx))
        # for each nonmissing genotype in that sample
        @inbounds for (loc_idx, geno) in enumerate(sample_row)
            # if the genotype is missing, fill in missing for all alleles of that locus
            if ismissing(geno) 
                @inbounds [all_alleles[loc_idx][i][idx] = missing for i in keys(all_alleles[loc_idx])] 
                continue
            else
                # otherwise, if genotype is present
                @inbounds @simd for allele in geno
                    # +1 to its counter
                    @inbounds all_alleles[loc_idx][allele][idx] += Int8(1)
                end
            end
        end
    end
    # concatenate the alleles into a matrix and those into one big matrix
    return reduce(hcat, [reduce(hcat, values(sort(i))) for i in all_alleles])
end


"""
    _freqmatrix_zero(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by zeros.
"""
function _freqmatrix_zero(data::PopData)
    counts = _countmatrix(data)
    # replace occurences of missing with 0
    @inbounds counts[findall(ismissing, counts)] .= 0
    # divide each row (sample) by the ploidy of that sample
    return counts ./ data.sampleinfo.ploidy
end

"""
    _freqmatrix_mean(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by the global mean allele frequency.
"""
function _freqmatrix_mean(data::PopData)
    counts = Matrix{Union{Missing, Float64}}(_countmatrix(data))
    # get the mean for each allele for each locus
    # iterate over rows (samples)
    @inbounds for (_sample,_counts) in enumerate(eachrow(counts))
        countidx = findall(!ismissing, _counts)
        # divide non-missing values in the row by the sample's ploidy and replace the original values
        #@inbounds _counts[countidx] /= data.sampleinfo.ploidy[_sample]
        @inbounds _counts[countidx] .= _counts[countidx] ./ data.sampleinfo.ploidy[_sample]
    end
    # fill in missing values with the means
    @inbounds for _allele in eachcol(counts)
        @inbounds _allele[findall(ismissing, _allele)] .= mean(skipmissing(_allele))
    end
    return Matrix{Float64}(counts)
end

"""
    _freqmatrix_missing(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are kept as `missing`.
"""
function _freqmatrix_missing(data::PopData)
    counts = Matrix{Union{Missing, Float64}}(_countmatrix(data))
    # iterate over rows (samples)
    @inbounds for (_sample,_counts) in enumerate(eachrow(counts))
        countidx = findall(!ismissing, _counts)
        # divide non-missing values in the row by the sample's ploidy and replace the original values
        @inbounds _counts[countidx] .= _counts[countidx] ./ data.sampleinfo.ploidy[_sample]
    end
    return counts
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
