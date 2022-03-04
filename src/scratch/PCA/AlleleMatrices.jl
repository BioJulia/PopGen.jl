"""
    allele_matrix(data::PopData; by::String = "count", missings = "mean", scale = false, center = false)
Return a matrix of allele counts or frequencies per genotype where rows are samples
and columns are the occurence count or frequency of an allele for that locus in that sample.
Loci and alleles are sorted alphanumerically. Setting `scale` or `center` as `true` will
compute allele frequencies regardless of the `by` keyword.

### Keyword Arguments
- `by`: a `String` of `count` or `frequency` (default: `count`)
- `missings`: a `String` denoting how to handle missing values when outputting `frequency` (default: `mean`)
    - `"missing"`: fallback method to keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')

**Example**

```
julia> cats = @nancycats ;

julia> frq = allele_matrix(cats) ;  frq[1:5,1:6]
5×6 Matrix{Union{Missing, Int8}}:
  missing   missing   missing   missing   missing   missing
  missing   missing   missing   missing   missing   missing
 0         0         0         0         0         0
 0         0         0         0         0         0
 0         0         0         0         0         0

julia> frq = allele_matrix(cats, by = "frequency") ;  frq[1:5,1:6]
5×6 Matrix{Union{Missing, Float32}}:
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0

 julia> frq = allele_matrix(cats, by = "frequency", missings = "zero") ;  frq[1:5,1:6]
 5×6 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0

 julia> frq = allele_matrix(cats, missings = "mean", scale = true, center = true) ;  frq[1:5,1:6]
 5×6 Matrix{Float32}:
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
```
 """
function allele_matrix(data::PopData; by::String = "count", missings::String = "mean", scale::Bool = false, center::Bool = false)
    if occursin(lowercase(by), "counts") && !any([scale, center])
        return _count_matrix(data)
    elseif occursin(lowercase(by), "frequency") && !any([scale, center])
        return _freq_matrix(data, missings = missings)
    elseif any([scale, center])
        return _scaled_freq_matrix(data, missings = missings, scale = scale, center = center)
    else
        throw(ArgumentError("Choose from either \"count\" or \"frequency\" methods"))
    end
end


"""
    _scaled_freq_matrix(data::PopData; missings::String = "mean", scale::Bool = true, center::Bool = true)
Returns a Z-score scaled matrix of allele frequencies where rows are samples 
and columns are the frequency of an allele for that locus in that sample.
- `missings`: a `String` denoting how to handle missing values when outputting `frequency` (default: `mean`)
    - `"missing"`: fallback method to keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')
"""
function _scaled_freq_matrix(data::PopData; missings::String = "mean", scale::Bool = true, center::Bool = true)
    if lowercase(missings) ∈ ["0", "zero", "zeros", "mean"]
        freqs = _freq_matrix(data, missings = missings)
        mtx = standardize(ZScoreTransform, freqs, dims = 1, scale = scale, center = center)
        # replace almost-zero values caused by missing values with 0.0
        replace!(x ->  0 < x < (10^-9) ? 0.0 : x, mtx)
        return mtx
    else
        occursin(lowercase(missings), "missings") && throw(ArgumentError("\"missing\" method cannot be scaled"))
        !occursin(lowercase(missings), "missings") && throw(ArgumentError("use one of \"zero\" or \"mean\" for handling missing values"))
    end
end

"""
    _count_matrix(data::PopData)
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
This is differentiated from `allele_count_matrix` by preserving missing values as
`missing`
"""
function _count_matrix(data::PopData)
    geno_mtx = loci_matrix(data)
    n_samples, n_loci = size(data)
    keytype = eltype(data.loci.genotype |> eltype |> nonmissingtype)
    # create a vector of empty dicts
    all_alleles = [Dict{keytype, Vector{Union{Missing,Int8}}}() for i in 1:n_loci]
    # fill the dicts with keys for every allele for that locus, and the values are a vector of zeros
    for (i,j) in enumerate(eachcol(geno_mtx))
       [get!(all_alleles[i], k, zeros(Int8, n_samples)) for k in unique_alleles(j)]
    end
    # for every sample row
    @inbounds @sync for (idx,sample_row) in enumerate(eachrow(geno_mtx))
        #Base.Threads.@spawn begin
            # for each nonmissing genotype in that sample
            @inbounds for (loc_idx, geno) in enumerate(sample_row)
                # if the genotype is missing, fill in missing for all alleles of that locus
                if ismissing(geno) 
                    [all_alleles[loc_idx][i][idx] = missing for i in keys(all_alleles[loc_idx])] 
                    continue
                end
                # otherwise, if genotype is present
                @inbounds for allele in geno
                    # +1 to its counter
                    @inbounds all_alleles[loc_idx][allele][idx] += 1
                end
            end
        #end
    end
    # concatenate the alleles into a matrix and those into one big matrix
    return reduce(hcat, [reduce(hcat, values(sort(i))) for i in all_alleles])
end


"""
    _freq_matrix(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
"""
function _freq_matrix(data::PopData; missings::String = "zero")
    counts = _count_matrix(data)
    if lowercase(missings) ∈ ["0", "zero", "zeros"]
        # replace occurences of missing with 0
        @inbounds counts[findall(ismissing, counts)] .= 0
        # divide each row (sample) by the ploidy of that sample
        return counts ./ data.meta.ploidy
    elseif lowercase(missings) == "mean"
        # covert the matrix to Floats
        counts = Matrix{Union{Missing, Float32}}(counts)
        # get the mean for each allele for each locus
        # iterate over rows (samples)
        @inbounds @sync for (_sample,_counts) in enumerate(eachrow(counts))
            Base.Threads.@spawn begin
                # divide non-missing values in the row by the sample's ploidy and replace the original values
                @inbounds _counts[.!ismissing.(_counts)] .= _counts[.!ismissing.(_counts)] ./ data.meta.ploidy[_sample]
            end
        end
        #means = map(x -> mean(skipmissing(x)), eachcol(counts))
        # fill in missing values with the means
        @inbounds for (idx,_allele) in enumerate(eachcol(counts))
            @inbounds _allele[ismissing.(_allele)] .= mean(skipmissing(_allele))
        end
        return Matrix{Float32}(counts)
    else
        !occursin(lowercase(missings), "missings") && println("Incorrect missing data handling method, defaulting to \"missing\"\n")
        # convert the matrix to Floats
        counts = Matrix{Union{Missing, Float32}}(counts)
        # iterate over rows (samples)
        @inbounds @sync for (_sample,_counts) in enumerate(eachrow(counts))
            Base.Threads.@spawn begin
                # divide non-missing values in the row by the sample's ploidy and replace the original values
                @inbounds _counts[.!ismissing.(_counts)] .= _counts[.!ismissing.(_counts)] ./ data.meta.ploidy[_sample]
            end
        end
        return counts
    end
end


#=  this generates the locus_allele names
    #Base.Iterators.flatten(all_alleles) |> length
    #loc_names = loci(data)
    #map(zip(loc_names, all_alleles)) do (l,a)
    #    l * "_" .* string.(a)
=#
