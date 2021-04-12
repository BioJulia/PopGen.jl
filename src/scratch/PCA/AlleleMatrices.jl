"""
    allele_matrix(data::PopData; by::String = "count")
Create a matrix of allele counts or frequencies per genotype where rows are samples
and columns are the occurence count or frequency of an allele for that locus in that sample.
Loci and alleles are sorted alphanumerically.

**Example**

```
julia> cats = @nancycats ;

julia> allele_matrix(cats)
237×108 Matrix{Int8}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  …  0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  0  1     0  0  0  0  2  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  1  0  0  0  0     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  0  1  …  0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  2  0  0  0  0     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  0  1     0  0  0  1  0  1  0  0  0  0  0  0
 ⋮              ⋮              ⋮        ⋱              ⋮              ⋮     
 0  0  0  0  0  0  0  1  0  0  0  1  0  …  0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0  0  1  0     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  1  0  0  0  1  0     0  0  0  1  0  0  0  0  0  1  0  0
 0  0  0  0  0  0  0  1  0  0  0  0  1     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  1  0     0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  1  0  0  1  …  0  0  0  2  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  1  0  0  1  0     0  0  0  2  0  0  0  0  0  0  0  0

julia> allele_matrix(cats, by = "frequency")
 237×108 Matrix{Float32}:
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.5  0.0  0.0  0.0  0.0  0.0  0.0
  ⋮                        ⋮              ⋱                      ⋮         
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5     0.0  0.0  0.0  0.0  0.5  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.5     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
 """
function allele_matrix(data::PopData; by::String = "count")
    if occursin(lowercase(by), "counts")
        _count_matrix(data)
    elseif occursin(lowercase(by), "frequency")
        _freq_matrix(data)
    else
        throw(ArgumentError("Choose from either \"count\" or \"frequency\" methods"))
    end
end


"""
    _count_matrix(data::PopData)
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
"""
function _count_matrix(data::PopData)
    geno_mtx = loci_matrix(data)
    n_samples, n_loci = size(data)
    keytype = eltype(x.loci.genotype |> eltype |> nonmissingtype)
    # create a vector of empty dicts
    all_alleles = [Dict{keytype, Vector{Int8}}() for i in 1:n_loci]
    # fill the dicts with keys for every allele for that locus, and the values are a vector of zeros
    for (i,j) in enumerate(eachcol(geno_mtx))
       [get!(all_alleles[i], k, zeros(Int8, n_samples)) for k in unique_alleles(j)]
    end
    # for every sample row
    @inbounds @sync for (idx,sample_row) in enumerate(eachrow(geno_mtx))
        Base.Threads.@spawn begin
            # for each nonmissing genotype in that sample
            @inbounds for locus in eachindex(skipmissing(sample_row))
                # for each allele in the nonmising genotype
                @inbounds for allele in sample_row[locus]
                    # +1 to its counter
                    @inbounds all_alleles[locus][allele][idx] += 1
                end
            end
        end
    end
    # concatenate the alleles into a matrix and those into one big matrix
    return reduce(hcat, [reduce(hcat, values(sort(i))) for i in all_alleles])
end


"""
    _freq_matrix(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
"""
function _freq_matrix(data::PopData)
    _count_matrix(data) ./ data.meta.ploidy
end


#=  this generates the locus_allele names
    #Base.Iterators.flatten(all_alleles) |> length
    #loc_names = loci(data)
    #map(zip(loc_names, all_alleles)) do (l,a)
    #    l * "_" .* string.(a)
=#