function allele_matrix(data::PopData)
    geno_mtx = loci_matrix(data)
    n_samples, n_loci = size(data)
    # create a vector of empty dicts
    all_alleles = [Dict{eltype(x.loci.genotype |> eltype |> nonmissingtype), Vector{Int8}}() for i in 1:n_loci]
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
    #Base.Iterators.flatten(all_alleles) |> length
    #loc_names = loci(data)
    #map(zip(loc_names, all_alleles)) do (l,a)
    #    l * "_" .* string.(a)
    #end        
end