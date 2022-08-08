"""
    _fst_permute(population_1::T, population_2::T) where T<:AbstractMatrix
Returns two matrices with rows (samples) shuffled between them. Respects the
number of rows of the original matrices (i.e. population sizes).
"""
function _fst_permute(data::T, n1::Integer, n2::Integer) where T<:AbstractMatrix
    perm_matrix = @inbounds view(data, shuffle(1:(n1 +n2)), :)
    #partition shuffled matrix into original sample sizes
    new_pop_1, new_pop_2 = partitionarray(perm_matrix, [n1, n2])
    return new_pop_1, new_pop_2 
end

#= deprecated? #
function _fst_permutation(data::T) where T<:AbstractMatrix
    @inbounds data[shuffle(1:size(data,1)), :]
end
=#
"""
    _fst_permution(data::PopData, method::Function, iterations::Int64)
Returns a `PairwiseFST` object containing a dataframe of Pairwise FST calculations. The contained
dataframe has FST values below the diagonal and P values above it. This method is used internally
and wrapped by the public API provided in `pairwisefst()`.
"""
function _fst_permutation(data::PopData, method::Function, iterations::Int64)
    method == AMOVA && return _amovafst_permutation(data, iterations)
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = data.metadata.populations
    n_loci = data.metadata.loci
    results = zeros(Float64, npops, npops)
    #perm_vector = zeros(Int64, iterations-1)
    pbar = ProgressBar(;refresh_rate=90, transient = true)
    job = addjob!(pbar; description= "$(string(method)) FST: ", N = Int64((npops * (npops-1))/2))
    start!(pbar)
    @inbounds for i in 2:npops
        pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
        n_pop1 = size(pop1, 1)
        @inbounds for j in 1:(i-1)
            pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
            n_pop2 = size(pop2, 1)
            merged = vcat(pop1, pop2)
            fst_val = method(pop1,pop2)
            pval = 0
            @inbounds @sync for iter in 1:iterations-1
                Base.Threads.@spawn begin
                    @inbounds perm_p1, perm_p2 = _fst_permute(merged, n_pop1, n_pop2)
                    pval += fst_val <= method(perm_p1, perm_p2)
                end
            end
            @inbounds results[i,j] = fst_val
            @inbounds results[j,i] = (pval + 1) / iterations 
            update!(job)
        end
    end
    stop!(pbar)
    println("Below diagonal: FST values | Above diagonal: P values")
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "$method (with p-values)")
end

function _amovafst_permutation(data::PopData, iterations::Int64)
    distmtx = pairwise(SqEuclidean(), matrix(data), dims=1)
    grpcol = data.metadata.sampleinfo.population
    levels = unique(grpcol)
    npops = data.metadata.populations
    # create vector of vectors containing indices for each unique group
    groupidx = map(x -> findall(==(x),grpcol), levels)
    #permNs = length.(groupidx)
    results = zeros(Float64, npops, npops)
    df_among = 1
    pbar = ProgressBar(;refresh_rate=90, transient = true)
    job = addjob!(pbar; description= "AMOVA-based FST: ", N = Int64((npops * (npops-1))/2))
    start!(pbar)
    for i in 2:npops
        pop1 = groupidx[i]
        n1 = length(pop1)
        dw1 = reduce(+, view(distmtx, pop1, pop1))
        SS_within1 = dw1 / 2n1
        for j in 1:(i-1)
            pop2 = groupidx[j]
            n2 = length(pop2)
            N = n1 + n2
            dw2 = reduce(+, view(distmtx, pop2, pop2))
            da = reduce(+, view(distmtx, pop1, pop2))
            df_within = N - 2
            # possibly needs to be include 2nd pop?
            SS_among1 = ((da + dw1) / 2N) - (dw1 / 2n1)
            SS_among2 = ((da + dw2) / 2N) - (dw2 / 2n2)
            SS_among = SS_among1 + SS_among2
            σ²_within = (SS_within1 + (dw2 / 2n2)) / df_within
            n_c = (N - ((n1^2 + n2^2)/N)) #/ df_among
            #return n_c
            #n_c = (N - (reduce(+, length.(groupidx[[i,j]]).^2) / N)) / df_among
            σ²_among = ((SS_among / df_among) - σ²_within) / n_c
            FST = σ²_among / (σ²_among + σ²_within)
            #@inbounds results[i,j] =  results[j,i] = FST
            @inbounds results[i,j] = FST
            
            pval = 0
            @inbounds @sync for iter in 1:iterations-1
                Base.Threads.@spawn begin
                    p1, p2 = partitionarray(shuffle(vcat(pop1, pop2)), [n1, n2])
                    i_dw1 = reduce(+, view(distmtx, p1, p1))
                    i_SS_within1 = i_dw1 / 2n1
                    i_dw2 = reduce(+, view(distmtx, p2, p2))
                    i_da = reduce(+, view(distmtx, p1, p2))
                    i_SS_among1 = ((i_da + i_dw1) / 2N) - (i_dw1 / 2n1)
                    i_SS_among2 = ((i_da + i_dw2) / 2N) - (i_dw2 / 2n2)
                    i_SS_among = i_SS_among1 + i_SS_among2
                    i_σ²_within = (i_SS_within1 + (i_dw2 / 2n2)) / df_within
                    i_n_c = (N - ((n1^2 + n2^2)/N)) #/ df_among
                    i_σ²_among = ((i_SS_among / df_among) - i_σ²_within) / i_n_c
                    i_FST = i_σ²_among / (i_σ²_among + i_σ²_within)
                    pval += FST <= i_FST
                end
            end
            @inbounds results[j,i] = (pval + 1) / iterations
            update!(job)
        end

    end
    return PairwiseFST(DataFrame(results, Symbol.(levels)), "AMOVA-based (with p-values)")
end