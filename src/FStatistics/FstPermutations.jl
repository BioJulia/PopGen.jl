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
            pbar.update!(job)
        end
    end
    stop!(pbar)
    println("Below diagonal: FST values | Above diagonal: P values")
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "$method estimator (with p-values)")
end
