function _permuted_WeirCockerham(data::PopData, iterations::Int64)
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    n_loci = size(data)[2]
    results = zeros(npops, npops)
    perm_vector = zeros(iterations-1)
    p = Progress(sum(1:(npops - 1)) * (iterations - 1), dt = 1, color = :blue, barglyphs = BarGlyphs("|=> |"), barlen = 30)
    @inbounds for i in 2:npops
        @inbounds for j in 1:(i-1)
            pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
            pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
            wc_val = weircockerham_fst(pop1,pop2)
            @inbounds @sync for iter in 1:iterations-1
                Base.Threads.@spawn begin
                    perm_p1, perm_p2 = wc_permutation(pop1, pop2)
                    perm_vector[iter] = weircockerham_fst(perm_p1, perm_p2)
                    pops_text = string(pops[i]) * " & " * string(pops[j])
                    ProgressMeter.next!(p; showvalues = [(:Populations, pops_text), (:Iteration, "$iter")])
                end
            end
            @inbounds results[i,j] = wc_val
            @inbounds results[j,i] = (sum(wc_val .<= perm_vector) + 1) / iterations
        end
    end
    println("Below diagonal: FST values | Above diagonal: P values")
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Weir & Cockerham")
end

function wc_permutation(population_1::T, population_2::T) where T<:AbstractMatrix
    merged = vcat(population_1, population_2)
    # get n for each population
    p1_size = size(population_1, 1)
    p2_size = size(population_2, 1)
    # total number of samples
    n_samples = p1_size + p2_size
    # generate permutation indices for the first "population"
    perm1 = sample(Xoroshiro128Star(), 1:n_samples, p1_size, replace = false)
    # generate permutation indices for the second "population"
    perm2 = @inbounds (1:n_samples)[Not(perm1)]
    # index the merged matrix with the permuted indices
    new_pop_1 = @inbounds @view merged[perm1,:]
    new_pop_2 = @inbounds @view merged[perm2,:]
    return new_pop_1, new_pop_2 
end