"""
    _fst_permutation(population_1::T, population_2::T) where T<:AbstractMatrix
Returns two matrices with rows (samples) shuffled between them. Respects the
number of rows of the original matrices (i.e. population sizes).
"""
function _fst_permutation(data::T, n1::Integer, n2::Integer) where T<:AbstractMatrix
    # total number of samples
    n_range = 1:(n1 +n2)
    # generate permutation indices for the first "population"
    perm_matrix = @inbounds data[shuffle(n_range), :]
    #partition shuffled matrix into original sample sizes
    new_pop_1, new_pop_2 = partitionarray(perm_matrix, [n1, n2])
    return new_pop_1, new_pop_2 
end

function _fst_permutation(data::T) where T<:AbstractMatrix
    @inbounds data[shuffle(1:size(data,1)), :]
end


function _permuted_Hudson(data::PopData, iterations::Int64)
    !isbiallelic(data) && throw(error("Data must be biallelic to use the Hudson estimator"))
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = data.metadata.populations
    n_loci = data.metadata.loci
    results = zeros(npops, npops)
    perm_vector = zeros(iterations-1)
    p = Progress(sum(1:(npops - 1)) * (iterations - 1), dt = 1, color = :blue, barglyphs = BarGlyphs("|=> |"), barlen = 30)
    @inbounds for i in 2:npops
        pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
        n_pop1 = size(pop1, 1)
        @inbounds for j in 1:(i-1)
            pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
            n_pop2 = size(pop2, 1)
            merged = vcat(pop1, pop2)
            fst_val = _hudson_fst(pop1,pop2)
            @inbounds @sync for iter in 1:iterations-1
                Base.Threads.@spawn begin
                    perm_p1, perm_p2 = _fst_permutation(merged, n_pop1, n_pop2)
                    @inbounds perm_vector[iter] = _hudson_fst(perm_p1, perm_p2)
                    pops_text = string(pops[i]) * " & " * string(pops[j])
                    ProgressMeter.next!(p; showvalues = [(:Populations, pops_text), (:Iteration, "$iter")])
                end
            end
            @inbounds results[i,j] = fst_val
            @inbounds results[j,i] = (sum(fst_val .<= perm_vector) + 1) / iterations
        end
    end
    println("Below diagonal: FST values | Above diagonal: P values")
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Hudson et al. 1992 (with p-values)")
end


function _permuted_Nei(data::PopData, iterations::Int64)
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = data.metadata.populations
    n_loci = data.metadata.loci
    results = zeros(npops, npops)
    perm_vector = zeros(iterations-1)
    p = Progress(sum(1:(npops - 1)) * (iterations - 1), dt = 1, color = :blue, barglyphs = BarGlyphs("|=> |"), barlen = 30)
    @inbounds for i in 2:npops
        pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
        n_pop1 = size(pop1, 1)
        @inbounds for j in 1:(i-1)
            pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
            n_pop2 = size(pop2, 1)
            merged = vcat(pop1, pop2)
            fst_val = _nei_fst(pop1,pop2)
            @inbounds @sync for iter in 1:iterations-1
                Base.Threads.@spawn begin
                    perm_p1, perm_p2 = _fst_permutation(merged, n_pop1, n_pop2)
                    @inbounds perm_vector[iter] = _nei_fst(perm_p1, perm_p2)
                    pops_text = string(pops[i]) * " & " * string(pops[j])
                    ProgressMeter.next!(p; showvalues = [(:Populations, pops_text), (:Iteration, "$iter")])
                end
            end
            @inbounds results[i,j] = fst_val
            @inbounds results[j,i] = (sum(fst_val .<= perm_vector) + 1) / iterations
        end
    end
    println("Below diagonal: FST values | Above diagonal: P values")
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Nei 1987 (with p-values)")
end


function _permuted_WeirCockerham(data::PopData, iterations::Int64)
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = data.metadata.populations
    n_loci = data.metadata.loci
    results = zeros(npops, npops)
    perm_vector = zeros(iterations-1)
    p = Progress(sum(1:(npops - 1)) * (iterations - 1), dt = 1, color = :blue, barglyphs = BarGlyphs("|=> |"), barlen = 30)
    @inbounds for i in 2:npops
        pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
        n_pop1 = size(pop1, 1)
        @inbounds for j in 1:(i-1)
            pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
            n_pop2 = size(pop2, 1)
            merged = vcat(pop1, pop2)
            fst_val = _weircockerham_fst(pop1,pop2)
            @inbounds @sync for iter in 1:iterations-1
                Base.Threads.@spawn begin
                    perm_p1, perm_p2 = _fst_permutation(merged, n_pop1, n_pop2)
                    @inbounds perm_vector[iter] = _weircockerham_fst(perm_p1, perm_p2)
                    pops_text = string(pops[i]) * " & " * string(pops[j])
                    ProgressMeter.next!(p; showvalues = [(:Populations, pops_text), (:Iteration, "$iter")])
                end
            end
            @inbounds results[i,j] = fst_val
            @inbounds results[j,i] = (sum(fst_val .<= perm_vector) + 1) / iterations
        end
    end
    println("Below diagonal: FST values | Above diagonal: P values")
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Weir & Cockerham (with p-values)")
end
