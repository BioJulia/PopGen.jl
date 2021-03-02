export pairwise_fst

# generate a wrapper struct so we can nicely print the results
struct PairwiseFST
    results::DataFrame
    method::String
end

# pretty-printing of FST results
function Base.show(io::IO, data::PairwiseFST)
    show(
        io,
        data.results,
        show_row_number=false,
        rowlabel = Symbol(" "),
        eltypes = false,
        row_names = names(data.results),
        title = "Pairwise FST: " * data.method
    )
end


function FST_global(data::AbstractDataFrame)
    # observed/expected het per locus per pop
    het_df = DataFrames.combine(
        groupby(data, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => hetero_e => :het_exp,
        :genotype => allele_freq => :alleles
    )
    # collapse down to retrieve averages and counts
    n_df = DataFrames.combine(
        groupby(het_df, :locus),
        :n => count_nonzeros => :count,
        :n => (n -> count_nonzeros(n) / reciprocal_sum(n)) => :mn,
        [:het_obs, :het_exp, :n] => ((o,e,n) -> mean(skipmissing(gene_diversity_nei87.(e, o, count_nonzeros(n) / reciprocal_sum(n))))) => :HS,
        :het_obs => (o -> mean(skipmissing(o)))=> :Het_obs,
        :alleles => (alleles ->  sum(values(avg_allele_freq(alleles, 2))))=> :avg_freq
        )

    Ht = @inbounds 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    DST = @inbounds Ht .- n_df.HS
    DST′ = @inbounds n_df.count ./ (n_df.count .- 1) .* DST
    HT′ = @inbounds n_df.HS .+ DST′

    round(safemean(DST′) / safemean(HT′), digits = 4)
end


function _pairwise_Nei(data::PopData)
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    results = zeros(npops, npops)
    for i in 2:npops
        for j in 1:(i-1)
            results[i,j] = FST_global(DataFrame(idx_pdata[[j,i]]))
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Nei 1987")
end

function _pairwise_WeirCockerham(data::PopData)
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    results = zeros(npops, npops)
    for i in 2:npops
        for j in 1:(i-1)
            results[i,j] = weircockerham1984(DataFrame(idx_pdata[[j,i]]))
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Weir-Cockerham")
end

#TODO add methods and complete docstring
"""
    pairwise_fst(data::PopData; method::String)
Calculate pairwise FST between populations in a `PopData` object.
#### Methods:
- `"Nei87"`: Nei (1987) method
- `"WC84"` : Weir-Cockerham (1984) method (default)
"""
function pairwise_fst(data::PopData; method::String = "WC")
    occursin("nei", lowercase(method)) && return _pairwise_Nei(data)
    occursin("wc", lowercase(method)) && return _pairwise_WeirCockerham(data)
end


function FST_global(data::PopData)
    # observed/expected het per locus per pop
    het_df = DataFrames.combine(
        groupby(data, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => (i -> hetero_e(i)) => :het_exp,
        :genotype => allele_freq => :alleles
    )
    # collapse down to retrieve averages and counts
    n_df = DataFrames.combine(
        groupby(het_df, :locus),
        :n => (n -> sum(.!iszero.(n)))=> :count,
        :n => (n -> sum(.!iszero.(n)) ./ sum(reciprocal.(n))) => :mn,
        [:het_obs, :het_exp, :n] => ((o,e,n) -> mean(skipmissing(gene_diversity_nei87.(e,o,sum(.!iszero.(n)) ./ sum(reciprocal.(n)))))) => :HS,
        :het_obs => (o -> mean(skipmissing(o)))=> :Het_obs,
        :alleles => (alleles ->  sum(values(avg_allele_freq(alleles, 2))))=> :avg_freq
        )

    Ht = @inbounds 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    DST = @inbounds Ht .- n_df.HS
    DST′ = @inbounds n_df.count ./ (n_df.count .- 1) .* DST
    HT′ = @inbounds n_df.HS .+ DST′

    HT = safemean(HT)
    DST = safemean(DST)
    DST′ = safemean(DST′)
    HT′ = safemean(HT′)

    return Dict{Symbol, Float64}(
        :FST => DST / HT,
        :FST′ => DST′ / HT′
    )
end

function FST_locus(data::PopData)
    # observed/expected het per locus per pop
    het_df = DataFrames.combine(
        groupby(data.loci, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => (i -> hetero_e(i)) => :het_exp,
        :genotype => allele_freq => :alleles
    )
    # collapse down to retrieve averages and counts
    n_df = DataFrames.combine(
        [:het_obs, :het_exp, :n, :alleles] => (o,e,n,alleles) -> (
            count = sum(.!iszero.(n)),
            mn = sum(.!iszero.(n)) ./ sum(reciprocal.(n)),
            HS = mean(skipmissing(gene_diversity_nei87.(e,o,sum(.!iszero.(n)) ./ sum(reciprocal.(n))))),
            Het_obs = mean(skipmissing(o)),
            avg_freq = sum(values(avg_allele_freq(alleles)).^2)
            ),
        groupby(het_df, :locus)
    )

    HT = 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    DST = HT .- n_df.HS
    DST′ = n_df.count ./ (n_df.count .- 1) .* DST
    HT′ = n_df.HS .+ DST′

    return Dict{Symbol, Vector{Float64}}(
        :FST => DST ./ HT,
        :FST′ => DST′ ./ HT′
        )
end



function f_stat_p(data::PopData; nperm::Int = 1000)
    observed = FST_global(data)
    popnames = data.meta.population    
    perm_dict = Dict{Symbol,Vector{Float64}}(
        :FST => Vector{Float64}(undef, nperm-1),
        :FST′ => Vector{Float64}(undef, nperm-1)
        )
    p = Progress(nperm-1, dt = 1, color = :blue)
    @inbounds @sync for n in 1:nperm-1
        Base.Threads.@spawn begin
            tmp = FST_global(permute_samples!(copy(data.loci), popnames))
            perm_dict[:FST][n] = tmp[:FST]
            perm_dict[:FST′][n] = tmp[:FST′]
            next!(p)
        end
    end
    @inbounds for (k,v) in perm_dict
        observed[k] = (1 + sum(observed[k] .<= v))/nperm
    end

    return observed
end

