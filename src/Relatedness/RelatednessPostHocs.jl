export relatedness_posthoc

function sig_within(data::PopData, results::DataFrame, population::String, iterations::Int = 20000)
    # add extra columns of population names
    pop_names = select(data.meta, :name => :sample_1 ,:population => :pop_1)
    tmp_res = innerjoin(results, pop_names, on = :sample_1)
    select!(pop_names, :sample_1 => :sample_2,:pop_1 => :pop_2)
    tmp_res = innerjoin(tmp_res, pop_names, on = :sample_2)
    tmp_res.groups = collect(zip(tmp_res.pop_1, tmp_res.pop_2))
    tmp_res.group = tmp_res.pop_1 .* tmp_res.pop_2
    select!(tmp_res, Not([:sample_1, :sample_2, :n_loci, :pop_1, :pop_2]))
    size(tmp_res,1) == 0 && error("Samples in the relatedness results not found in the provided PopData.")
    estimators = Symbol.(names(tmp_res)[begin:end-2])

    # extract the values for within-population estimates
    within_coeff = [tmp_res[tmp_res.group .== population^2, i] for i in estimators]
    # get the mean of those values
    mean_within = mean.(within_coeff)
    # extract the values of the between-population estimates
    among_coeff = [tmp_res[(population .∈ tmp_res.groups) .& (tmp_res.group .!= population^2), i] for i in estimators]
    # get the mean of those values
    mean_among = mean.(among_coeff)
    # find the difference
    differ = mean_within - mean_among

    # setup bootstrap
    n_within = length(within_coeff[1])
    n_among = length(among_coeff[1])
    n_tot = n_within + n_among
    n_range = collect(1:n_tot)
    n_iter = iterations - 1
    
    bootstrapped = [Vector{Float64}(undef, n_iter) for i in 1:length(estimators)]
    p = Progress(iterations * length(estimators), dt = 1, color = :blue)
    
    for i in 1:length(estimators)
        est = estimators[i]
        iter_count = 0 
        @sync @inbounds for j in 1:n_iter
            Base.Threads.@spawn begin
                iter_count += 1
                idx_within = sample(Xoroshiro128Star(), 1:n_tot, n_within, replace = false)
                idx_among = n_range[Not(idx_within)]
                
                mu_within = mean(results[:, est][idx_within])
                mu_among = mean(results[:, est][idx_among])
                
                bootstrapped[i][j] = mu_within - mu_among
                
                pair_text =  population * "  $(iter_count)" * "/" * "$(iterations)"
                ProgressMeter.next!(p; showvalues = [(:Population, pair_text), (:Method, string(est))])
            end
        end
    end
    [(sum(differ[i] .<= bootstrapped[i]) + 1) / iterations for i in 1:length(estimators)]
    #(sum(differ .<= bootstrapped) + 1) / (iterations)
end


"""
    relatedness_posthoc(::PopData, results::DataFrame; iterations::Int)

Performs a posthoc analysis using the resulting DataFrame or NamedTuple
from relatedness(). This analysis uses permutations to test if a population has
significantly higher within-population relatedess than between-population relatedness.
The `results` object must have been generated from the provided `PopData`. Use `iterations = `
to specify the number of iterations for the permutation tests (default = `20000`). **Recommended**
that you use `MultipleTesting.jl` to correct resulting P-values.

**Example**
```
julia> cats = @nancycats ;

julia> rel_out = relatedness(cats, method = [Ritland, Moran], iterations = 100);

julia> relatedness_posthoc(cats, rel_out)
17x3 DataFrame
 Row │ population  Ritland_P  Moran_P
     │ String      Float64    Float64
─────┼────────────────────────────────
   1 │ 1              5.0e-5   5.0e-5
   2 │ 2              5.0e-5   5.0e-5
   3 │ 3              5.0e-5   5.0e-5
   4 │ 4              5.0e-5   5.0e-5
   5 │ 5              5.0e-5   5.0e-5
   6 │ 6              5.0e-5   5.0e-5
   7 │ 7              5.0e-5   5.0e-5
   8 │ 8              5.0e-5   5.0e-5
   9 │ 9              5.0e-5   5.0e-5
  10 │ 10             5.0e-5   5.0e-5
  11 │ 11             5.0e-5   5.0e-5
  12 │ 12             5.0e-5   5.0e-5
  13 │ 13             5.0e-5   5.0e-5
  14 │ 14             5.0e-5   5.0e-5
  15 │ 15             5.0e-5   5.0e-5
  16 │ 16             5.0e-5   5.0e-5
  17 │ 17             5.0e-5   5.0e-5
```
"""
function relatedness_posthoc(data::PopData, results::DataFrame; iterations::Int = 20000)
    all_pops = unique(data.meta.population)
    estimators = Symbol.(names(results)[names(results) .∉ Ref(["sample_1", "sample_2", "n_loci"])] .* "_P")
    sigs = map(pop -> sig_within(data, results, pop, iterations), all_pops)

    DataFrame(
        :population => all_pops,
        [estimators[i] => getindex.(sigs,i) for i in 1:length(estimators)]...
    )
end

function relatedness_posthoc(data::PopData, results::NamedTuple; iterations::Int = 20000)
    estimators = keys(results)
    coeffs = [results[i][:, estimators[i]] for i in 1:length(estimators)]
    df = select(results[1], "sample_1", "sample_2", "n_loci")
    [df[:, estimators[i]] = coeffs[i] for i in 1:length(estimators)]
    relatedness_posthoc(data, df, iterations = iterations)
end