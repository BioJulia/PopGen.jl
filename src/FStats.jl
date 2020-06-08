"""
    f_stat_p(data::PopObj, nperm::Int64)
Performs a permutation test (permuting individuals among populations) to calculate
p-values for F-statistics
"""
function f_stat_p(data::PopData, nperm::Int64 = 1000)
    tmp = deepcopy(data)

    observed = f_stats(data, "global")
    sims = deepcopy(observed)
    for boot in 1:(nperm-1)
        shuffle!(tmp.samples.population)
        append!(sims, f_stats(tmp, "global"))
    end
    return map(i -> sum((observed .<= sims)[!,i])/nperm, 1:ncol(observed))
end
