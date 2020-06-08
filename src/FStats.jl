"""
    f_stat_p(data::PopObj, nperm::Int64)
Performs a permutation test (permuting individuals among populations) to calculate
p-values for F-statistics
"""
function f_stat_p(data::PopData, nperm::Int = 1000)
    tmp = deepcopy(data)
    observed = summary(data, by = "global")
    sims = deepcopy(observed)
    for permutation in 1:(nperm-1)
        permute_samples!(tmp)
        summary(tmp, by = "global")
        #append!(sims, summary(tmp, by = "global"))
    end
    return
    return map(i -> sum((observed .<= sims)[!,i])/nperm, 1:ncol(observed))
end
