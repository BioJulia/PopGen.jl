function fstat(data::PopData)
    # observed/expected het per locus per pop
    het_df = DataFrames.combine(
        groupby(data.loci, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => (i -> hetero_e(i)) => :het_exp,
        :genotype => allele_freq => :alleles
    )
    transform!(het_df,
        #[:het_obs, :het_exp, :n] => ((o,e,n) -> (e .- o ./ 2 ./ n)) => :Hs,
        #[:het_obs, :het_exp, :n] => ((o,e,n) -> (e .- o ./ 2 ./ n) .* (n ./ (n .- 1))) => :Hs_cor,
        [:het_obs, :het_exp, :n] => ((o,e,n) -> gene_diversity_nei87.(e,o,n, corr = false)) => :Hs,
        [:het_obs, :het_exp, :n] => ((o,e,n) -> gene_diversity_nei87.(e,o,n)) => :Hs_nei,
    )
    transform!(
        het_df,
        [:het_obs, :Hs] => ((Ho, Hs) -> 1 .- (Ho./Hs)) => :FIS
    )

    # number of populations each locus appears in
    n_df = DataFrames.combine(
        groupby(het_df, :locus),
        :n => (i -> sum(.!iszero.(i))) => :count,
        :n => (i -> sum(reciprocal.(i))) => :N,
        :alleles => (i -> sum(values(avg_allele_freq(i)).^2)) => :avg_freq,
        :alleles => (i -> sum(values(avg_allele_freq(i)).^2)) => :avg_freq
    )
    transform!(
        n_df,
        [:count, :N] =>((np, N) -> np ./ N) => :mn
    )
    transform!(
        n_df,
        [:]
    )

    return n_df
end



z = DataFrames.combine(
    groupby(x.loci, [:locus, :population]),
    :genotype => allele_freq => :alleles
)

DataFrames.combine(
       groupby(z, :locus),
       :mp2 => (i -> sum(i.^2)) => :mp
       )


function avg_allele_freq(allele_dicts::AbstractVector{T}) where T<:Dict{Int16,Float32}
   sum_dict = Dict{Int16, Tuple{Float32, Int}}()
   all_alleles = keys.(allele_dicts) |> Base.Iterators.flatten |> collect |> unique
   @inbounds for allele in all_alleles
       for allele_dict in allele_dicts
           sum_dict[allele] = get!(sum_dict, allele, (0., 0)) .+ (get!(allele_dict, allele, 0.), 1)
       end
   end
   avg_dict = Dict{Int16, Float32}()
   @inbounds for (key, value) in sum_dict
       freq_sum, n = value
       if !iszero(freq_sum)
           @inbounds avg_dict[key] = freq_sum / n
       end
   end
   return avg_dict
end
