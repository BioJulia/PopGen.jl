using BenchmarkTools
data = nancycats();
data2 = gulfsharks();

struct Locus
    alleles::Vector{Union{Int16,In8}}
    pops::Vector{String}
    freqs::Vector{Union{Missing,Dict{String,Float64}}}
end

"""
    allele_freq(locus::SubArray{<:Union{Missing,NTuple{N,<:Integer}}}) where N
Calculate allele counts for a single locus of a `PopObj` split by population
using `group()`. Returns a `Dict` of allele's and their frequencies.
"""
function allele_freq_sub(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
    d = Dict{String,Float64}()
    flat_alleles = Base.Iterators.flatten(locus |> skipmissing) |> collect
    uniq_alleles = unique(flat_alleles)
    allele_counts = [count(i -> i == j, flat_alleles) for j in uniq_alleles]
    total = sum(allele_counts)
    freq = allele_counts ./ total
#=
    for (i,j) in enumerate(uniq_alleles)
        d[string(j)] = freq[i]
    end
    return d
=#
end


## heterozygosity of a locus in a single pop
function het_ob(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}},1}) where N
    a = geno_freq(locus)  # get genotype freqs at locus
    tmp = 0
    genos = keys(a) |> collect
    @inbounds for geno in genos
        length(unique(geno)) == 1 && continue
        tmp += a[geno]
    end
    return tmp
end

y_2 = deepcopy(ncats.loci)
insertcols!(y_2, 1, :population => ncats.samples.population)

### this might potentially work!
a = aggregate(y, :population, het_ob)

# count ind
c = aggregate(y, :population, i -> length(skipmissing(i) |> collect))


b = map(mean, eachcol(a[!,:2:end]))


y_2 = deepcopy(ncats.loci)
insertcols!(y_2, 1, :population => ncats.samples.population)
a2 = aggregate(y, :population, het_ob)
c2 = aggregate(y_2, :population, i -> length(skipmissing(i) |> collect))
