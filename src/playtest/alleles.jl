using BenchmarkTools
data = nancycats();
data2 = gulfsharks();

struct Locus
    alleles::Vector{Union{Int16,In8}}
    pops::Vector{String}
    freqs::Vector{Union{Missing,Dict{String,Float64}}}
end

"""
    allele_freq_sub(locus::SubArray{<:Union{Missing,NTuple{N,<:Integer}}}) where N
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
    for (i,j) in enumerate(uniq_alleles)
        d[string(j)] = freq[i]
    end
    return d
end


## heterozygosity of a locus in a single pop
function het_obs(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}},1}) where N
    all(ismissing.(locus)) == true && return missing
    a = geno_freq(locus)  # get genotype freqs at locus
    tmp = 0
    genos = keys(a) |> collect
    @inbounds for geno in genos
        length(unique(geno)) == 1 && continue
        tmp += a[geno]
    end
    return tmp
end


y = deepcopy(data.loci)
insertcols!(y, 1, :population => data.samples.population)
a = aggregate(y, :population, het_obs);

# count ind
c = aggregate(y, :population, i -> length(skipmissing(i) |> collect))


b = map(mean, eachcol(a[!,:2:end]))


y_2 = deepcopy(data2.loci);
insertcols!(y_2, 1, :population => data2.samples.population);
c2 = aggregate(y_2, :population, i -> length(skipmissing(i) |> collect))

loc = data.loci.fca45
all(ismissing.(locus)) == true && return missing
loc_corr = loc |> skipmissing |> collect
len = length(loc_corr)
decomposed_genos = length.(unique.(loc_corr))
sum(decomposed_genos .> 1) / len
