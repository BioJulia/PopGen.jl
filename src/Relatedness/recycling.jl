

#### Rework for efficiency ####

function QuellerGoodnight(loc::Symbol, geno1::Genotype, geno2::Genotype, alleles::T) where T <: NamedTuple
    a,b = geno1
    c,d = geno2

    n1 = sum((a == c, a == d, b == c, b == d)) - 2.0 * (alleles[loc][a] + alleles[loc][b])
    n2 = sum((a == c, a == d, b == c, b == d)) - 2.0 * (alleles[loc][c] + alleles[loc][d])

    d1 = 2.0 * (1.0 + (a==b) - alleles[loc][a] - alleles[loc][b])
    d2 = 2.0 * (1.0 + (c==d) - alleles[loc][c] - alleles[loc][d])
    return (n1, d1, n2, d2)
end

function Ritland(loc::Symbol, geno1::Genotype, geno2::Genotype, alleles::T) where T <: NamedTuple
    a,b = geno1
    c,d = geno2

    A = ((alleles[loc] |> length) - 1)

    R = 0.0
    for allele in keys(alleles[loc])
        # Individual locus relatedness value (eq 7 in paper)
        R += ((((a == allele) + (b == allele)) * ((c == allele) + (d == allele))) / (4.0 * alleles[loc][allele]))
    end
    R = (2 / A) * (R - 1.0)
    # return numerator,denominator
    return ((R * A), A, 0.0, 0.0)
end

function LynchRitland(loc::Symbol, geno1::Genotype, geno2::Genotype, alleles::T) where T <: NamedTuple
    a,b = geno1
    c,d = geno2
    A = ((alleles[loc] |> length) - 1)

    n1 = alleles[loc][a] * ((b == c) + (b == d)) + alleles[loc][b] * ((a == c) + (a == d)) - 4.0 * alleles[loc][a] * alleles[loc][b]
    n2 = alleles[loc][c] * ((d == a) + (d == b)) + alleles[loc][d] * ((c == a) + (c == b)) - 4.0 * alleles[loc][c] * alleles[loc][d]

    d1 = 2.0 * (1.0 + (a == b)) * (alleles[loc][a] + alleles[loc][b]) - 8.0 * alleles[loc][a] * alleles[loc][b]
    d2 = 2.0 * (1.0 + (c == d)) * (alleles[loc][c] + alleles[loc][d]) - 8.0 * alleles[loc][c] * alleles[loc][d]

    RL = (n1 / d1) + (n2 / d2)
    # return numerator, denominator
    return (RL, A, 0.0, 0.0)
end

"""
    relatedness_moment(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness by Ritland (1996).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 10 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
Ritland original citation: https://www.cambridge.org/core/journals/genetics-research/article/estimators-for-pairwise-relatedness-and-individual-inbreeding-coefficients/9AE218BF6BF09CCCE18121AA63561CF7
"""
function relatedness_moment(data::PopData, ind1::String, ind2::String; alleles::T, method::Vector{Function}) where T <: NamedTuple
    #NEED TO CHECK TO CONFIRM EQUATIONS
    #Extract the pair of interest's genotypes
    gen1 = get_genotypes(data, ind1)
    gen2 = get_genotypes(data, ind2)

    d = Dict{Symbol,Vector{Float64}}()

    for (loc,geno1,geno2) in zip(skipmissings(Symbol.(loci(data)), gen1, gen2)...)
        for mthd in method
            nu, denm, nu2, denm2 = mthd(loc, geno1, geno2, alleles)
            get!(d, Symbol(mthd), zeros(4))[1] += nu
            get!(d, Symbol(mthd), zeros(4))[2] += denm
            get!(d, Symbol(mthd), zeros(4))[3] += nu2
            get!(d, Symbol(mthd), zeros(4))[4] += denm2
        end
    end
    if haskey(d, :QuellerGoodnight)
        qg = d[:QuellerGoodnight]
        d[:QuellerGoodnight][1] = (qg[1]/qg[2]) + (qg[3]/qg[4])
        d[:QuellerGoodnight][2] = 2.0
    end
    return NamedTuple{Tuple(keys(d))}(getindex.(values(d), 1) ./ getindex.(values(d), 2))
end

# matrix version?
function QuellerGoodnight2(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple

    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    numerator2 = 0.0
    denominator1 = 0.0
    denominator2 = 0.0

    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    # drop missing values
    loc, geno1, geno2 = collect.(skipmissings(Symbol.(loci(data)), geno1, geno2))

    # convert to 2-dimensional matrix for each individual
    gen1 = hcat(getindex.(geno1, 1), getindex.(geno1, 2))
    gen2 = hcat(getindex.(geno2, 1), getindex.(geno2, 2))
    # ((a == c) + (a == d) + (b == c) + (b == d))
    id_mtx = reduce(*, (gen1 .== gen2) + (gen1 .== reverse(gen2, dims = 2)), dims = 2)
    a_b_frqs = [alleles[loc[i]][gen1[i,j]] for i in 1:size(gen1)[1], j in 1:2]
    c_d_frqs = [alleles[loc[i]][gen2[i,j]] for i in 1:size(gen2)[1], j in 1:2]
    
    numerator1 = reduce(+, id_mtx - (2 * reduce(+, a_b_frqs, dims = 2))) 
    numerator2 = reduce(+, id_mtx - (2 * reduce(+, c_d_frqs, dims = 2)))
    # reduce(+,)
    denominator1 = 2.0 * (1.0 + map(x ->x[1] == x[2],eachrow(gen1)) - foldl(-, eachcol(a_b_frqs)))
    denominator2 = 2.0 * (1.0 + map(x ->x[1] == x[2],eachrow(gen2)) - foldl(-, eachcol(a_b_frqs)))
    #return [alleles[loc[i]][gen1[i,j]] for i in 1:size(gen1)[1], j in 1:2]
    return numerator2
    a,b = gen1
    c,d = gen2

        ident = ((a == c) + (a == d) + (b == c) + (b == d))
        numerator1 += ident - 2.0 * (alleles[loc][a] + alleles[loc][b])
        numerator2 += ident - 2.0 * (alleles[loc][c] + alleles[loc][d])

        denominator1 += (2.0 * (1.0 + (a==b) - alleles[loc][a] - alleles[loc][b]))
        denominator2 += (2.0 * (1.0 + (c==d) - alleles[loc][c] - alleles[loc][d]))

    return (numerator1/denominator1 + numerator2/denominator2)/2.0
end

mapreduce(x->x^2, +, [1:3;])
mapreduce(isodd, *, a, dims=1)
z = (a .== b) .+ (a .== reverse(b, dims = 2))
mapreduce(x -> x, +, (a .== b) .+ (a .== reverse(b, dims = 2)), dims = 2)