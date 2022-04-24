function kinship_noboot(data::PopData; method::Function, kwargs...)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    vecs = [i for i in eachrow(locmtx)]
    if Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch]
        out = NamedArray(method.(vecs, permutedims(vecs)))
    elseif Symbol(method) ∈ [:Loiselle, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland]
        allelefrequencies = allelefreq(data)
        out = NamedArray(method.(vecs, permutedims(vecs), Ref(allelefrequencies), n_samples = length(ids)))
    else
        throw(ArgumentError("Method $method is not a valid method. See ?kinship for a list of options."))
    
    end
    setnames!(out, String.(ids),1)
    setnames!(out, String.(ids),2)
    return out
end 

_Blouin(geno1::NTuple{2,T}, geno2::NTuple{2,T}) where T<:Union{Int16, Int8} = @inbounds ((geno1[1] ∈ geno2) & (geno2[1] ∈ geno1)) + ((geno1[2] ∈ geno2) & (geno2[2] ∈ geno1))
_Blouin(geno1::Missing, geno2::Genotype) = missing
_Blouin(geno1::Genotype, geno2::Missing) = missing
_Blouin(geno1::Missing, geno2::Missing) = missing
Blouin(ind1::GenoArray, ind2::GenoArray)::Float64 = @inbounds mean(skipmissing(_Blouin.(ind1, ind2) ./ 2))

# TODO check math, diagonal should = 1
_LiHorvitz(geno1::NTuple{2,T}, geno2::NTuple{2,T}) where T<:Union{Int16, Int8} = sum(geno1 .∈ collect(geno2)')
_LiHorvitz(geno1::Missing, geno2::Genotype) = missing
_LiHorvitz(geno1::Genotype, geno2::Missing) = missing
_LiHorvitz(geno1::Missing, geno2::Missing) = missing
LiHorvitz(ind1::GenoArray, ind2::GenoArray)::Float64 = @inbounds mean(skipmissing(_LiHorvitz.(ind1, ind2) ./4))

_Lynch(geno1::NTuple{2,T}, geno2::NTuple{2,T}) where T<:Union{Int16, Int8} = @inbounds ((geno1[1] ∈ geno2) + (geno1[2] ∈ geno2) + (geno2[1] ∈ geno1) + (geno2[2] ∈ geno1))
_Lynch(geno1::Missing, geno2::Genotype) = missing
_Lynch(geno1::Genotype, geno2::Missing) = missing
_Lynch(geno1::Missing, geno2::Missing) = missing
Lynch(ind1::GenoArray, ind2::GenoArray)::Float64 = @inbounds mean(skipmissing(_Lynch.(ind1, ind2) ./ 4))


function _Loiselle(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{IntT, Float64}) where T<:Union{Int16, Int8}
    @inbounds sum(skipmissing([((sum(geno1 .== allele) / 2.0) - frqdict[allele]) * ((sum(geno2 .== allele) / 2.0) - frqdict[allele]) for allele in keys(frqdict)]))
end
_Loiselle(geno1::Missing, geno2::Genotype, frqdict::Dict) = missing
_Loiselle(geno1::Genotype, geno2::Missing, frqdict::Dict) = missing
_Loiselle(geno1::Missing, geno2::Missing, frqdict::Dict) = missing
function Loiselle(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: NamedTuple
    kw = Dict(kwargs)
    frqdict = collect(allelefrq)
    num = _Loiselle.(ind1, ind2, frqdict)
    denom = [sum(values(frqs) .* (1 .- values(frqs))) for frqs in frqdict]
    return sum(skipmissing(num)) / sum(skipmissing(denom)) + 2.0 / (2 * kw[:n_samples] - 1)
end


function _LynchLi(geno1::NTuple{2,T}, geno2::NTuple{2,T}) where T<:Union{Int16, Int8}
    a,b = geno1 ; c,d = geno2
    0.5 * (((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (a == b))) + ((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (c == d))))
end
_LynchLi(geno1::Missing, geno2::Genotype) = missing
_LynchLi(geno1::Genotype, geno2::Missing) = missing
_LynchLi(geno1::Missing, geno2::Missing) = missing
function LynchLi(ind1::T, ind2::T, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    Sxy = _LynchLi.(ind1, ind2)
    #TODO Change to unbiased formulation (eq 25)
    S0 = @inbounds [2.0 * sum(values(loc) .^ 2) - sum(values(loc) .^ 3) for loc in alleles]
    return @inbounds sum(skipmissing(Sxy .- S0)) / sum(1.0 .- S0)
end


function _LynchRitland(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    fq_a, fq_b, fq_c, fq_d = map(i -> frqdict[i], (a,b,c,d))
    n1 = fq_a * ((b == c) + (b == d)) + fq_b * ((a == c) + (a == d)) - 4.0 * fq_a * fq_b
    n2 = fq_c * ((d == a) + (d == b)) + fq_d * ((c == a) + (c == b)) - 4.0 * fq_c * fq_d
    d1 = 2.0 * (1.0 + (a == b)) * (fq_a + fq_b) - 8.0 * fq_a * fq_b
    d2 = 2.0 * (1.0 + (c == d)) * (fq_c + fq_d) - 8.0 * fq_c * fq_d
    WL1 = ((1 + (a == b)) * (fq_a + fq_b) - 4 * fq_a * fq_b) / (2 * fq_a * fq_b)
    WL2 = ((1 + (c == d)) * (fq_c + fq_d) - 4 * fq_c * fq_d) / (2 * fq_c * fq_d)
    numer = ((n1 / d1) * WL1 + (n2 / d2) * WL2)
    denom = WL1 + WL2
    return (numer, denom)
end
_LynchRitland(geno1::Missing, geno2::Genotype, frqdict::Dict) = missing
_LynchRitland(geno1::Genotype, geno2::Missing, frqdict::Dict) = missing
_LynchRitland(geno1::Missing, geno2::Missing , frqdict::Dict) = missing

function LynchRitland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: NamedTuple
    frqdict = collect(allelefrq)
    numdenom = _LynchRitland.(ind1, ind2, frqdict)
    numer = sum(getindex.(skipmissing(numdenom), 1))
    denom = sum(getindex.(skipmissing(numdenom), 2))
    return numer / (denom ./ 2.0)
end

function _Moran(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    num = 0.0 ; denom = 0.0
    for allele in keys(frqdict)
        fq = frqdict[allele]
        num += ((sum(geno1 .== allele) / 2.0) - fq) * ((sum(geno2 .== allele) / 2.0) - fq)
        denom += (((sum(geno1 .== allele) / 2.0) - fq)^2 + ((sum(geno2 .== allele) / 2.0) - fq)^2) / 2.0
    end
    return (num, denom)
end
_Moran(geno1::Missing, geno2::Genotype, frqdict::Dict) = missing
_Moran(geno1::Genotype, geno2::Missing, frqdict::Dict) = missing
_Moran(geno1::Missing, geno2::Missing , frqdict::Dict) = missing

function Moran(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    frqdict = collect(allelefrq)
    numdenom = _Moran.(ind1, ind2, frqdict)
    numer = sum(getindex.(skipmissing(numdenom), 1))
    denom = sum(getindex.(skipmissing(numdenom), 2))
    return numer/denom
end


function _QuellerGoodnight(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    ident = ((a == c) + (a == d) + (b == c) + (b == d))
    fq_a, fq_b, fq_c, fq_d = map(i -> frqdict[i], (a,b,c,d))

    num1 = ident - 2.0 * (fq_a + fq_b)
    num2 = ident - 2.0 * (fq_c + fq_d)

    denom1 = (2.0 * (1.0 + (a==b) - fq_a - fq_b))
    denom2 = (2.0 * (1.0 + (c==d) - fq_c - fq_d))
    return (num1, num2, denom1, denom2)
end
_QuellerGoodnight(geno1::Missing, geno2::Genotype, frqdict::Dict) = missing
_QuellerGoodnight(geno1::Genotype, geno2::Missing, frqdict::Dict) = missing
_QuellerGoodnight(geno1::Missing, geno2::Missing , frqdict::Dict) = missing

function QuellerGoodnight(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: NamedTuple
    frqdict = collect(allelefrq)
    numdenom = _QuellerGoodnight.(ind1, ind2, frqdict)
    numer1 = sum(getindex.(skipmissing(numdenom), 1))
    numer2 = sum(getindex.(skipmissing(numdenom), 2))
    denom1 = sum(getindex.(skipmissing(numdenom), 3))
    denom2 = sum(getindex.(skipmissing(numdenom), 4))
    return (numer1/denom1 + numer2/denom2)/2.0
end


function _Ritland(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    alle = (a,b,c,d)
    A = length(frqdict) - 1
    R = 0.0
    for i in unique(alle)
        # Individual locus relatedness value (eq 7 in paper)
        R += sum(i .== alle) / (4.0 * frqdict[i])
        #R += ((((a == i) + (b == i)) * ((c == i) + (d == i))) / (4.0 * alleles[loc][i]))
    end
    R = (2.0 / A) * (R - 1.0)
    # numerator for weighted combination of loci
    num = (R * A)
    # denominator for weighted combination of loci
    denom = A
    return (num, denom)
end
_Ritland(geno1::Missing, geno2::Genotype, frqdict::Dict) = missing
_Ritland(geno1::Genotype, geno2::Missing, frqdict::Dict) = missing
_Ritland(geno1::Missing, geno2::Missing , frqdict::Dict) = missing

function Ritland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: NamedTuple
    frqdict = collect(allelefrq)
    numdenom = _Ritland.(ind1, ind2, frqdict)
    numer = sum(getindex.(skipmissing(numdenom), 1))
    denom = sum(getindex.(skipmissing(numdenom), 2))
    return numer / denom
end