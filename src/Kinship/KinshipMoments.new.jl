########  Moments-based Estimator Methods ###########

function _blouin(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8} 
    @inbounds ((geno1[1] ∈ geno2) & (geno2[1] ∈ geno1)) + ((geno1[2] ∈ geno2) & (geno2[2] ∈ geno1))
end
function Blouin(ind1::GenoArray, ind2::GenoArray)::Float64
    res = 0.0
    n = 0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            res += _blouin(i1, i2)
            n += 1
        end
    end
    return res / n / 2.0
end

function _lihorvitz(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8}
    @inbounds (geno1[1] == geno2[1]) + (geno1[1] == geno2[2]) + (geno1[2] == geno2[1]) + (geno1[2] == geno2[2]) 
end
function LiHorvitz(ind1::GenoArray, ind2::GenoArray)::Float64
    res = 0.0
    n = 0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            res += _lihorvitz(i1, i2)
            n += 1
        end
    end
    res / n / 4.0
end

function _lynch(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8}
    @inbounds ((geno1[1] ∈ geno2) + (geno1[2] ∈ geno2) + (geno2[1] ∈ geno1) + (geno2[2] ∈ geno1))
end
function Lynch(ind1::GenoArray, ind2::GenoArray)::Float64
    res = 0.0
    n = 0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            res += _lynch(i1, i2)
            n += 1
        end
    end
    res / n / 4.0
end


function _loiselle_num(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64})::Float64 where T<:Union{Int16, Int8}
    sum(pairs(frqdict)) do (allele, frq)
        @inbounds ((((geno1[1] == allele) + (geno1[2] == allele)) / 2.0) - frq) * ((((geno2[1] == allele) + (geno2[2] == allele)) / 2.0) - frq)     
    end
end
function _loiselle_denom(freqs)::Float64
    sum(freqs) do fq
        fq * (1.0-fq)
    end
end
function Loiselle(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numer = 0.0
    denom = 0.0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            @inbounds frqs = allelefrq[i]
            fq = values(frqs)
            numer += _loiselle_num(i1, i2, frqs)
            @inbounds denom += _loiselle_denom(fq)
        end
    end
    numer / denom + 2.0 / (2.0 * kwargs[:n_samples] - 1.0)
end

##BUG THIS DOESNT AGREE WITH THE ORIGINAL
## CHECK AGAINST COANCESTRY
function _lynchli(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8}
    a,b = geno1 ; c,d = geno2
    0.5 * (((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (a == b))) + ((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (c == d))))
end
function _lynchliS0(alleles)::Float64
    res1 = 0.0
    res2 = 0.0
    @inbounds for i in alleles
        res1 += i^2
        res2 += i^3
    end
    return 2.0 * res1 - res2
end
function LynchLi(ind1::T, ind2::T, alleles::U; kwargs...) where T <: GenoArray where U <: Tuple
    #TODO Change to unbiased formulation (eq 25)
    numerator = 0.0
    denom = 0.0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        @inbounds loc = values(alleles[i])
        S0 = _lynchliS0(loc)
        if (i1 === missing) | (i2 === missing)
            continue
        else
            # this is Sxy - S0
            numerator += _lynchli(i1, i2) - S0
            denom += 1.0 - S0
        end
    end
    return numerator / denom
end

#BUG not consistent with Coancestry
function _lynchritland(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    fq_a = frqdict[a]
    fq_b = frqdict[b]
    fq_c = frqdict[c]
    fq_d = frqdict[d]
    n1 = fq_a * ((b == c) + (b == d)) + fq_b * ((a == c) + (a == d)) - 4.0 * fq_a * fq_b
    n2 = fq_c * ((d == a) + (d == b)) + fq_d * ((c == a) + (c == b)) - 4.0 * fq_c * fq_d
    d1 = 2.0 * (1.0 + (a == b)) * (fq_a + fq_b) - 8.0 * fq_a * fq_b
    d2 = 2.0 * (1.0 + (c == d)) * (fq_c + fq_d) - 8.0 * fq_c * fq_d
    WL1 = ((1.0 + (a == b)) * (fq_a + fq_b) - 4.0 * fq_a * fq_b) / (2.0 * fq_a * fq_b)
    WL2 = ((1.0 + (c == d)) * (fq_c + fq_d) - 4.0 * fq_c * fq_d) / (2.0 * fq_c * fq_d)
    numer = ((n1 / d1) * WL1 + (n2 / d2) * WL2)
    denom = WL1 + WL2
    return (numer, denom)
end
function LynchRitland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numer = 0.0
    denom = 0.0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        @inbounds freqs = allelefrq[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            n, d = _lynchritland(i1, i2, freqs)
            numer += n
            denom += d
        end
    end    
    return numer / (denom / 2.0)
end

function _moran(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    num = 0.0 ; denom = 0.0
    @inbounds for (allele, fq) in frqdict
        g1 = ((((geno1[1] == allele) + (geno1[2] == allele)) / 2.0) - fq)
        g2 = ((((geno2[1] == allele) + (geno2[2] == allele)) / 2.0) - fq)
        num += g1 * g2
        denom += ((g1^2) + (g2^2))
    end
    return (num, denom)
end
function Moran(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    numer = 0.0
    denom = 0.0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        @inbounds freqs = allelefrq[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            n, d = _moran(i1, i2, freqs)
            numer += n
            denom += d
        end
    end    
    return numer/(denom / 2.0)
end

#BUG not consistent with Coancestry
function _quellergoodnight(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    ident = ((a == c) + (a == d) + (b == c) + (b == d))
    fq_a = frqdict[a]
    fq_b = frqdict[b]
    fq_c = frqdict[c]
    fq_d = frqdict[d]

    num1 = ident - 2.0 * (fq_a + fq_b)
    num2 = ident - 2.0 * (fq_c + fq_d)

    denom1 = (2.0 * (1.0 + (a==b) - fq_a - fq_b))
    denom2 = (2.0 * (1.0 + (c==d) - fq_c - fq_d))
    return (num1, num2, denom1, denom2)
end
function QuellerGoodnight(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numer1 = 0.0
    denom1 = 0.0
    numer2 = 0.0
    denom2 = 0.0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        @inbounds freqs = allelefrq[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            @inbounds n1, n2, d1, d2 = _quellergoodnight(i1, i2, freqs)
            numer1 += n1
            denom1 += d1
            numer2 += n2
            denom2 += d2
        end
    end    
    return (numer1/denom1 + numer2/denom2)/2.0
end

#BUG Math not consistent with Coancestry
function _ritland(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    A = length(frqdict) - 1.0
    R = 0.0
    for (allele, frq) in frqdict
        R += ((((a == allele) + (b == allele)) * ((c == allele) + (d == allele))) / (4.0 * frq))
    end
    R = ((2.0 / A) * (R - 1.0)) * A
    return (R, A)
end
function Ritland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numer = 0.0
    denom = 0.0
    @inbounds for i in eachindex(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        @inbounds freqs = allelefrq[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            n, d = _ritland(i1, i2, freqs)
            numer += n
            denom += d
        end
    end    
    return numer / denom
end