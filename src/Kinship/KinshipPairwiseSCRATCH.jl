function kinship_noboot(data::PopData; method::Function, kwargs...)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    result = NamedArray(Array{ZZ(Float64, n, n))
    
    #result = NamedArray(zeros(Float64, n, n))
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    if Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch]
        for i in 1:n-1
            v1 = view(locmtx,i,:)
            for j in i+1:n
                result[i,j] = method(v1, view(locmtx,j,:))
            end
        end
    elseif Symbol(method) ∈ [:Loiselle, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland]
        allelefrequencies = allelefreqtuple(data)
        @inbounds for i in 1:n-1
            @inbounds v1 = view(locmtx,i,:)
            @inbounds for j in i+1:n
                @inbounds v2 = view(locmtx,j,:)
                @inbounds result[i,j] = method(v1, v2, allelefrequencies, n_samples = n)
            end
        end
    else
        throw(ArgumentError("Method $method is not a valid method. See ?kinship for a list of options."))    
    end
    return result
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


function _Loiselle(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64})::Float64 where T<:Union{Int16, Int8}
    @inbounds sum(skipmissing([((sum(geno1 .== allele) / 2.0) - frqdict[allele]) * ((sum(geno2 .== allele) / 2.0) - frqdict[allele]) for allele in keys(frqdict)]))
end
_Loiselle(geno1::Missing, geno2::Genotype, frqdict::Dict) = missing
_Loiselle(geno1::Genotype, geno2::Missing, frqdict::Dict) = missing
_Loiselle(geno1::Missing, geno2::Missing, frqdict::Dict) = missing
function Loiselle(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    kw = Dict(kwargs)
    num = @inbounds _Loiselle.(ind1, ind2, allelefrq)
    denom = @inbounds [sum(values(frqs) .* (1 .- values(frqs))) for frqs in allelefrq]
    return sum(skipmissing(num)) / sum(skipmissing(denom)) + 2.0 / (2.0 * kw[:n_samples] - 1.0)
end


function _LynchLi(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8}
    a,b = geno1 ; c,d = geno2
    0.5 * (((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (a == b))) + ((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (c == d))))
end
_LynchLi(geno1::Missing, geno2::Genotype) = missing
_LynchLi(geno1::Genotype, geno2::Missing) = missing
_LynchLi(geno1::Missing, geno2::Missing) = missing
function LynchLi(ind1::T, ind2::T, alleles::U; kwargs...) where T <: GenoArray where U <: Tuple
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

function LynchRitland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numdenom = skipmissing(_LynchRitland.(ind1, ind2, allelefrq))
    numer = sum(getindex.(numdenom, 1))
    denom = sum(getindex.(numdenom, 2))
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

function Moran(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    numdenom = skipmissing(_Moran.(ind1, ind2, allelefrq))
    numer = sum(getindex.(numdenom, 1))
    denom = sum(getindex.(numdenom, 2))
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

function QuellerGoodnight(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numdenom = skipmissing(_QuellerGoodnight.(ind1, ind2, allelefrq))
    numer1 = sum(getindex.(numdenom, 1))
    numer2 = sum(getindex.(numdenom, 2))
    denom1 = sum(getindex.(numdenom, 3))
    denom2 = sum(getindex.(numdenom, 4))
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

function Ritland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numdenom = skipmissing(_Ritland.(ind1, ind2, allelefrq))
    numer = sum(getindex.(numdenom, 1))
    denom = sum(getindex.(numdenom, 2))
    return numer / denom
end


"""
    kinshiptotable(kinshipresults::T) where T<:NamedMatrix
Converts the `NamedMatrix` result from the `kinship()` function into a `DataFrame`.

**Example**
```julia`
julia> cats = @nancycats ; kin = kinship(cats, method = Moran) ;

julia> kinshiptotable(a)
22366×3 DataFrame
   Row │ sample1  sample2  kinship      
       │ String   String   Float64      
───────┼────────────────────────────────
     1 │ cc_001   cc_002    0.00688008
     2 │ cc_001   cc_003   -0.0286812
     3 │ cc_001   cc_005   -0.000749142
     4 │ cc_001   cc_007    0.0516361
     5 │ cc_001   cc_008    0.0261128
     6 │ cc_001   cc_009   -0.00187027
     7 │ cc_001   cc_010    0.0182852
   ⋮   │    ⋮        ⋮          ⋮
 22361 │ seg_028  seg_029  -0.0472928
 22362 │ seg_028  seg_030  -0.0172853
 22363 │ seg_028  seg_031  -0.00240921
 22364 │ seg_029  seg_030  -0.0278483
 22365 │ seg_029  seg_031   0.0297876
 22366 │ seg_030  seg_031  -0.0371295
                      22353 rows omitted
```
"""
function kinshiptotable(kinshipresults::T) where T<:NamedMatrix
    ids = names(kinshipresults)[1]
    n = size(kinshipresults,1)
    vals = [kinshipresults[i, j] for i in 1:n-1 for j in i+1:n]
    idpairs = collect(pairwisepairs(ids))
    DataFrame(:sample1 => first.(idpairs), :sample2 => getindex.(idpairs, 2), :kinship => vals)
end


#TODO move to popgencore
function allelefreqtuple(data::PopData)
    Tuple(DataFrames.combine(
        groupby(data.genodata, :locus),
        :genotype => allelefreq => :frq
        )[:, :frq]
    )
end