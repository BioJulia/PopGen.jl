function kinship_noboot(data::PopData; method::Function, kwargs...)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    result = NamedArray(zeros(Float64, n, n))
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    if Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch]
        @inbounds for i in 1:n-1
            @inbounds v1 = view(locmtx,i,:)
            @inbounds for j in i+1:n
                @inbounds v2 = view(locmtx,j,:)
                est = method(v1, v2)
                @inbounds result[j,i] = result[i,j] = est
            end
        end
    elseif Symbol(method) ∈ [:Loiselle, :Loiselle2, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland]
        allelefrequencies = @inbounds (allelefreq(i) for i in eachcol(locmtx)) |> Tuple
        @inbounds for i in 1:n-1
            @inbounds v1 = view(locmtx,i,:)
            @inbounds for j in i+1:n
                @inbounds v2 = view(locmtx,j,:)
                est = method(v1, v2, allelefrequencies, n_samples = n)
                @inbounds result[j,i] = result[i,j] = est
            end
        end
    else
        throw(ArgumentError("Method $method is not a valid method. See ?kinship for a list of options."))    
    end
    return result
end 

using OnlineStats
using Term.progress

function kinship_boot(data::PopData; method::Function, iterations::Int, kwargs...)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    nloc = size(locmtx, 2)
    idxrange = 1:nloc
    result = NamedArray(zeros(Float64, n, n))
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    bresult = Matrix{Float64}(undef, n,n)
    b_sdev = similar(bresult)
    b_ci = Matrix{Vector{Float64}}(undef, n,n)
    #if Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch]
    pbar = ProgressBar(refresh_rate=30)
    samplejob = addjob!(pbar, N = n)
    job = addjob!(pbar; N= Int64((n * (n-1))/2))
    start!(pbar)
    @inbounds @sync  for i in 1:n-1
        Base.Threads.@spawn begin
            @inbounds v1 = view(locmtx,i,:)
            @inbounds for j in i+1:n
                boot_idx = Vector{Int64}(undef, nloc)
                @inbounds v2 = view(locmtx,j,:)
                est = method(v1, v2)
                @inbounds result[j,i] = result[i,j] = est
                bootstats = Series(Variance(), Quantile([0.025, 0.975]))
                @inbounds for bt_idx in 1:iterations
                    @inbounds boot_idx .= rand(idxrange, nloc)
                    v1b = @inbounds view(locmtx, i,boot_idx)
                    v2b = @inbounds view(locmtx, j, boot_idx)
                    b_est = method(v1b, v2b)
                    isnan(b_est) ? continue : fit!(bootstats, b_est)
                end
                @inbounds bresult[i,j] = bootstats.stats[1].μ
                # median
                @inbounds b_sdev[i,j] = sqrt(bootstats.stats[1].σ2)
                # Stdev
                @inbounds b_ci[i,j] = value(bootstats.stats[2])
                progress.update!(job)
            end
            progress.update!(samplejob)
        end
    end
    stop!(pbar)
    ci = uppertri2vec(b_ci)
    cilow = getindex.(ci, 1)
    cihi = getindex.(ci,2)
    out = kinshiptotable(result, Symbol(method))
    insertcols!(out, :bootmean => uppertri2vec(bresult), :std => uppertri2vec(b_sdev), :CIlower => cilow, :CLupper => cihi)
    return out    
end


########  Estimator Methods ###########



function _blouin(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8} 
    @inbounds ((geno1[1] ∈ geno2) & (geno2[1] ∈ geno1)) + ((geno1[2] ∈ geno2) & (geno2[2] ∈ geno1))
end
function Blouin(ind1::GenoArray, ind2::GenoArray)::Float64
    l = length(ind1)
    res = 0.0
    n = 0
    @inbounds for i in 1:l
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

# TODO check math, should diag = 1?
function _lihorvitz(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8}
    @inbounds (geno1[1] == geno2[1]) + (geno1[1] == geno2[2]) + (geno1[2] == geno2[1]) + (geno1[2] == geno2[2]) 
end
function LiHorvitz(ind1::GenoArray, ind2::GenoArray)::Float64
    l = length(ind1)
    res = 0.0
    n = 0
    @inbounds for i in 1:l
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
    l = length(ind1)
    res = 0.0
    n = 0
    @inbounds for i in 1:l
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
        ((((geno1[1] == allele) + (geno1[2] == allele)) / 2.0) - frq) * ((((geno2[1] == allele) + (geno2[2] == allele)) / 2.0) - frq)     
    end
end
function _loiselle_denom(freqs)::Float64
    sum(freqs) do fq
        fq * (1.0-fq)
    end
end
function Loiselle(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    l = length(ind1)
    numer = 0.0
    denom = 0.0
    #n = 0
    @inbounds for i in 1:l
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        if (i1 === missing) | (i2 === missing)
            continue
        else
            @inbounds frqs = allelefrq[i]
            fq = values(frqs)
            numer += _loiselle_num(i1, i2, frqs)
            @inbounds denom += _loiselle_denom(fq)
            #n += 1
        end
    end
    numer / denom + 2.0 / (2.0 * kwargs[:n_samples] - 1.0)
end


function _lynchli(geno1::NTuple{2,T}, geno2::NTuple{2,T})::Float64 where T<:Union{Int16, Int8}
    a,b = geno1 ; c,d = geno2
    0.5 * (((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (a == b))) + ((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (c == d))))
end
function _lynchliS0(alleles)::Float64
    res1 = 0.0
    res2 = 0.0
    for i in alleles
        res1 += i^2
        res2 += i^3
    end
    return 2.0 * res1 - res2
end
function LynchLi(ind1::T, ind2::T, alleles::U; kwargs...) where T <: GenoArray where U <: Tuple
    #TODO Change to unbiased formulation (eq 25)
    Sxy = 0.0
    S0 = 0.0
    @inbounds for i in 1:length(ind1)
        @inbounds i1 = ind1[i]
        @inbounds i2 = ind2[i]
        @inbounds loc = values(alleles[i])
        S0tmp = _lynchliS0(loc)
        S0 += 1.0 - S0tmp
        if (i1 === missing) | (i2 === missing)
            continue
        else
            Sxy += _lynchli(i1, i2) - S0tmp
        end
    end
    return Sxy / S0
end


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
    @inbounds for i in 1:length(ind1)
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

# FIX _MORAN MATH
function _moran(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    num = 0.0 ; denom = 0.0
    @inbounds for (allele, fq) in pairs(frqdict)
        num += (((geno1[1] == allele) + (geno1[2] == allele)) / 2.0) * (((geno2[1] == allele) + (geno2[2] == allele)) / 2.0)
        #num += @inbounds ((sum(geno1 .== allele) / 2.0) - fq) * ((sum(geno2 .== allele) / 2.0) - fq)
        denom += ((((geno1[1] == allele) + (geno1[2] == allele)) / 2.0) - fq)^2 + ((((geno2[1] == allele) + (geno2[2] == allele)) / 2.0) - fq)^2
        #denom += @inbounds (((sum(geno1 .== allele) / 2.0) - fq)^2 + ((sum(geno2 .== allele) / 2.0) - fq)^2)
    end
    return (num, denom)
end
function Moran(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    numer = 0.0
    denom = 0.0
    @inbounds for i in 1:length(ind1)
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
    @inbounds for i in 1:length(ind1)
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


function _Ritland(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    a,b = geno1
    c,d = geno2
    A = length(frqdict) - 1.0
    R = 0.0
    for (allele, frq) in pairs(frqdict)
        R += ((((a == allele) + (b == allele)) * ((c == allele) + (d == allele))) / (4.0 * frq))
    end
    R = ((2.0 / A) * (R - 1.0)) * A
    return (R, A)
end
function Ritland(ind1::GenoArray, ind2::GenoArray, allelefrq::U; kwargs...) where U <: Tuple
    numer = 0.0
    denom = 0.0
    @inbounds for i in 1:length(ind1)
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
function kinshiptotable(kinshipresults::T, mthd::Symbol) where T<:NamedMatrix
    n = size(kinshipresults,1)
    n == size(kinshipresults,2) || throw(DimensionMismatch("The input matrix must be symmetrical, but has size $(size(kinshipresults))"))
    ids = names(kinshipresults)[1]
    vals = uppertri2vec(kinshipresults)
    idpairs = pairwisepairs(ids)
    DataFrame(:sample1 => first.(idpairs), :sample2 => getindex.(idpairs, 2), mthd => vals)
end

function uppertri2vec(x, diag::Bool = false)
    n = size(x,1)
    if !diag
        [x[i, j] for i in 1:n-1 for j in i+1:n]
    else
        [x[i, j] for i in 1:n-1 for j in i:n]
    end
end

function lowertri2vec(x, diag::Bool = false)
    n = size(x,1)
    if !diag
        [x[j,i] for i in 1:n-1 for j in i+1:n]
    else
        [x[j,i] for i in 1:n-1 for j in i:n]
    end
end