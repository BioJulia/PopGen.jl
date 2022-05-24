function kinship_new(data::PopData; method::Function, iterations::Int = 0, interval::Vector{Float64} = [0.025, 0.975])
    # sanity checks
    data.metadata.ploidy != 2 && error("kinship analyses currently only support diploid samples")
    length(interval) != 2 && throw(ArgumentError("Keyword argument \`interval\` must be a vector with 2 elements."))
    
    if (iterations == 0) && (Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch])
        _kinship_noboot_nofreq(data, method)
    elseif (iterations == 0) && (Symbol(method) ∈ [:Loiselle, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland])
        _kinship_noboot_freq(data, method)
    elseif (iterations > 0) && (Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch])
        _kinship_boot_nofreq(data, method, iterations, interval)
    elseif (iterations > 0) && (Symbol(method) ∈ [:Loiselle, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland])
        _kinship_boot_freq(data, method, iterations, interval)
    else
        throw(ArgumentError("Invalid method provided: $method. See the docstring \`?kinship\` for usage information."))
    end
end

function kinship_new(data::PopData, samplenames::AbstractVector{T}; method::Function, iterations::Int = 0, interval::Vector{Float64} = [0.025, 0.975]) where T<:AbstractString
    # sanity checks
    data.metadata.ploidy != 2 && error("kinship analyses currently only support diploid samples")
    length(interval) != 2 && throw(ArgumentError("Keyword argument \`interval\` must be a vector with 2 elements."))
    newdata = data[data.genodata.name .∈ Ref(samplenames)]
    if (iterations == 0) && (Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch])
        _kinship_noboot_nofreq(newdata, method)
    elseif (iterations == 0) && (Symbol(method) ∈ [:Loiselle, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland])
        _kinship_noboot_freq(newdata, method)
    elseif (iterations > 0) && (Symbol(method) ∈ [:Blouin, :LiHorvitz, :Lynch])
        _kinship_boot_nofreq(newdata, method, iterations, interval)
    elseif (iterations > 0) && (Symbol(method) ∈ [:Loiselle, :LynchLi, :LynchRitland, :Moran, :QuellerGoodnight, :Ritland])
        _kinship_boot_freq(newdata, method, iterations, interval)
    else
        throw(ArgumentError("Invalid method provided: $method. See the docstring \`?kinship\` for usage information."))
    end
end

### Internal implementations ###

#### Non-Bootstrapped ####
# Returns a NamedMatrix
function _kinship_noboot_nofreq(data::PopData, method::Function)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    result = NamedArray{Float64}(n, n)
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    @inbounds for i in 1:n-1
        @inbounds v1 = view(locmtx,i,:)
        @inbounds for j in i+1:n
            @inbounds v2 = view(locmtx,j,:)
            est = method(v1, v2)
            @inbounds result[j,i] = result[i,j] = est
        end
    end
    return result
end 

function _kinship_noboot_freq(data::PopData, method::Function)
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    result = NamedArray{Float64}(n, n)
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    allelefrequencies = @inbounds Tuple(allelefreq(i) for i in eachcol(locmtx))
    @inbounds for i in 1:n-1
        @inbounds v1 = view(locmtx,i,:)
        @inbounds for j in i+1:n
            @inbounds v2 = view(locmtx,j,:)
            @inbounds result[j,i] = result[i,j] = method(v1, v2, allelefrequencies, n_samples = n)
        end
    end
    return result
end 

#### Bootstrapped ####
# Returns a DataFrame
# Includes a progress bar from Term.jl
# Uses OnlineStats.jl to calculate mean/variance/CI without additional allocations
function _kinship_boot_nofreq(data::PopData, method::Function, iterations::Int, interval::Vector{Float64} = [0.025, 0.975])
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    nloc = size(locmtx, 2)
    idxrange = 1:nloc
    result = NamedArray{Float64}(n, n)
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    bresult = Matrix{Float64}(undef, n,n)
    b_sdev = similar(bresult)
    b_ci = Matrix{Vector{Float64}}(undef, n,n)
    pbar = ProgressBar(;refresh_rate=90, transient = true)
    job = addjob!(pbar; description= "Kinship: ", N= Int64((n * (n-1))/2))
    start!(pbar)
    @inbounds @sync for i in 1:n-1
        Base.Threads.@spawn begin
            @inbounds v1 = view(locmtx,i,:)
            boot_idx = Vector{Int64}(undef, nloc)
            sizehint!(boot_idx, nloc)
            @inbounds for j in i+1:n
                @inbounds v2 = view(locmtx,j,:)
                @inbounds result[j,i] = result[i,j] = method(v1, v2)
                bootstats = Series(Variance(), Quantile(interval))
                k = 1
                while k <= iterations
                    k += 1
                    @inbounds boot_idx .= rand(idxrange, nloc)
                    v1b = @inbounds view(locmtx, i,boot_idx)
                    v2b = @inbounds view(locmtx, j, boot_idx)
                    b_est = method(v1b, v2b)
                    isnan(b_est) ? continue : fit!(bootstats, b_est)
                end
                @inbounds bresult[i,j] = bootstats.stats[1].μ
                @inbounds b_sdev[i,j] = sqrt(bootstats.stats[1].σ2)
                @inbounds b_ci[i,j] = value(bootstats.stats[2])
                progress.update!(job)
            end
        end
    end
    stop!(pbar)
    ci = uppertri2vec(b_ci)
    cilow = getindex.(ci, 1)
    cihi = getindex.(ci,2)
    out = kinshiptotable(result, Symbol(method))
    insertcols!(out, :bootmean => uppertri2vec(bresult), :std => uppertri2vec(b_sdev), :CI_lower => cilow, :CI_upper => cihi)
    return out    
end

function _kinship_boot_freq(data::PopData, method::Function, iterations::Int, interval::Vector{Float64} = [0.025, 0.975])
    locmtx = locimatrix(data)
    ids = samplenames(data)
    n = length(ids)
    nloc = size(locmtx, 2)
    idxrange = 1:nloc
    result = NamedArray{Float64}(n, n)
    setnames!(result, String.(ids),1)
    setnames!(result, String.(ids),2)
    bresult = Matrix{Float64}(undef, n,n)
    b_sdev = similar(bresult)
    b_ci = Matrix{Vector{Float64}}(undef, n,n)
    pbar = ProgressBar(;refresh_rate=90, transient = true)
    job = addjob!(pbar; description= "Kinship: ", N= Int64((n * (n-1))/2))
    start!(pbar)
    allelefrequencies = @inbounds Tuple(allelefreq(i) for i in eachcol(locmtx))
    @inbounds @sync for i in 1:n-1
        Base.Threads.@spawn begin
            @inbounds v1 = view(locmtx,i,:)
            boot_idx = Vector{Int64}(undef, nloc)
            sizehint!(boot_idx, nloc)
            @inbounds for j in i+1:n
                @inbounds v2 = view(locmtx,j,:)
                @inbounds result[j,i] = result[i,j] = method(v1, v2, allelefrequencies, n_samples = n)
                bootstats = Series(Variance(), Quantile(interval))
                k = 1
                while k <= iterations
                    k += 1
                    @inbounds boot_idx .= rand(idxrange, nloc)
                    v1b = @inbounds view(locmtx, i,boot_idx)
                    v2b = @inbounds view(locmtx, j, boot_idx)
                    b_est = method(v1b, v2b, allelefrequencies, n_samples = n)
                    isnan(b_est) ? continue : fit!(bootstats, b_est)
                end
                @inbounds bresult[i,j] = bootstats.stats[1].μ
                @inbounds b_sdev[i,j] = sqrt(bootstats.stats[1].σ2)
                @inbounds b_ci[i,j] = value(bootstats.stats[2])
                progress.update!(job)
            end
        end
    end
    stop!(pbar)
    ci = uppertri2vec(b_ci)
    cilow = getindex.(ci, 1)
    cihi = getindex.(ci,2)
    out = kinshiptotable(result, Symbol(method))
    insertcols!(out, :bootmean => uppertri2vec(bresult), :std => uppertri2vec(b_sdev), :CI_lower => cilow, :CI_upper => cihi)
    return out    
end

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

# TODO check math, should diag = 1?
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
    Sxy = 0.0
    S0 = 0.0
    @inbounds for i in eachindex(ind1)
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

# FIX _MORAN MATH
function _moran(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
    num = 0.0 ; denom = 0.0
    @inbounds for (allele, fq) in frqdict
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


function _Ritland(geno1::NTuple{2,T}, geno2::NTuple{2,T}, frqdict::Dict{T, Float64}) where T<:Union{Int16, Int8}
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


"""
    kinshiptotable(kinshipresults::T, methd::Symbol) where T<:NamedMatrix
Converts the `NamedMatrix` result from the non-bootstrapped `kinship()` results into a `DataFrame`.
The second positonal argument (`methd`) is the name of the value column (default: `kinship`). For
better analysis workflow, it would be useful to specify the method for this column, to
keep track of which estimator was used (e.g., `Blouin`, `LynchLi`, etc.)

**Example**
```julia`
julia> cats = @nancycats ; kin = kinship(cats, method = Moran) ;

julia> kinshiptotable(kin, :Moran)
22366×3 DataFrame
   Row │ sample1  sample2  Moran      
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
function kinshiptotable(kinshipresults::T, mthd::Symbol=:nothing) where T<:NamedMatrix
    n = size(kinshipresults,1)
    n == size(kinshipresults,2) || throw(DimensionMismatch("The input matrix must be symmetrical, but has size $(size(kinshipresults))"))
    ids = names(kinshipresults)[1]
    vals = uppertri2vec(kinshipresults)
    idpairs = pairwisepairs(ids)
    meth = mthd == :nothing ? :kinship : mthd 
    DataFrame(:sample1 => first.(idpairs), :sample2 => getindex.(idpairs, 2), meth => vals)
end


#TODO update kinshipposthoc for this new interface