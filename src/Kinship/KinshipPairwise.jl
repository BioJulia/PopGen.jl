#=
function testkin(x, y)
    sym = Symbol(y)
    old = kinship(x, method = y)
    new = kinshiptotable(kinship_new(x, method = y))
    old[:, sym] = round.(old[:, sym], digits = 5)
    new[:, :kinship] = round.(new.kinship, digits = 5)
    if old[:, sym] == new.kinship
        return true
    else
        insertcols!(old, :new => new.kinship)
        return old[old[:, sym] .!= old.new, :]
    end
end
=#
function kinship(data::PopData; method::Function, iterations::Int = 0, interval::Vector{Float64} = [0.025, 0.975])
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

function kinship(data::PopData, samplenames::AbstractVector{T}; method::Function, iterations::Int = 0, interval::Vector{Float64} = [0.025, 0.975]) where T<:AbstractString
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
                update!(job)
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
                update!(job)
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


"""
    kinshiptotable(kinshipresults::T, methd::Symbol) where T<:NamedMatrix
Converts the `NamedMatrix` result from the non-bootstrapped `kinship()` results into a `DataFrame`.
The second positonal argument (`methd`) is the name of the value column (default: `kinship`). For
better analysis workflow, it would be useful to specify the method for this column, to
keep track of which estimator was used (e.g., `Blouin`, `LynchLi`, etc.)

**Example**
```julia
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