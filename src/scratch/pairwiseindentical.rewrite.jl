function pwhelper(x::T,y::U)::Float64 where T<:AbstractArray where U<:AbstractArray
    mean(skipmissing(x .== y))
end

function pw2(data::PopData)
    locmtx = locimatrix(data)
    ids = unique(data.genodata.name)
    vecs = [i for i in eachrow(locmtx)]
    out = NamedArray(pwhelper.(vecs, permutedims(vecs)))
    setnames!(out, String.(ids),1)
    setnames!(out, String.(ids),2)
    return out
end

function pwsimiarhelper(x::Int8, y::Int8)
    (x == -1 | y == -1) && return missing
    (x == 0 & y == 0) && return false
    x == y && return true
    (x != 0 & y != 0) && return true

end

function pwsimilar(data::PopData)
    allelemtx = _allelematrix(data, by = "counts")
end