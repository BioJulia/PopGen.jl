
function pwsimilarhelper(x::Int8, y::Int8)
    (x == -1 | y == -1) && return missing
    (x == 0 & y == 0) && return false
    x == y && return true
    (x != 0 & y != 0) && return true

end

function pwsimilar(data::PopData)
    allelemtx = matrix(data, "count")
end