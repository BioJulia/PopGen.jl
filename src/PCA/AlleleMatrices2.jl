
function _freqmatrix_mean2(data::PopData)
    counts = @inbounds _countmatrix(data) ./ data.sampleinfo.ploidy
    map(eachcol(counts)) do alcol
        colmean = mean([x for x in alcol if x >= 0])
        replace!(x -> x < 0 ? colmean : x, alcol)
    end
    return counts
end


function _freqmatrix_missing(data::PopData)
    replace(_countmatrix(data), -1 => missing) ./ data.sampleinfo.ploidy

    counts = Matrix{Union{Missing, Float64}}(_countmatrix(data))
    # iterate over rows (samples)
    @inbounds for (_sample,_counts) in enumerate(eachrow(counts))
        countidx = findall(!ismissing, _counts)
        # divide non-missing values in the row by the sample's ploidy and replace the original values
        @inbounds _counts[countidx] .= _counts[countidx] ./ data.sampleinfo.ploidy[_sample]
    end
    return counts
end


missing/2