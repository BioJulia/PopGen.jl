function _onehot_biallelic(data::PopData)
    data.metadata.biallelic == false && throw(ArgumentError("Data must be biallelic. If practical to do so, use dropmultiallelic() to remove non-biallelic loci"))
    gmtx = locimatrix(data)
    mapreduce(hcat, eachcol(gmtx)) do locus
        alle = uniquebialleles(locus)
        d = Dict{Union{Missing, NTuple{2, eltype(alle)}}, Int8}(
            (alle[1], alle[1]) => Int8(0),
            (alle[1], alle[2]) => Int8(1),
            (alle[2], alle[2]) => Int8(2),
            missing => Int8(-1)
        )
        Int8[d[i] for i in locus]
    end
end

function randomforest(data::PopData)
    missing ∈ data.genodata.genotype && throw(error("Unfortunately, random forest analysis cannot work with missing data. You may try to filter out loci and/or samples with missing data or use a genotype inputation method (e.g. fastPhase, BEAGLE)."))
    input = _onehot_biallelic(data)


end
