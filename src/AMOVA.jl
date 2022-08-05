function allelicdistance(x, y)
    map(zip(x,y)) do (i,j)
        mapreduce(+,i) do alle
            count(∉(alle), j)
        end
    end
end

function _missinglocusfilter(data::PopData, cutoff::Float64)
    miss = missingdata(data, by = "locus")
    passing = findall(<=(cutoff), miss.percent)
    @info "removing $(data.metadata.loci - length(passing)) loci with >$cutoff percent missingness"
    data[data.genodata.locus .∈ Ref(loci(data)[passing])]
end

function amova(data::PopData; hierarchy::String, missing_cutoff::Float64 = 0.05)
    # remove any whitespace and split by /
    groups = Symbol.(split(replace(hierarchy, " " => ""), r"/"))
    # pull out sample name and stratification columns
    # include error handling for missing columns using try...catch
    strata = try
        data.metadata.sampleinfo[:, [:name,groups...]]
    catch
        erroridx = groups .∉ Ref(propertynames(data.metadata.sampleinfo))
        missingcolumns = join(groups[erroridx], ", ")
        throw(ArgumentError("Columns {$missingcolumns} not found in metadata.sampleinfo dataframe."))
    end
    # filter out loci with too much missingness
    newdat = _missinglocusfilter(data, missing_cutoff)
    # matrix of allele frequencies
    mtx = matrix(newdat)  
    # Squared Eucledian distance matrix based on allele freqs
    #2 .* pairwise(SqEuclidean(), mtx, dims=1)
    distmtx = pairwise(SqEuclidean(), mtx, dims=1)

    group = groups[1]
    grpcol = @view strata[:,group]
    levels = unique(grpcol)
    # create vector of vectors containing indices for each unique group
    groupidx = map(x -> findall(==(x),grpcol), levels)
    #popstats = Dict{Symbol, Float64}[Dict{Symbol, Float64}() for i in 1:length(levels)]
    N = length(levels) - 1
    df_within = 0
    SS_within = 0.0
    SS_among = 0.0
    for j in groupidx
        n = length(j)
        dw = reduce(+, view(distmtx, j, j))
        da = reduce(+, view(distmtx, j, setdiff(1:data.metadata.samples, j)))
        df_within += n - 1
        SS_within += dw / 2n
        SS_among += ((da + dw) / 2N) - (dw / 2n)
    end
    df_among = N - 1
    df_total = df_among + df_within
    SS_total = SS_within + SS_among
    σ²_within = SS_within / df_within
    n_c = (N - (reduce(+, length.(groupidx).^2) / N)) / df_among
    σ²_among = ((SS_among / df_among) - σ²_within) / n_c
    FST = σ²_among / (σ²_among + σ²_within)
end
