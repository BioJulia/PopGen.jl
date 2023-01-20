function allelicdistance(x, y)
    map(zip(x,y)) do (i,j)
        mapreduce(+,i) do alle
            count(∉(alle), j)
        end
    end
end

function _missinglocusfilter(data::PopData, cutoff::Float64)
    miss = missingdata(data, by = "locus")
    #passing = findall(<=(cutoff), miss.percent)
    failing = findall(>=(cutoff), miss.percent)
    if isempty(failing)
        return data
    else
        @info "removing $(length(failing)) loci with >$(cutoff*100)% missing data"
        #return data[data.genodata.locus .∈ Ref(loci(data)[passing])]
        return data[data.genodata.locus .∉ Ref(loci(data)[failing])]
    end
end

# this is needed for pairwisefst dispatch
AMOVA(data::PopData; kwargs...) = amova(data, kwargs...)


function amova(data::PopData; hierarchy::String, missing_cutoff::Union{Nothing,Float64} = 0.05)
    # remove any whitespace and split by /
    groups = Symbol.(split(replace(hierarchy, " " => ""), r"/"))
    # pull out sample name and stratification columns
    # include error handling for missing columns using try...catch
    strata = try
        data.metadata.sampleinfo[:, [:name,groups...]]
    catch
        erroridx = groups .∉ Ref(propertynames(data.metadata.sampleinfo))
        missingcolumns = join(groups[erroridx], ", ")
        throw(ArgumentError("Columns {$missingcolumns} not found in the PopData.sampleinfo dataframe."))
    end
    # filter out loci with too much missingness
    data = isnothing(missing_cutoff) ? data : _missinglocusfilter(data, missing_cutoff)
    # Squared Eucledian distance matrix based on allele freqs
    #2 .* pairwise(SqEuclidean(), mtx, dims=1)
    distmtx = pairwise(SqEuclidean(), matrix(data), dims=1)

    group = first(groups)
    grpcol = @view strata[:,group]
    levels = unique(grpcol)
    # create vector of vectors containing indices for each unique group
    groupidx = map(x -> findall(==(x),grpcol), levels)
    #popstats = Dict{Symbol, Float64}[Dict{Symbol, Float64}() for i in 1:length(levels)]
    df_among = length(levels) - 1
    df_within = 0
    N = data.metadata.samples
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
    df_total = df_among + df_within
    SS_total = SS_within + SS_among
    σ²_within = SS_within / df_within
    n_c = (N - (reduce(+, length.(groupidx).^2) / N)) / df_among
    σ²_among = ((SS_among / df_among) - σ²_within) / n_c
    FST = σ²_among / (σ²_among + σ²_within)
    return AMOVAResult(
        ["Total", "Among $(group)s", "Within $(group)s"],
        [df_total, df_among, df_within],
        [SS_total, SS_among, SS_within],
        [SS_among / df_among, σ²_within],
        [σ²_among, σ²_within],
        FST
        )
end

#=
function _amova_within_ind(mtx::AbstractMatrix)
    N = size(mtx, 1)
    df_among = N - 1
    df_within = 0
    SS_within = 0.0
    SS_among = 0.0
    n = 1
    for i in 1:Ntot
        dw = distmtx[i, i]
        da = reduce(+, view(distmtx, j, setdiff(1:N, i)))
        #df_within += n - 1
        SS_within += dw / 2n
        SS_among += ((da + dw) / 2N) - (dw / 2n)
    end
    df_total = df_among
    SS_total = SS_within + SS_among
    σ²_within = SS_within / df_within
    n_c = (N - (reduce(+, length.(groupidx).^2) / N)) / df_among
    σ²_among = ((SS_among / df_among) - σ²_within) / n_c
    FST = σ²_among / (σ²_among + σ²_within)
    
end
=#

struct AMOVAResult
    source::Vector{String}
    df::Vector{Int}
    SS::Vector{Float64}
    MS::Vector{Float64}
    σ²::Vector{Float64}
    FST::Float64
end

# Borrowed heavily from StatsModels.jl to have consistent style with expected
# model outputs https://github.com/JuliaStats/StatsModels.jl/blob/6f19ecca344b4dc99c41d0e70976f306a3dd9e72/src/lrtest.jl#L130
function Base.show(io::IO, result::AMOVAResult)
    source = result.source
    df = result.df
    ss = result.SS
    ms = result.MS
    σ² = result.σ²
    nc = 6
    nr = length(result.source)
    outrows = Matrix{String}(undef, nr+1, nc)

    outrows[1, :] = ["Source", "DF", "SS", "MS", "σ²", "FST"]

    outrows[2, :] = [
        source[1],
        Printf.@sprintf("%.0d", first(df)),
        Printf.@sprintf("%.4f", first(ss)),
        Printf.@sprintf("%.4f", first(ms)),
        Printf.@sprintf("%.4f", first(σ²)),
        Printf.@sprintf("%.4f", result.FST)
    ]

    for i in 2:nr
        outrows[i+1, :] = [
            source[i],
            Printf.@sprintf("%.0d", df[i]),
            Printf.@sprintf("%.0d", ss[i]),
            i > length(ms) ? " " : Printf.@sprintf("%.4f", ms[i]),
            i > length(σ²) ? " " : Printf.@sprintf("%.4f", σ²[i]),
            " "
        ]
    end
    colwidths = length.(outrows)
    max_colwidths = [maximum(view(colwidths, :, i)) for i in 1:nc]
    totwidth = sum(max_colwidths) + 2*(nc-1)

    println(io, "Analysis of Molecular Variance")
    println(io, '─'^totwidth)

    for r in 1:nr+1
        for c in 1:nc
            cur_cell = outrows[r, c]
            cur_cell_len = length(cur_cell)

            padding = " "^(max_colwidths[c]-cur_cell_len)
            if c > 1
                padding = "  "*padding
            end

            print(io, padding)
            print(io, cur_cell)
        end
        print(io, "\n")
        r == 1 && println(io, '─'^totwidth)
    end
    print(io, '─'^totwidth)
end