

function genepop(
    infile::String;
    ndigits::Int64 = 3,
    popsep::Any = "POP",
    marker = "snp",
)

    if lowercase(marker) == "snp"
        geno_type = Int8
    else
        geno_type = Int16
    end
    gpop = split(open(readlines, infile)[2:end], popsep)
    gpop = split(open(readlines, "C:/Users/pdime/Desktop/PopGen.jl/data/data/gulfsharks.gen")[2:end], "POP")
    gpop = split(open(readlines, "/home/pdimens/PopGen.jl/data/data/gulfsharks.gen")[2:end], "POP")
    @info "\n$(abspath(infile)) \n$(length(gpop)-1) populations detected using \"$popsep\" seperator \n$(length(gpop[1])) loci detected \n$(sum(length.(gpop[2:end]))) samples detected"

    if length(gpop[1]) == 1     # loci horizontally stacked
        locinames = strip.(split(gpop[1] |> join, ",") |> Vector{String})
        replace!(locinames, "." => "_")
    else                        # loci vertically stacked
        locinames = replace(gpop[1], "." => "_")
    end
    d = Dict()
    popid = Vector{String}()
    indnames = Vector{String}()
    for i = 2:length(gpop)
        append!(popid, fill(i - 1, length(gpop[i])))
        for j = 1:length(gpop[i])
            phasedloci = []
            samplerow = split(strip(gpop[i][j]), r"\,|\s|\t") |> Vector{String}
            push!(indnames, samplerow[1])
            unphasedloci = samplerow[3:end]
            # just in case -9 = missing
            replace!(unphasedloci, "-9" => "0"^ndigits)
            # phase the loci into tuples
            for locus in unphasedloci
                phasedlocus = parse.(
                    geno_type,
                    [join(i) for i in Iterators.partition(locus, ndigits)],
                ) |> sort |> Tuple
                push!(phasedloci, phasedlocus)
            end
            for (loc, geno) in zip(Symbol.(locinames), phasedloci)
                push!(d[loc], geno)
            end
        end
    end

#TODO
    loc2 = Symbol.(locinames)
    genos = (; zip(loc2, Array{Union{Tuple,Missing},1}(d[i]) for i in loc2)...) # convert to a giant named tuple
    ploidy = length.(d[locinames[1]])   # lazy finding of ploidy from single locus
    for (loc, ploid) in zip(locinames, ploidy)
        miss_geno = fill(0, ploid) |> Tuple
        msat_miss_geno = ("0")
        replace!(d[loc], miss_geno => missing)
        replace!(d[loc], msat_miss_geno => missing)
    end
    # typesafe genotype DataFrame
    #loci_df = DataFrame([Symbol(i) => Array{Union{Tuple,Missing},1}(d[i]) for i in locinames])
#=
    loci_table = JuliaDB.table([Array{Union{Tuple,Missing},1}(d[i]) for i in locinames], names = Symbol.(locinames))
    samples_df = DataFrame(
        name = string.(indnames),
        population = string.(popid),
        ploidy = Int8.(ploidy),
        longitude = fill(missing, length(indnames)),
        latitude = fill(missing, length(indnames)),
    )
    samples_table = JuliaDB.table(
        name = string.(indnames),
        population = string.(popid),
        ploidy = Int8.(ploidy),
        longitude = fill(missing, length(indnames)),
        latitude = fill(missing, length(indnames)),
    )
    PopObj(samples_df, loci_df)
=#
end
