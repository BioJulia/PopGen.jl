using NamedTupleTools, JuliaDB

geno_type = Int8
ndigits = 3
function phase(loc::String, type::DataType, digit::Int)
    phasedlocus = parse.(type,
        [join(i) for i in Iterators.partition(loc, digit)]) |> sort |> Tuple
end

function genepop(
    infile::String;
    ndigits::Int = 3,
    popsep::Any = "POP",
    marker::String = "snp",
)
    println("\n", "Input File : ", abspath(infile))
    if lowercase(marker) == "snp"
        geno_type = Int8
    else
        geno_type = Int16
    end
    gpop = split(open(readlines, infile)[2:end], popsep)
    #gpop = split(open(readlines, "C:/Users/pdime/Desktop/PopGen.jl/data/data/gulfsharks.gen")[2:end], "POP")
    gpop = split(open(readlines, "/home/pdimens/PopGen.jl/data/data/gulfsharks.gen")[2:end], "POP")
    @info "\n$(abspath(infile)) \n$(length(gpop)-1) populations detected using \"$popsep\" seperator \n$(length(gpop[1])) loci detected \n$(sum(length.(gpop[2:end]))) samples detected"
    if length(gpop[1]) == 1     # loci horizontally stacked.
        locinames = strip.(split(gpop[1] |> join, ",") |> Vector{String})
        replace!(locinames, "." => "_")
    else                        # loci vertically stacked
        locinames = replace(gpop[1], "." => "_")
    end
    table_fields = append!(
                        ["name",
                        "population",
                        "ploidy",
                        "longitude",
                        "latitude"],
                        locinames
                        )
    sample_meta = [Vector{String}(),
                   Vector{Int}(),
                   Vector{Int8}(),
                   Vector{Union{Missing,Float32}}(),
                   Vector{Union{Missing, Float32}}()]
    table_types = append!(sample_meta, fill(Vector{Union{Missing, Tuple}}(), length(locinames)))
    d = namedtuple(table_fields, table_types)
    #d = Dict(Symbol(i) => Vector{Union{Missing, Tuple}}() for i in locinames)
    #popid = Vector{Int}()
    #indnames = Vector{String}()
    #ploidy = Vector{Int8}()
    for i = 2:length(gpop)
        append!(d.population, fill(i - 1, length(gpop[i])))
        for j = 1:length(gpop[i])
            #println(Base.Threads.threadid()) ########
            #phasedloci = Vector{Union{Missing, Tuple}}()
            samplerow = split(strip(gpop[i][j]), r"\,\s|\s|\t") |>Vector{String}
            push!(d.name, popfirst!(samplerow))
            # just in case -9 = missing
            replace!(samplerow, "-9" => "0"^ndigits)
            # phase the loci into tuples
            for (locname, locus) in zip(locinames,samplerow)
                phasedlocus = parse.(
                                geno_type,
                                [join(i) for i in Iterators.partition(locus, ndigits)]) |> sort |> Tuple
                push!(d[Symbol(locname)], phasedlocus)
            end

            #=phasedloci = map(x -> phase(x, geno_type, ndigits), samplerow)
            return
        end
    end
    #= infer ploidy by taking the mode of # alleles per locus per indiv for
       the first 20% of loci =#
       #map(length., [d[i] for i in locinames[1:end÷5]])
    #ploidy = mode(length.(phasedloci[1:end÷5] |> skipmissing))
#TODO
    loc2 = Symbol.(locinames)
    #genos = (; zip(loc2, Vector{Union{Tuple,Missing}}(d[i]) for i in loc2)...) # convert to a giant named tuple
    ploidy = length.(d[locinames[1]])   # lazy finding of ploidy from single locus
    for (loc, ploid) in zip(locinames, ploidy)
        miss_geno = fill(0, ploid) |> Tuple
        msat_miss_geno = ("0")
        replace!(d[loc], miss_geno => missing)
        replace!(d[loc], msat_miss_geno => missing)
    end
    sampleinfo = (name = indnames,
                population = string.(popid),
    #ploidy = Int8.(ploidy),
                longitude = fill(missing, length(indnames)) |> Vector{Union{Missing, T} where T <: Real},
                latitude = fill(missing, length(indnames)) |> Vector{Union{Missing, T} where T <: Real},
    )
    JDB_table = table(sampleinfo)

 length(indnames)) |> Vector{Union{Missing, T} where T <: Real},
        (i = d[i] for i in locinames)
    ))
    #PopObj(samples_df, loci_df)

end
