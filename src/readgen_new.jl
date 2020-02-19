using CSV, JuliaDB, CategoricalArrays, JuliaDBMeta, BenchmarkTools

infile = "data/data/nancycats.gen"
infile2 = "data/data/gulfsharks.gen"

@inline function find_ploidy(genotypes::T where T<:SubArray)
    collect(skipmissing(genotypes))[1] |> length |> Int8
end

##### read_genepop long
function longen(
    infile::String;
    digits::Int = 3,
    popsep::String = "POP",
    marker::String = "snp",
    diploid::Bool = true,
    delim::String = "\t"
)
    # establish types for genotypes
    if lowercase(marker) == "snp"
        geno_type = Int8
    else
        geno_type = Int16
    end

    # find the #samples per population
    pop_idx = Vector{Int}()
    for (line, gtext) in enumerate(readlines(infile))
        if gtext == popsep
            push!(pop_idx, line)
        end
    end
    length(pop_idx) == 0 && error("No populations found in $infile using separator \"$popsep\". Please check the spelling and try again.")
    # create a theoretical place were another POP would be (preserve last population for counting)
    push!(pop_idx, countlines(infile) + 1)
    popcounts = (pop_idx[2:end] .- pop_idx[1:end-1]) .- 1

    # check for horizontal formatting, where popsep would appear on the third line
    if pop_idx[1] <= 3
        gpop = open(infile, "r")
        # skip first line
        readline(gpop)
        # second line should have all the loci
        locus_name_raw = readline(gpop)
        close(gpop)
        locinames = strip.(split(locus_name_raw |> join, ","))
        map(i -> replace!(i, "." => "_"), locinames)
    else
        #standard  vertical formatting
        locinames = CSV.read(
                    infile,
                    header = 0,
                    datarow = 2,
                    copycols = false,
                    limit = pop_idx[1]-2,
                    type = String
        )
        locinames = String.(locinames.Column1)
    end

    # load in samples and genotypes
    coln = append!(["name"], locinames)
    if diploid == true
        geno_parse = CSV.File(
                infile,
                delim = delim,
                header = coln,
                datarow = pop_idx[1] + 1,
                comment = popsep,
        )
    else
        geno_parse = CSV.File(
                    infile,
                    delim = delim,
                    header = coln,
                    datarow = pop_idx[1] + 1,
                    comment = popsep,
                    type = String
        )
    end
    # initiate with empty table
    outtable = table((
        name=Vector{String}(),
        locus=Vector{String}(),
        genotype=Vector{Union{Missing, NTuple{N, geno_type}} where N}()
        ))

    for samplerow in geno_parse
        vals = values(samplerow)
        tmp = table((
                name = fill(vals[1],length(locinames)),
                locus = locinames,
                genotype = vals[2:end]
                ))
        outtable = merge(outtable, tmp)
    end

    ## create population column
    popnames = string.(collect(1:length(popcounts)))
    popid_array = [fill(i, j) for (i,j) in zip(popnames,popcounts .* length(locinames))] |>
        Iterators.flatten |>
        collect |>
        CategoricalArray

    # add population names to genotype table
    outtable = insertcolsafter(outtable, 1, :population => (popid_array))

    # do some table formatting
    ## remove commas from names, and format types for names and loci
    outtable = transform(outtable,
                    :name => replace.(outtable.columns.name, "," => "") |> CategoricalArray,
                    :locus => CategoricalArray(outtable.columns.locus),
                    )

    ## phase the genotypes
    if diploid == true
        outtable = transform(outtable, :genotype => map(i -> phase_dip.(i.genotype, geno_type, digits), outtable))
    else
        outtable = transform(outtable, :genotype => map(i -> phase.(i.genotype, geno_type, digits), outtable))
    end

    # take a piece of the genotype table out and create a new table with the ploidy
    sample_table = @groupby outtable (:name, :population) {ploidy = find_ploidy(:genotype)}

    # add missing long and lat columns
    sample_table = insertcolsafter(sample_table, 3, :longitude => fill(missing, length(sample_table.columns.name)) |> Vector{Union{Missing, Float32}},
                        :latitude => fill(missing, length(sample_table.columns.name)) |> Vector{Union{Missing, Float32}})

    return sample_table
    sample_table = transform(sample_table, :name => string.(sample_table.columns.name),  :population => string.(sample_table.columns.population))#, :population => string.(:population))

    PopObj(sample_table, outtable)
end
