export genepop

"""
    genepop(infile::String; kwargs...)
Load a Genepop format file into memory as a PopData object.
- `infile::String` : path to Genepop file

### Keyword Arguments
- `digits::Integer`: number of digits denoting each allele (default: `3`)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `diploid::Bool`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent::Bool`   : whether to print file information during import (default: `true`)

### File must follow standard Genepop formatting:
- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default `POP`) must delimit populations
- Sample name is immediately followed by a *comma*
- File is *tab or space delimted* (but not both!)

### Genepop file example:
wasp_hive.gen: Wasp populations in New York \n
Locus1\n
Locus2\n
Locus3\n
pop\n
Oneida_01,  250230  564568  110100\n
Oneida_02,  252238  568558  100120\n
Oneida_03,  254230  564558  090100\n
pop\n
Newcomb_01, 254230  564558  080100\n
Newcomb_02, 000230  564558  090080\n
Newcomb_03, 254230  000000  090100\n
Newcomb_04, 254230  564000  090120\n

## Example
```
waspsNY = genepop("wasp_hive.gen", digits = 3, popsep = "pop")
```
"""

infile = "data/data/gulfsharks.gen"

function genepop2(
    infile::String;
    digits::Int = 3,
    popsep::String = "POP",
    diploid::Bool = true,
    silent::Bool = false,
)
    # open the file as lines of strings to suss out loci names, pop idx, and popcounts
    gpop_readlines = readlines(infile)

    # find the #samples per population
    pop_idx = Vector{Int}()
    for (line, gtext) in enumerate(gpop_readlines)
        if gtext == popsep
            push!(pop_idx, line)
        end
    end
    length(pop_idx) == 0 &&
    error("No populations found in $infile using separator \"$popsep\". Please check the spelling and try again.")

    # create a theoretical place were the last popsep would be (preserve last population for counting)
    push!(pop_idx, countlines(infile) + 1)
    popcounts = (pop_idx[2:end] .- pop_idx[1:end-1]) .- 1

    # check for the delimiter type in the first sample record
    firstrecord = gpop_readlines[pop_idx[1]+1]
    if occursin("\t", firstrecord) & occursin(" ", firstrecord)
        error("$infile contains both tab and space delimiters. Please format the file so it uses either one or the other.")
    elseif occursin("\t", firstrecord)
        delim = "\t"
        delim_txt = "tab"
    elseif occursin(" ", firstrecord)
        delim = " "
        delim_txt = "space"
    else
        error("Please format $infile to be either tab or space delimited")
    end

    # check for horizontal formatting, where popsep would appear on the third line
    if pop_idx[1] <= 3

        gpop = open(infile, "r") ; readline(gpop)  # skip first line
        # second line should have all the loci
        locus_name_raw = readline(gpop) ; close(gpop)
        locinames = strip.(split(locus_name_raw |> join, ","))
        categorical(map(i -> replace!(i, "." => "_"), locinames), compress = true)
        horizontal = true
    else
        #standard  vertical formatting
        locinames = CSV.read(
            infile,
            header = 0,
            datarow = 2,
            copycols = false,
            limit = pop_idx[1] - 2,
            type = String,
        )
        locinames = categorical(locinames.Column1, compress = true)
        horizontal = false
    end

    if !silent
        @info "\n$(abspath(infile))\n$(delim_txt) delimiter detected\nhorizontal formatting: $(horizontal)\n$(sum(popcounts)) samples detected\n$(length(locinames)) loci detected"
    end

    # load in samples and genotypes
    coln = append!(["name"], locinames)
    geno_parse = CSV.read(
        infile,
        delim = delim,
        header = coln,
        datarow = pop_idx[1] + 1,
        comment = popsep,
        type = String
    )

    popnames = string.(collect(1:length(popcounts)))
    popnames = fill.(popnames,popcounts) |> Base.Iterators.flatten |> collect
    insertcols!(geno_parse, 2, :population => categorical(popnames, compress = true))
    geno_parse = DataFrames.stack(geno_parse, DataFrames.Not([:name, :population]))
    geno_parse.name = categorical(map(i -> replace(i, "," => ""), geno_parse.name), compress = true)

    #geno_type = determine_marker(infile, geno_parse, digits)
    geno_type = Int8
    geno_parse.value = map(i -> phase.(i, geno_type, digits), geno_parse.value)
    return geno_parse
    #=
    loci_table = table(
                    (name = geno_parse.name,
                    population = geno_parse.population,
                    locus = geno_parse.variable,
                    genotype = geno_parse.value),
                    pkey = [:name, :population]
                )

                =#    #=
    # add population names to genotype table
    loci_table = table((
        name = categorical(name, compress = true),
        population = categorical(popid_array, compress = true),
        locus = categorical(loci, compress = true),
        genotype = genotype,
    ), pkey = :locus)

    # make sure levels are sorted by order of appearance
    levels!(loci_table.columns.locus, unique(loci_table.columns.locus))
    levels!(loci_table.columns.name, unique(loci_table.columns.name))

    ploidy = (@groupby loci_table :name {
        ploidy = find_ploidy(:genotype),
    }).columns.ploidy

    # take a piece of the genotype table out and create a new table with the ploidy
    sample_table = table((
        name = levels(loci_table.columns.name),
        population =
                [fill(i, j) for (i, j) in zip(popnames, popcounts)] |>
                Base.Iterators.flatten |> collect,
        ploidy = ploidy,
        latitude =
                fill(missing, length(levels(loci_table.columns.name))) |>
                Vector{Union{Missing,Float32}},
        longitude =
                fill(missing, length(levels(loci_table.columns.name))) |>
                Vector{Union{Missing,Float32}},
    ), pkey= :population)

    PopData(sample_table, loci_table)
=#
end
