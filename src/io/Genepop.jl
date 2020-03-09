#=
This file handles the import/export of Genepop file format
=#

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

## Example
```
waspsNY = genepop("wasp_hive.gen", digits = 3, popsep = "POP")
```
### Genepop file example:
Wasp populations in New York \n
Locus1\n
Locus2\n
Locus3\n
POP\n
Oneida_01,  250230  564568  110100\n
Oneida_02,  252238  568558  100120\n
Oneida_03,  254230  564558  090100\n
POP\n
Newcomb_01, 254230  564558  080100\n
Newcomb_02, 000230  564558  090080\n
Newcomb_03, 254230  000000  090100\n
Newcomb_04, 254230  564000  090120\n
"""
function genepop(
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
        map(i -> replace!(i, "." => "_"), locinames)
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
        locinames = String.(locinames.Column1)
        horizontal = false
    end

    if !silent
        @info "\n$(abspath(infile))\n$(delim_txt) delimiter detected\nhorizontal formatting: $(horizontal)\n$(sum(popcounts)) samples detected\n$(length(locinames)) loci detected"
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
            type = String,
        )
    end

    geno_type = determine_marker(infile, geno_parse, digits)

    # initiate with empty vectors for samplename, locus, and geno
    name = Vector{String}()
    loci = Vector{String}()
    if diploid == true
        geno_raw = Vector{Int32}()
    else
        geno_raw = Vector{String}()
    end

    @inbounds for samplerow in geno_parse
        vals = values(samplerow)
        corr_name = replace(vals[1], "," => "")
        append!(name, fill(corr_name, length(locinames)))
        append!(loci, locinames)
        append!(geno_raw, vals[2:end])
    end

    ## phase the genotypes
    if diploid == true
        genotype = map(i -> phase_dip.(i, geno_type, digits), geno_raw)
    else
        genotype = map(i -> phase.(i, geno_type, digits), geno_raw)
    end

    ## create population vector for loci table
    popnames = string.(collect(1:length(popcounts)))
    popid_array =
        [
            fill(i, j)
            for (i, j) in zip(popnames, popcounts .* length(locinames))
        ] |> Base.Iterators.flatten |> collect

    # add population names to genotype table
    loci_table = table((
        name = categorical(name, true),
        population = categorical(popid_array, true),
        locus = categorical(loci, true),
        genotype = genotype,
    ))

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
    ))

    PopData(sample_table, loci_table)
end
