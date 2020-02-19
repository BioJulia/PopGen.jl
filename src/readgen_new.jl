using CSV, JuliaDB, CategoricalArrays, JuliaDBMeta, BenchmarkTools

infile = "data/data/nancycats.gen"
infile2 = "data/data/gulfsharks.gen"

@inline function find_ploidy(genotypes::T where T<:SubArray)
    for i in genotypes
        i !== missing && return Int8(length(i))
    end
end

"""
    phase(loc::String, type::DataType, digit::Int)
Takes a String of numbers returns a typed locus appropriate for PopGen.jl as used in the
`genepop` and `csv` file parsers. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

## Examples
```
ph_locus = phase("128114", Int16, 3)
map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
"""
@inline function phase(loc::T, type::DataType, digit::Int) where T<:AbstractString
    loc == "-9" || loc == "0"^length(loc) && return missing
    phased = map(i -> parse(type, join(i)), Iterators.partition(loc, digit))
    sort!(phased)
    tupled = Tuple(phased)
    return tupled
end

"""
    phase_dip(loc::String, type::DataType, digit::Int)
A diploid-optimized variant of `phase()` that uses integer division to split the alleles
of a locus into a tuple of `type` Type. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

## Examples
```
ph_locus = phase("128114", Int16, 3)
map(i -> phase_dip(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase_dip(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
"""
@inline function phase_dip(loc::T, type::DataType, digit::Int) where T<:Signed
    loc == -9 || iszero(loc) && return missing
    units = 10^digit
    allele1 = loc ÷ units |> type
    allele2 = loc % units |> type
    return [allele1, allele2] |> sort |> Tuple
end

@inline function phase_dip(loc::Missing, type::DataType, digit::Int)
    return missing
end

@inline function phase_dip(loc::String, type::DataType, digit::Int)
    loc == "-9" || loc == "0"^length(loc) && return missing
    newloc = parse(Int32, loc)
    units = 10^digit
    allele1 = newloc ÷ units |> type
    allele2 = newloc % units |> type
    return [allele1, allele2] |> sort |> Tuple
end



##### read_genepop long
"""
    genepop(infile::String; kwargs...)
Load a Genepop format file into memory as a PopObj object.
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

    # check for the marker type and assign Int8 or Int16 depending on max value
    # there's no real reason 40 was chosen other than being a reasonable buffer
    # for edge cases and missing data
    # since the genotypes are sorted anyway, we can look at the last (largest) value
    sample_geno = phase.(split(firstrecord, delim)[2:end], Int16, digits)
    if maximum(map(i -> i[end], sample_geno |> skipmissing)) < 40
        geno_type = Int8    # is a snp
    else
        geno_type = Int16   # is a msat
    end

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

    PopObj(sample_table, loci_table)
end
