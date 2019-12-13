"""
    phase(loc::String, type::DataType, digit::Int)
Takes a String of numbers returns a typed locus appropriate for PopGen.jl as used in the
`genepop` and `csv` file parsers. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

## Examples
`ph_locus = phase("128114", Int16, 3)`

`map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])`

`[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]`
"""
function phase(loc::String, type::DataType, digit::Int)
    phased = [join(k) for k in Iterators.partition(loc, digit)]
    typed = parse.(type, phased) |> sort!
    tupled = Tuple(typed)
    iszero(tupled |> collect) && return missing
    return tupled
end

### CSV parsing ###

"""
    csv(infile::String; delim::Union{Char,String,Regex}, digits::Int64 = 3)
Load a CSV-type file into memory as a PopObj object. *There should be no empty cells in your
CSV file*
### Arguments
- `infile` : path to CSV file
- `delim` : delimiter characters. By default auto-parses comma, space, or tab (use ONLY one!)
- `digits` : number of digits denoting each allele (default = 3)
- `marker` : "snp" (default) or "msat"

### File formatting:
- First row is column names in this order:
    1. name
    2. population
    3. longitude
    4. latitude
    5. locus_1_name
    6. locus_2_name
    7. etc...

#### Formatting location data
If location data is missing for a sample (which is ok!), make sure the value is written
as *0*, otherwise there will be transcription errors!

## Example
`lizardsCA = Read.csv("CA_lizards.csv", digits = 3);`

### Formatting example:
name,population,long,lat,Locus1,Locus2,Locus3   \n
sierra_01,mountain,11.11,-22.22,001001,002002,001001   \n
sierra_02,mountain,11.12,-22.21,001001,001001,001002   \n
snbarb_03,coast,0,0,001001,001001,001002 \n
snbarb_02,coast,11.14,-22.24,001001,001001,001001 \n
snbarb_03,coast,11.15,0,001002,001001,001001 \n
"""
function csv(
    infile::String;
    delim::Union{Char,String,Regex} = r"\,|\t|\s",
    digits::Int = 3,
    marker = "snp",
)

    # instantiate empty vectors
    indnames = Vector{String}()
    popid = Vector{Union{String,Int32}}()
    ploidy = Vector{Int8}()
    locx = Vector{Union{Missing, Float32}}()
    locy = Vector{Union{Missing, Float32}}()
    sample_geno_array = Array{Vector{Union{Missing, Tuple}},1}()

    # establish types for genotypes
    if lowercase(marker) == "snp"
        geno_type = Int8
    else
        geno_type = Int16
    end

    file = open(infile, "r")

    # get header information
    col_names = map(join, split(readline(file), delim))
    replace!(col_names, "." => "_")
    locinames = col_names[5:end]

    @inbounds for ln in eachline(file)
        tmp = map(join, split(ln, delim))
        # just in case -9 = missing
        replace!(tmp, "-9" => "0"^digits)
        push!(indnames, popfirst!(tmp))
        push!(popid, popfirst!(tmp))
        push!(locx, parse(Float32, popfirst!(tmp)))
        push!(locy, parse(Float32, popfirst!(tmp)))

        # phase the loci
        phasedloci = map(i -> phase(i, geno_type, digits), tmp)
        push!(sample_geno_array, phasedloci)

        # get the ploidy via # of alleles in the first non-missing locus
        sample_ploidy = (skipmissing(phasedloci) |> collect)[1]  |> length
        push!(ploidy, sample_ploidy)

    end
    close(file)

    @info "\n$(abspath(infile)) \n$(length(indnames)) samples detected \n$(length(unique(popid))) populations detected \n$(length(locinames)) loci detected"

    # convert the entire thing into a matrix
    geno_matrix = hcat(sample_geno_array...)

    # convert is back into an array of arrays, transposed across loci
    # instead of samples
    loci_array = collect.(eachrow(geno_matrix))

    # replace 0.0 with missing in location
    replace!(locx, 0.0 => missing)
    replace!(locy, 0.0 => missing)

    # genotype DataFrame
    loci_df = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), loci_array)])

    samples_df = DataFrame(
        name = indnames,
        population = popid,
        ploidy = ploidy,
        longitude = locx,
        latitude = locy,
    )

    return PopObj(samples_df, loci_df)
end

### GenePop parsing ###
"""
    genepop(infile::String; digits::Int = 3, popsep::String = "POP", marker::String = "snp")
Load a Genepop format file into memory as a PopObj object.
### Arguments
- `infile` : path to Genepop file
- `digits` : number of digits denoting each allele
- `popsep` : word that separates populations in `infile` (default: "POP")
- `marker` : "snp" (default) or "msat"
### File must follow standard Genepop formatting:
- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default "POP") must delimit populations
- Sample name is followed by a comma and space or tab (", ")
- File is tab or space delimted

## Example

 `waspsNY = genepop("wasp_hive.gen", digits = 3, popsep = "POP")`

### Genepop file example:
Wasp populations in New York \n
Locus1\n
Locus2\n
Locus3\n
POP\n
Oneida_01,  250230 564568 110100\n
Oneida_02,  252238 568558 100120\n
Oneida_03,  254230 564558 090100\n
POP\n
Newcomb_01,  254230 564558 080100\n
Newcomb_02,  000230 564558 090080\n
Newcomb_03,  254230 000000 090100\n
Newcomb_04,  254230 564000 090120\n
"""
function genepop(
    infile::String;
    digits::Int = 3,
    popsep::String = "POP",
    marker::String = "snp",
)
    if lowercase(marker) == "snp"
        geno_type = Int8
    else
        geno_type = Int16
    end

    gpop = split(open(readlines, infile)[2:end], popsep)

    @info "\n$(abspath(infile)) \n$(sum(length.(gpop[2:end]))) samples detected \n$(length(gpop)-1) populations detected using \"$popsep\" seperator \n$(length(gpop[1])) loci detected"

    if length(gpop[1]) == 1     # loci horizontally stacked.
        locinames = strip.(split(gpop[1] |> join, ",") |> Vector{String})
        replace!(locinames, "." => "_")
    else                        # loci vertically stacked
        locinames = replace(gpop[1], "." => "_")
    end

    sample_meta_names = [:name, :population, :ploidy, :latitude, :longitude]
    sample_meta_array = [
        Vector{String}(),
        Vector{Union{String,Int}}(),
        Vector{Int8}(),
        Vector{Union{Missing,Float32}}(),
        Vector{Union{Missing, Float32}}()
    ]

    sample_geno_array = Vector{Vector{Union{Missing, Tuple}}}()

    for i = 2:length(gpop)
        for j in gpop[i]
            samplerow = map(join, split(strip(j), r"\,\t|\,\s|\s|\t"))

            # just in case missing genotypes are coded as -9 (msats)
            replace!(samplerow, "-9" => "0"^digits)

            # sample name and population
            push!(sample_meta_array[1], popfirst!(samplerow))
            push!(sample_meta_array[2], i-1)

            # phase the loci into tuples
            converted_row = map(k -> phase(k, geno_type, digits), samplerow)
            push!(sample_geno_array, converted_row)

            # get the ploidy via # of alleles in the first non-missing locus
            push!(sample_meta_array[3], length((skipmissing(converted_row) |> collect)[1]))

            # long and lat
            push!(sample_meta_array[4], missing) ; push!(sample_meta_array[5], missing)
        end
    end

    # convert the entire thing into a matrix
    geno_matrix = hcat(sample_geno_array...)

    # convert back into an array of arrays of loci instead of samples
    loci_array = collect.(eachrow(geno_matrix))

    # create the two dataframes necessary for a PopObj
    samples = DataFrame([j => k for (j,k) in zip(sample_meta_names, sample_meta_array)])
    loci = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), loci_array)])

    return PopObj(samples, loci)
end

### VCF parsing ###

"""
    vcf(infile::String)
Load a VCF file into memory as a PopObj object. Population and [optional]
location information need to be provided separately.
- `infile` : path to VCF file
"""
function vcf(infile::String)
    vcf_file = VCF.Reader(open(infile, "r"))

    # get sample names from header
    sample_names = header(vcf_file).sampleID

    # fill in pop/lat/long with missing
    population = fill(missing, length(sample_names))
    lat = fill(missing, length(sample_names))
    long = fill(missing, length(sample_names))

    # get loci names
    locinames = Vector{String}()

    ## array of genotypes
    locus_array = Vector{Vector{Union{Missing, Tuple}}}()

    # get genotypes
    for record in vcf_file
        # fix locus names to be syntax-safe
        chr_safe = replace(VCF.chrom(record), "." => "_")
        chr_safer = replace(chr_safe, "|" => "_")
        pos = VCF.pos(record) |> string
        push!(locinames, chr_safer*"_"*pos)

        # get the genotype information
        geno_raw = [split(i, ('/', '|')) for i in VCF.genotype(record, :, "GT")] |> sort

        # change missing data "." to "-1"
        geno_corr_miss = map(i -> replace(i, "." => "-1"), geno_raw)

        # convert everything to an integer
        geno_int = map(i -> parse.(Int8, i), geno_corr_miss)

        # add 1 to shift genos so 0 is 1 and -1 is 0
        geno_shift = map(i -> i .+ Int8(1), geno_int)
        geno_final = [replace(i, 0 => missing) for i in geno_shift]
        geno_tuple = [Tuple(i) for i in geno_final]
        push!(locus_array, geno_tuple)
    end

    # intelligently scan for ploidy
    ploidy = Vector{Int8}()
    @inbounds for i in 1:length(sample_names)
        @inbounds for j in 1:length(locinames)
            locus_array[j][i] === missing && continue   # if missing, go to next locus
            ploid = length(locus_array[j][i])   # if not, get the # of alleles
            push!(ploidy, ploid)    # push that to the ploidy vector
            break   # break out of the loop and begin next sample
        end
    end

    @info "\n$(abspath(infile)) \n$(length(sample_names)) samples detected \n$(length(locinames)) loci detected"

    # create loci dataframe
    loci_df = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), locus_array)])

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = population,
        ploidy = ploidy,
        latitude = lat,
        longitude = long
    )
    PopObj(samples_df, loci_df)
end
