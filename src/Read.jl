

### Variant Call Format parsing ###

"""
    bcf(infile::String)
Load a BCF file into memory as a PopObj object. Population and [optional]
location information need to be provided separately.
- `infile` : path to BCF file
"""
function bcf(infile::String)
    bcf_file = BCF.Reader(open(infile, "r"))

    # get sample names from header
    sample_names = header(bcf_file).sampleID

    # fill in pop/lat/long with missing
    population = fill("missing", length(sample_names))
    loc_xy = Vector{Union{Missing, Float32}}()
    append!(loc_xy, fill(missing, length(sample_names)))

    # get loci names
    locinames = Vector{String}()

    ## array of genotypes
    locus_array = Vector{Vector{Union{Missing, Tuple{Vararg}}}}()

    # get genotypes
    for record in bcf_file
        # fix locus names to be syntax-safe
        chr_safe = replace(BCF.chrom(record), r"\.|\-|\=|\/" => "_")
        chr_safer = replace(chr_safe, "|" => "_")
        pos = BCF.pos(record) |> string
        push!(locinames, chr_safer*"_"*pos)

        # get the genotype information
        geno_raw = [split(i, ('/', '|')) for i in BCF.genotype(record, :, "GT")] |> sort

        # change missing data "." to "-1"
        geno_corr_miss = map(i -> replace(i, "." => "-1"), geno_raw)

        # convert everything to an integer
        geno_int = map(i -> parse.(Int8, i), geno_corr_miss)

        # add 1 to shift genos so 0 is 1 and -1 is 0 etc.
        geno_shift = map(i -> i .+ Int8(1), geno_int)
        geno_final = [replace(i, 0 => missing) for i in geno_shift]
        geno_tuple = [Tuple(i) for i in geno_final]
        push!(locus_array, geno_tuple)
    end
    close(vcf_file)

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

    @info "\n$(abspath(infile))
$(length(sample_names)) samples detected
population info must be added <---
$(length(locinames)) loci detected"

    # create loci dataframe
    loci_df = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), locus_array)])

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = population,
        ploidy = ploidy,
        latitude = loc_xy,
        longitude = loc_xy
    )
    PopObj(samples_df, loci_df)
end

### delimited parsing ###

"""
    delimited(infile::String; delim::Union{Char,String,Regex}, digits::Int64 = 3)
Load a delimited-type file into memory as a PopObj object. *There should be no empty cells
in your file*
### Arguments
- `infile` : path to file
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
```
lizardsCA = Read.delimited("CA_lizards.csv", digits = 3);
```

### Formatting example:
name,population,long,lat,Locus1,Locus2,Locus3   \n
sierra_01,mountain,11.11,-22.22,001001,002002,001001   \n
sierra_02,mountain,11.12,-22.21,001001,001001,001002   \n
snbarb_03,coast,0,0,001001,001001,001002 \n
snbarb_02,coast,11.14,-22.24,001001,001001,001001 \n
snbarb_03,coast,11.15,0,001002,001001,001001 \n
"""
function delimited(
    delim::Union{String} = ",",
    marker::String = "snp",
    digits::Int = 3,
    diploid::Bool = true
    )
    if lowercase(marker) == "snp"
        geno_type = Int8
    else
        geno_type = Int16
    end

    # read in CSV
    csv_in = CSV.read(infile, delim = delim, header = 1, copycols = true, type = String)

    # pull out loci
    loci_df = select(csv_in, Not(:1,:2,:3,:4))

    # pull out sample meta
    select!(csv_in, [:1,:2,:3,:4])
    rename!(csv_in, [:name, :population, :longitude, :latitude])
    csv_in.name = String.(csv_in.name)
    csv_in.population = String.(csv_in.population)

    # print some basic file info
    @info "\n$(abspath(infile))
$(size(csv_in,1)) samples detected
$(length(unique(csv_in.population))) populations detected
$(size(loci_df,2)) loci detected"

    # handle location data
    csv_in.longitude = map(csv_in.longitude) do val
        floatval = parse(Float64, val)
        if iszero(floatval) == true
            return missing
        else
            return floatval
        end
    end |> Vector{Union{Missing, Float64}}

    csv_in.latitude = map(csv_in.latitude) do val
        floatval = parse(Float64, val)
        if iszero(floatval) == true
            return missing
        else
            return floatval
        end
    end |> Vector{Union{Missing, Float64}}

    # phase the locus genotypes
    if diploid == true
        phased_df = map(eachcol(loci_df)) do locus
                phase_dip.(locus, geno_type, digits)
            end |> DataFrame

        ploidy = fill(2, length(csv_in.name))

    else
        phased_df = map(eachcol(loci_df)) do locus
                phase.(locus, geno_type, digits)
            end |> DataFrame

        ## ploidy finding
        ploidy = Vector{Int8}()

        @inbounds for i in 1:length(csv_in.name)
            @inbounds for j in eachcol(phased_df)
                j[i] === missing && continue   # if missing, go to next locus
                push!(ploidy, length(j[i]))   # if not, get ploidy and push to vector
                break   # break out of the loop and begin next sample
            end
        end
    end

    rename!(phased_df, names(loci_df))
    insertcols!(csv_in, 3, :ploidy => ploidy)

    return PopObj(csv_in, phased_df)

end

const csv = delimited

### GenePop parsing ###
"""
    genepop(infile::String; kwargs...)
Load a Genepop format file into memory as a PopObj object.
- `infile::String` : path to Genepop file

### Keyword Arguments
- `digits::Integer`: number of digits denoting each allele (default:  3)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `marker::String` : "snp" (default) or "msat"
- `delim::String` : file is tab `"\t"` (default) or space `" "` delimited
- `diploid::Bool`  : whether samples are diploid for parsing optimizations (default: true)
### File must follow standard Genepop formatting:
- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default "POP") must delimit populations
- Sample name is followed by a *comma and tab* (",  ")
- File is *tab or space delimted*

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
        #locinames = replace.(locinames, "." => "_")
    else
        #standard  vertical formatting
        locinames = CSV.read(infile,
                    header = 0,
                    datarow = 2,
                    copycols = false,
                    limit = pop_idx[1]-2,
                    type = String
                    ).Column1
        locinames = String.(locinames)
    end

    # load in samples and genotypes
    if diploid == true
        geno_df = CSV.read(infile,
                delim = delim,
                header = 0,
                datarow = pop_idx[1] + 1,
                copycols = false,
                comment = popsep,
                )
    else
        geno_df = CSV.read(infile,
                    delim = delim,
                    header = 0,
                    datarow = pop_idx[1] + 1,
                    copycols = false,
                    comment = popsep,
                    type = String
                    )
    end

    # split out sample name
    sample_df = select(geno_df, 1)
    # remove trailing commas from sample names
    sample_df[!,1] = replace.(string.(sample_df[!,1]), "," => "")
    # name the column
    rename!(sample_df, [:name])
    # get the number of samples
    n_samples = size(sample_df,1)

    # Count the number of popseps appearing in the file and name the pops 1...n in order of appearance
    popnames = string.(collect(1:length(popcounts)))
    popid_array = [fill(i, j) for (i,j) in zip(popnames,popcounts)]
    flat_popid = Iterators.flatten(popid_array) |> collect
    insertcols!(sample_df, 2, :population => flat_popid)
    # add remaining columns needed for PopObj.samples
    insertcols!(sample_df,3, :longitude => fill(missing, n_samples) |> Vector{Union{Missing, Float32}})
    insertcols!(sample_df,4, :latitude => fill(missing, n_samples) |> Vector{Union{Missing, Float32}})

    ## print input info
    @info "\n$(abspath(infile))
$(n_samples) samples detected
$(length(popnames)) populations detected using \"$popsep\" seperator
$(length(locinames)) loci detected"

    select!(geno_df, Not(1))
    # check columns to trim off ends if there were trailing whitespaces creating columns
    # full of only missing values
    while length(locinames) != size(geno_df,2)
        if ismissing.(geno_df[!, size(geno_df,2)]) |> all == true
            select!(geno_df, Not(size(geno_df,2)))
        end
    end

    # phasing the loci depending on diploid or not
    if diploid == true
        loc_df = map(eachcol(geno_df)) do locus
            map(i -> phase_dip(i, geno_type, digits), locus)
            #phase_dip.(locus, geno_type, digits)
        end |> DataFrame
        ploidy = fill(Int8(2), n_samples)
    else
        loc_df = map(eachcol(geno_df)) do locus
            map(i -> phase(i, geno_type, digits), locus)
            #phase.(locus, geno_type, digits)
        end |> DataFrame

        ## ploidy finding
        ploidy = Vector{Int8}()
        @inbounds for i in (1:n_samples)
            @inbounds for j in eachcol(loc_df)
                j[i] === missing && continue   # if missing, go to next locus
                push!(ploidy, length(j[i]))   # if not, get ploidy and push to vector
                break   # break out of the loop and begin next sample
            end
        end
    end

    # rename to locus dataframe columns with the loci names
    rename!(loc_df, locinames)
    # add the ploidy to the samples dataframe
    insertcols!(sample_df, 3, :ploidy => ploidy)

    return PopObj(sample_df, loc_df)
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
    population = fill("missing", length(sample_names))
    loc_xy = Vector{Union{Missing, Float32}}()
    append!(loc_xy, fill(missing, length(sample_names)))

    # get loci names
    locinames = Vector{String}()

    ## array of genotypes
    locus_array = Vector{Vector{Union{Missing, Tuple{Vararg}}}}()

    # get genotypes
    for record in vcf_file
        # fix locus names to be syntax-safe
        chr_safe = replace(VCF.chrom(record), r"\.|\-|\=|\/" => "_")
        chr_safer = replace(chr_safe, "|" => "_")
        pos = VCF.pos(record) |> string
        push!(locinames, chr_safer*"_"*pos)

        # get the genotype information
        geno_raw = [split(i, ('/', '|')) for i in VCF.genotype(record, :, "GT")] |> sort

        # change missing data "." to "-1"
        geno_corr_miss = map(i -> replace(i, "." => "-1"), geno_raw)

        # convert everything to an integer
        geno_int = map(i -> parse.(Int8, i), geno_corr_miss)

        # add 1 to shift genos so 0 is 1 and -1 is 0 etc.
        geno_shift = map(i -> i .+ Int8(1), geno_int)
        geno_final = [replace(i, 0 => missing) for i in geno_shift]
        geno_tuple = [Tuple(i) for i in geno_final]
        push!(locus_array, geno_tuple)
    end
    close(vcf_file)

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

    @info "\n$(abspath(infile))
$(length(sample_names)) samples detected
population info must be added <---
$(length(locinames)) loci detected"

    # create loci dataframe
    loci_df = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), locus_array)])

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = population,
        ploidy = ploidy,
        latitude = loc_xy,
        longitude = loc_xy
    )
    PopObj(samples_df, loci_df)
end


"""
    PopGen.read(infile::String; kwargs...)
Wraps `csv()`, `genepop()`, `bcf()`, and `vcf()` to read a file in as a `PopObj`. File type is
inferred from the file extension:
- delimited: .csv | .tsv | .txt
- genepop: .gen | .genepop
- variant call format: .vcf
This function uses the same keyword arguments (and defaults) as the file importing
functions it wraps; please see their respective docstrings in the Julia help console.
(e.g. `?genepop`) for specific usage details.

## Examples
```
read("cavernous_assfish.gen", markers = "msat", digits = 3)

read("bos_tauros.csv")`

read("juglans_nigra.vcf")
```
"""
function Base.read(infile::String; kwargs...)
    ext = split(infile, ".")[end]
    if ext == "gen" || ext == "genepop"
        return genepop(infile;kwargs...)

    elseif ext == "csv" || ext == "txt" || ext == "tsv"
        return delimited(infile; kwargs...)

    elseif ext == "vcf" || ext == "bcf"
        ext == "vcf" && return vcf(infile)
        ext == "bcf" && return bcf(infile)

    else
        @error "file type not recognized by filename extension \n delimited: .csv | .tsv | .txt \n genepop: .gen | .genepop \n variant call formant: .bcf | .vcf"
    end
end
