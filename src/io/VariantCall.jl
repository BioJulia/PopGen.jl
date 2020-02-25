#=
This file handles the import/export of Variant Call Format files
=#

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
