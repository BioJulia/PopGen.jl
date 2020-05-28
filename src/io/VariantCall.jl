#=
This file handles the import/export of Variant Call Format files
=#

export bcf, vcf

infile = "/home/pdimens/PopGen.jl/data/source/filtered_oyster.vcf"

"""
    bcf(infile::String)
Load a BCF file into memory as a PopObj object. Population and [optional]
location information need to be provided separately.
- `infile` : path to BCF file
"""
function bcf(infile::String)
    vcf_file = BCF.Reader(open(infile, "r"))

    # get sample names from header
    sample_names = header(vcf_file).sampleID

    # fill in pop/lat/long with missing
    population = fill("missing", length(sample_names))
    loc_xy = Vector{Union{Missing,Float32}}(undef, length(sample_names))

    # get loci names
    locinames = Vector{String}()

    ## array of genotypes
    locus_array = Vector{PoolGenotypeArray}()

    # get genotypes
    for record in vcf_file
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
        push!(locus_array, Tuple.(geno_final))
    end
    close(vcf_file)
    @info "\n$(abspath(infile))\n$(length(sample_names)) samples detected\npopulation info must be added <---\n$(length(locinames)) loci detected"

    # create loci dataframe
    loci_df = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), locus_array)])
    populations = fill("missing", length(sample_names))
    insertcols!(loci_df, 1, :name => sample_names, :population => populations)
    geno_parse = DataFrames.stack(loci_df, DataFrames.Not(1:2))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    categorical!(geno_parse, [:name, :population, :locus], compress = true)

    ploidy = DataFrames.combine(
        groupby(dropmissing(geno_parse), :name),
        :genotype => find_ploidy => :ploidy
    ).ploidy

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = populations,
        ploidy = ploidy,
        latitude = loc_xy,
        longitude = loc_xy
    )
    PopData(samples_df, geno_parse)
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
    loc_xy = Vector{Union{Missing,Float32}}(undef, length(sample_names))

    # get loci names
    locinames = Vector{String}()

    ## array of genotypes
    locus_array = Vector{PoolGenotypeArray}()

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
        push!(locus_array, Tuple.(geno_final))
    end
    close(vcf_file)
    @info "\n$(abspath(infile))\n$(length(sample_names)) samples detected\npopulation info must be added <---\n$(length(locinames)) loci detected"

    # create loci dataframe
    loci_df = DataFrame([j => k for (j,k) in zip(Symbol.(locinames), locus_array)])
    populations = fill("missing", length(sample_names))
    insertcols!(loci_df, 1, :name => sample_names, :population => populations)
    geno_parse = DataFrames.stack(loci_df, DataFrames.Not(1:2))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    categorical!(geno_parse, [:name, :population, :locus], compress = true)

    ploidy = DataFrames.combine(
        groupby(dropmissing(geno_parse), :name),
        :genotype => find_ploidy => :ploidy
    ).ploidy

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = populations,
        ploidy = ploidy,
        latitude = loc_xy,
        longitude = loc_xy
    )
    PopData(samples_df, geno_parse)
end
