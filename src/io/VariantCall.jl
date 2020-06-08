#=
This file handles the import/export of Variant Call Format files
=#

export bcf, vcf

infile = "/home/pdimens/PopGen.jl/data/source/filtered_oyster.vcf"

"""
    bcf(infile::String; silent::Bool = false)
Load a BCF file into memory as a PopObj object. Population and [optional]
location information need to be provided separately. Use `silent=true` to supress
printing during file loading.
- `infile` : path to VCF file
"""
function bcf(infile::String; silent::Bool = false)
    vcf_file = BCF.Reader(open(infile, "r"))

    # get sample names from header
    sample_names = header(vcf_file).sampleID

    # fill in pop/lat/long with missing
    population = fill("missing", length(sample_names))
    loc_xy = Vector{Union{Missing,Float32}}(undef, length(sample_names))

    # get loci names
    locinames = Vector{String}()

    ## genotype dataframe
    geno_df = DataFrame()

    # get genotypes
    for record in vcf_file
        # fix locus names to be syntax-safe
        chr_safe = replace(BCF.chrom(record), r"\.|\-|\=|\/" => "_")
        chr_safer = replace(chr_safe, "|" => "_")
        pos = VCF.pos(record) |> string
        push!(locinames, chr_safer*"_"*pos)

        # get the genotype information
        geno_raw = [split(i, ('/', '|')) for i in BCF.genotype(record, :, "GT")] |> sort

        # change missing data "." to "-1"
        geno_corr_miss = map(i -> replace(i, "." => "-1"), geno_raw)

        # convert everything to an integer
        geno_int = map(i -> parse.(Int8, i), geno_corr_miss)

        # add 1 to shift genos so 0 is 1 and -1 is 0 etc.
        geno_shift = map(i -> i .+ Int8(1), geno_int)
        geno_out = Vector{Union{Missing, Genotype}}()
        @inbounds for i in geno_shift
            if all(iszero.(i)) == true
                push!(geno_out, missing)
            else
                push!(geno_out, Tuple(i))
            end
        end
        insertcols!(geno_df, Symbol(chr_safer*"_"*pos) => geno_out)
    end
    close(vcf_file)
    if !silent
        @info "\n$(abspath(infile))\n$(length(sample_names)) samples detected\npopulation info must be added <---\n$(length(locinames)) loci detected"
    end

    # create loci dataframe
    insertcols!(geno_df, 1, :name => sample_names, :population => population)
    geno_parse = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    categorical!(geno_parse, [:name, :population, :locus], compress = true)

    # make sure levels are sorted by order of appearance
    levels!(geno_parse.locus, unique(geno_parse.locus))
    levels!(geno_parse.name, unique(geno_parse.name))
    ploidy = DataFrames.combine(
        groupby(geno_parse, :name),
        :genotype => (i -> find_ploidy(i[i .!== missing])) => :ploidy
    ).ploidy

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = population,
        ploidy = ploidy,
        latitude = loc_xy,
        longitude = loc_xy
    )
    PopData(samples_df, geno_parse)
end

### VCF parsing ###

"""
    vcf(infile::String; silent::Bool = false)
Load a VCF file into memory as a PopObj object. Population and [optional]
location information need to be provided separately. Use `silent=true` to supress
printing during file loading.
- `infile` : path to VCF file
"""
function vcf(infile::String, silent::Bool = false)
    vcf_file = VCF.Reader(open(infile, "r"))

    # get sample names from header
    sample_names = header(vcf_file).sampleID

    # fill in pop/lat/long with missing
    population = fill("missing", length(sample_names))
    loc_xy = Vector{Union{Missing,Float32}}(undef, length(sample_names))

    # get loci names
    locinames = Vector{String}()

    ## genotype dataframe
    geno_df = DataFrame()

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
        geno_out = Vector{Union{Missing, Genotype}}()
        @inbounds for i in geno_shift
            if all(iszero.(i)) == true
                push!(geno_out, missing)
            else
                push!(geno_out, Tuple(i))
            end
        end
        insertcols!(geno_df, Symbol(chr_safer*"_"*pos) => geno_out)
    end
    close(vcf_file)
    if !silent
        @info "\n$(abspath(infile))\n$(length(sample_names)) samples detected\npopulation info must be added <---\n$(length(locinames)) loci detected"
    end

    # create loci dataframe
    insertcols!(geno_df, 1, :name => sample_names, :population => population)
    geno_parse = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    categorical!(geno_parse, [:name, :population, :locus], compress = true)

    # make sure levels are sorted by order of appearance
    levels!(geno_parse.locus, unique(geno_parse.locus))
    levels!(geno_parse.name, unique(geno_parse.name))
    ploidy = DataFrames.combine(
        groupby(geno_parse, :name),
        :genotype => (i -> find_ploidy(i[i .!== missing])) => :ploidy
    ).ploidy

    # create samples df
    samples_df = DataFrame(
        name = sample_names,
        population = population,
        ploidy = ploidy,
        latitude = loc_xy,
        longitude = loc_xy
    )
    PopData(samples_df, geno_parse)
end
