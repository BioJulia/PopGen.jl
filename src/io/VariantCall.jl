export bcf, vcf

"""
    bcf(infile::String; ; rename_snp::Bool, silent::Bool, allow_monomorphic::Bool)
Load a BCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to BCF file (can be gzipped)
- `rename_snp` : true/false of whether to simplify marker names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:

|Base| | Allele |
|:---:|:---:|:---:|
| A | => | 1 |
| T | => | 2 |
| C | => | 3 |
| G | => | 4 |
"""
function bcf(infile::String; rename_snp::Bool = false, silent::Bool = false, allow_monomorphic::Bool = false)
    bases = (A = Int8(1), T = Int8(2), C = Int8(3), G = Int8(4), miss = Int8(0))
    stream = BCF.Reader(openbcf(infile))
    nmarkers = countlines(openbcf(infile)) - length(BCF.header(stream)) - 1
    sample_ID = header(stream).sampleID
    nsamples = length(sample_ID)
    loci_names = fill("marker", nmarkers)
    geno_df = DataFrame(fill((Int8(0), Int8(0)), nsamples, nmarkers))
    @info "\n$(abspath(infile))\n$nsamples samples detected\n$nmarkers markers detected\npopulation info must be added <---"
    for (idx,record) in enumerate(stream)
        ref_alt = Dict(-1 => "miss", 0 => VCF.ref(record), [i => j for (i,j) in enumerate(VCF.alt(record))]...)
        raw_geno = VCF.genotype(record, 1:nsamples, "GT")
        conv_geno = map(raw_geno) do rg
            tmp = replace.(rg, "." => "-1")
            ig = collect(parse.(Int8, split(tmp, r"\/|\|")))
            [bases[Symbol(ref_alt[i])] for i in ig] |> sort |> Tuple
        end
        geno_df[:, idx] = conv_geno
        loci_names[idx] = VCF.chrom(record) * "_" * string(VCF.pos(record))
    end
    close(stream)
    if rename_snp
        rename!(geno_df, [Symbol.("snp_" * i) for i in string.(1:nmarkers)])
    else
        rename!(geno_df, Symbol.(loci_names))
    end
    insertcols!(geno_df, 1, :name => sample_ID, :population => "missing")
    stacked_geno_df = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(stacked_geno_df, [:name, :population, :locus, :genotype])
    select!(stacked_geno_df, :name, :population, :locus, :genotype => (i -> replace(i, (0,0) => missing)) => :genotype)
    sort!(stacked_geno_df, [:name, :locus])
    # ploidy finding
    gdf = DataFrames.groupby(stacked_geno_df, :name)
    meta_df = combine(gdf,
        :genotype => (i -> Int8(mode(length.(skipmissing(i) |> collect)))) => :ploidy    
    )
    insertcols!(meta_df, 2, :population => "missing")
    insertcols!(meta_df, 4, :longitude => Vector{Union{Missing, Float32}}(undef, nsamples), :latitude => Vector{Union{Missing, Float32}}(undef, nsamples))
    if allow_monomorphic 
        pd_out = PopData(meta_df, stacked_geno_df)
    else
        pd_out = drop_monomorphic!(PopData(meta_df, stacked_geno_df))
    end
    return pd_out
end

### VCF parsing ###

"""
    vcf(infile::String; ; rename_snp::Bool, silent::Bool, allow_monomorphic::Bool)
Load a VCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to VCF file (can be gzipped)
- `rename_snp` : true/false of whether to simplify marker names to `snp_#` (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`)
Alleles are recoded according to the following schema:

|Base| | Allele |
|:---:|:---:|:---:|
| A | => | 1 |
| T | => | 2 |
| C | => | 3 |
| G | => | 4 |
"""
function vcf(infile::String; rename_snp::Bool = false, silent::Bool = false, allow_monomorphic::Bool = false)
    bases = (A = Int8(1), T = Int8(2), C = Int8(3), G = Int8(4), miss = Int8(0))
    #TODO @require for GeneticVariation and GZip
    stream = VCF.Reader(openvcf(infile))
    nmarkers = countlines(openvcf(infile)) - length(VCF.header(stream)) - 1
    sample_ID = header(stream).sampleID
    nsamples = length(sample_ID)
    loci_names = fill("marker", nmarkers)
    geno_df = DataFrame(fill((Int8(0), Int8(0)), nsamples, nmarkers))
    @info "\n$(abspath(infile))\n$nsamples samples detected\n$nmarkers markers detected\npopulation info must be added <---"
    for (idx,record) in enumerate(stream)
        ref_alt = Dict(-1 => "miss", 0 => VCF.ref(record), [i => j for (i,j) in enumerate(VCF.alt(record))]...)
        raw_geno = VCF.genotype(record, 1:nsamples, "GT")
        conv_geno = map(raw_geno) do rg
            tmp = replace.(rg, "." => "-1")
            ig = collect(parse.(Int8, split(tmp, r"\/|\|")))
            [bases[Symbol(ref_alt[i])] for i in ig] |> sort |> Tuple
        end
        geno_df[:, idx] = conv_geno
        loci_names[idx] = VCF.chrom(record) * "_" * string(VCF.pos(record))
    end
    close(stream)
    if rename_snp
        rename!(geno_df, [Symbol.("snp_" * i) for i in string.(1:nmarkers)])
    else
        rename!(geno_df, Symbol.(loci_names))
    end
    insertcols!(geno_df, 1, :name => sample_ID, :population => "missing")
    stacked_geno_df = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(stacked_geno_df, [:name, :population, :locus, :genotype])
    # replace (0,0) as missing
    select!(stacked_geno_df, :name, :population, :locus, :genotype => (i -> replace(i, (0,0) => missing)) => :genotype)
    sort!(stacked_geno_df, [:name, :locus])
    # ploidy finding
    #TODO different ploidy finding method
    gdf = DataFrames.groupby(stacked_geno_df, :name)
    meta_df = combine(gdf,
        :genotype => (i -> Int8(mode(length.(skipmissing(i) |> collect)))) => :ploidy    
    )
    insertcols!(meta_df, 2, :population => "missing")
    insertcols!(meta_df, 4, :longitude => Vector{Union{Missing, Float32}}(undef, nsamples), :latitude => Vector{Union{Missing, Float32}}(undef, nsamples))
    if allow_monomorphic 
        pd_out = PopData(meta_df, stacked_geno_df)
    else
        pd_out = drop_monomorphic!(PopData(meta_df, stacked_geno_df))
    end
    return pd_out
end



"""
    openvcf(::String)
Open VCF file (`.vcf` or `.vcf.gz`) and return an `IO` stream in reading mode `"r"`.
Adapted from OpenMendel/VCFTools.jl
https://github.com/OpenMendel/VCFTools.jl/blob/master/src/gtstats.jl#L169
"""
function openvcf(infile::String)
    if endswith(infile, ".vcf")
        return open(infile, "r")
    elseif endswith(infile, ".vcf.gz")
        return GZip.open(infile, "r")
    else
        throw(ArgumentError("The filename must end with .vcf or .vcf.gz"))
    end
end

"""
    openbcf(::String)
Open BCF file (`.bcf` or `.bcf.gz`) and return an `IO` stream in reading mode `"r"`.
Adapted from OpenMendel/VCFTools.jl
https://github.com/OpenMendel/VCFTools.jl/blob/master/src/gtstats.jl#L169
"""
function openbcf(infile::String)
    if endswith(infile, ".bcf")
        return open(infile, "r")
    elseif endswith(infile, ".bcf.gz")
        return GZip.open(infile, "r")
    else
        throw(ArgumentError("The filename must end with .bcf or .bcf.gz"))
    end
end