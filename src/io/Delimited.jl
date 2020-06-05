#=
This file handles the import/export of delmited file formats
=#

export delimited, csv
#TODO update API docs
"""
    delimited(infile::String; delim::Union{Char,String,Regex} = "auto", digits::Int64 = 3, silent::Bool = false)
Load a delimited-type file into memory as a PopData object. *There should be no empty cells
in your file*
### Arguments
- `infile` : path to file

### Keyword Arguments
- `delim` : delimiter characters. By default uses auto-parsing of `CSV.File`
- `digits` : number of digits denoting each allele (default: `3`)
- `diploid`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent`   : whether to print file information during import (default: `true`)

### File formatting:
- The first row is column names (names don't matter)
- The columns must be in this order:
    1. sample name
    2. population name
    3. longitude
    4. latitude
    5. locus_1 genotypes
    6. locus_2 genotypes
    7. etc...

#### Missing data
##### Genotypes
Missing genotypes can be formatted as all-zeros `000000`, left empty, or negative-nine `-9`

##### Location data
If location data is missing for a sample (which is ok!), make sure the value is
blank, otherwise there will be transcription errors! (look at line 3 in the example below)

## Example
```
lizardsCA = Read.delimited("CA_lizards.csv", digits = 3);
```

### Formatting example:
```
name,population,long,lat,Locus1,Locus2,Locus3   \n
sierra_01,mountain,11.11,-22.22,001001,-9,001001   \n
sierra_02,mountain,11.12,-22.21,001001,001001,001002   \n
snbarb_01,coast,,,001001,000000,001002 \n
snbarb_02,coast,11.14,-22.24,001001,001001,001001 \n
snbarb_03,coast,11.15,,001002,001001,001001 \n
```
"""
function delimited(
    infile::String;
    delim::Union{Char,String,Regex} = "auto",
    digits::Int = 3,
    diploid::Bool = true,
    silent::Bool = false
    )

    diploid ? type = nothing : type = String
    dlm = delim == "auto" ? nothing : delim

    file_parse = CSV.read(infile, delim = dlm, missingstrings = ["-9", ""], type = type) |> DataFrame

    locinames = names(file_parse)[5:end]

    meta = select(file_parse, 1:4)
    rename!(meta, [:name, :population, :longitude, :latitude])
    select!(file_parse, Not(3:4))
    geno_type = determine_marker(file_parse, digits)
    geno_parse = DataFrames.stack(file_parse, DataFrames.Not(1:2))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    categorical!(geno_parse, [:name, :population, :locus], compress = true)
    transform!(geno_parse, :genotype => (i -> phase.(i, geno_type, digits)) => :genotype)

    if !silent
        @info "\n$(abspath(infile))\n$(length(meta[!, 1])) samples across $(length(unique(meta[!,2]))) populations detected\n$(length(locinames)) loci detected"
    end

    # make sure levels are sorted by order of appearance
    levels!(geno_parse.locus, unique(geno_parse.locus))
    levels!(geno_parse.name, unique(geno_parse.name))

    ploidy = DataFrames.combine(
        groupby(dropmissing(geno_parse), :name),
        :genotype => find_ploidy => :ploidy
    ).ploidy
    DataFrames.insertcols!(meta, 3, :ploidy => ploidy)
    PopData(meta, geno_parse)
end

const csv = delimited

"""
    popdata2delimited(data::PopData; filename::String, delim::String = ",", digits::Integer = 3, format::String = "wide")
Write PopData to a text-delimited file. 
### Keyword Arguments
- `filename`: a `String` of the output filename
- `digits` : an `Integer` indicating how many digits to format each allele as (e.g. `(1, 2)` => `001002` for `digits = 3`)
- `format` : a `String` indicating whether to output in`"wide"` or `"long"` (aka `"tidy"`) format 
- `delim` : the `String` delimiter to use for writing the file. 

### Example
```julia
cats = nancycats();
fewer_cats = omit_samples(cats, samples(cats)[1:10]);
popdata2delimited(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "wide", delim = " ")
```
"""
function popdata2delimited(data::PopData; filename::String, delim::String = ",", digits::Integer = 3, format::String = "wide")
    unphased_df = select(data.loci, 1:3, 4 => (i -> unphase.(i, digits = digits)) => :genotype)
    if format == "wide"
        wide_df = DataFrames.unstack(unphased_df, :locus, :genotype)
        insertcols!(wide_df, 3, :longitude => data.meta.longitude, :latitude => data.meta.latitude)
        CSV.write(filename, wide_df, delim = delim) ;
        return
    else
        out_df = sort(
                    transform(
                        unphased_df, 1:2, 
                        3 => (i -> Vector{Union{Float32, Missing}}(undef, length(i))) => :longitude, 
                        4 => (i -> Vector{Union{Float32, Missing}}(undef, length(i))) => :latitude, 3:4
                        )
                    )
        for i in 1:length(data.meta.name)
            out_df[out_df.name .== data.meta.name[i], :longitude] .= data.meta.longitude[i]
            out_df[out_df.name .== data.meta.name[i], :latitude] .= data.meta.latitude[i]
        end
        CSV.write(filename, out_df, delim = delim) ;
        return
    end
end





for (i,j) in zip(samples, populations)
    data.meta[data.meta.name .== i, :population] .= j
    data.loci[data.loci.name .== i, :population] .= j
end