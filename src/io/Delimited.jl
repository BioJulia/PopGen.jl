#=
This file handles the import/export of delmited file formats
=#

export delimited, csv

"""
    delimited(infile::String; delim::Union{Char,String,Regex} = "auto", digits::Int64 = 3, silent::Bool = false)
Load a delimited-type file into memory as a PopData object. *There should be no empty cells
in your file*
### Arguments
- `infile` : path to file
- `delim` : delimiter characters. By default uses auto-parsing of `CSV.File`
- `digits` : number of digits denoting each allele (default: `3`)
- `diploid`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent`   : whether to print file information during import (default: `true`)

### File formatting:
- First row is column names in this order:
    1. name
    2. population
    3. longitude
    4. latitude
    5. locus_1_name
    6. locus_2_name
    7. etc...

#### Missing data
##### Genotypes
Missing genotypes can be formatted as all-zeros `000000` or negative-nine `-9`

##### Location data
If location data is missing for a sample (which is ok!), make sure the value is written
as *0*, otherwise there will be transcription errors!

## Example
```
lizardsCA = Read.delimited("CA_lizards.csv", digits = 3);
```

### Formatting example:
```
name,population,long,lat,Locus1,Locus2,Locus3   \n
sierra_01,mountain,11.11,-22.22,001001,002002,001001   \n
sierra_02,mountain,11.12,-22.21,001001,001001,001002   \n
snbarb_03,coast,0,0,001001,001001,001002 \n
snbarb_02,coast,11.14,-22.24,001001,001001,001001 \n
snbarb_03,coast,11.15,0,001002,001001,001001 \n
```
"""
function delimited(
    infile::String;
    delim::Union{Char,String,Regex} = "auto",
    digits::Int = 3,
    diploid::Bool = true,
    silent::Bool = false
    )

    if delim == "auto"
        geno_parse = CSV.File(infile)
    else
        geno_parse = CSV.File(infile, delim = delim)
    end

    locinames = string.(keys(geno_parse[1])[5:end])

    # initiate with empty vectors for samplename, locus, and geno
    name = Vector{String}()
    loci = Vector{String}()
    popid_meta = Vector{String}()
    popid_loci = Vector{String}()
    long_array = Vector{Union{Missing,Float32}}()
    lat_array = Vector{Union{Missing,Float32}}()
    if diploid == true
        geno_raw = Vector{Int32}()
    else
        geno_raw = Vector{String}()
    end

    geno_type = determine_marker(infile, geno_parse, digits)

    if geno_type == Int16
        marker_txt = "Microsatellite"
    else
        marker_txt = "SNP"
    end

    @inbounds for samplerow in geno_parse
        vals = values(samplerow)
        append!(name, fill(vals[1], length(locinames)))
        append!(popid_loci, fill(vals[2], length(locinames)))
        append!(popid_meta, fill(vals[2],1))
        append!(long_array, vals[3])
        append!(lat_array, vals[4])
        append!(loci, locinames)
        append!(geno_raw, vals[5:end])
    end

    if !silent
        @info "\n$(abspath(infile))\n$(marker_txt) markers detected\n$(countlines(infile) - 1) samples detected\n$(length(locinames)) loci detected"
    end

    ## phase the genotypes
    if diploid == true
        genotype = map(i -> phase_dip.(i, geno_type, digits), geno_raw)
    else
        genotype = map(i -> phase.(i, geno_type, digits), geno_raw)
    end

    loci_table = table((
        name = categorical(name, true),
        population = categorical(popid_loci, true),
        locus = categorical(loci, true),
        genotype = genotype
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
        population = popid_meta,
        ploidy = ploidy,
        longitude = replace(long_array, 0.0 => missing),
        latitude = replace(lat_array, 0.0 => missing)
    ))

        return PopData(sample_table, loci_table)
end

const csv = delimited
