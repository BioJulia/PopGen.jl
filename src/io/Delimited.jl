#=
This file handles the import/export of delmited file formats
=#

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
