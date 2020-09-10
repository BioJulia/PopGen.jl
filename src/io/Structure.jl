export structure

#TODO add to docs (API)
"""
	phase_structure(datatype::DataType, args...)
Takes a DataType (such as `Int8`) and a series of integers to return
a sorted Tuple of those integers converted to that DataType. i.e. takes
a series of alleles and returns a genotype. Returns `missing` if args are
`missing`. Used internally in PopGen.structure file reader.

**Example**
```
phase_structure(Int8, 1,2,3,4,3,4,6,1)
(1, 1, 2, 3, 3, 4, 4, 6)

phase_structure(Int16, missing, missing)
missing
```
"""
function phase_structure(datatype::DataType, args...)
    all(ismissing.(args)) && return missing
    return Tuple(datatype.(sort([args...])))
end

#TODO add to docs (own page + API)
"""
    structure(infile::String; kwargs...)
Load a Structure format file into memory as a PopData object.
- `infile::String` : path to Structure file

### Keyword Arguments
- `extracols::Integer`: how many additional optional columns there are beyond Stucture's POPDATA the reader needs to ignore (default: `0`)
    - these include POPFLAG, LOCDATA, or anything else you might have added
- `extrarows::Integer` : how many additional optional rows there are beyond the first row of locus names (default: `0`)
- `missingval::String`  : the value used to identify missing values in the data (default: `"-9"`)
- `silent::Bool`   : whether to print file information during import (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci in the dataset (default: `false`)
- `faststructure::Bool`: whether the file is fastStructure format (default: `false`)

### File must follow this Structure format:
- the file is `tab` or `space` delimited **but not both**
- first row is locus names separated by the delimiter
    - leading/trailing whitespaces are tolerated
    - optional rows allowed **after** the locus names
- number of rows per sample = ploidy
    - e.g. if diploid, that sample would have 2 rows
    - multi-column variant not supported
- first data column is sample name
- second data column is population ID
    - optional columns allowed **after** the population ID (2nd) column
- remaining columns are the genotype for that individual for that locus

### Structure file example:
```
locus_1	locus_2	locus_3	locus_4	locus_5
walnut_01	1	-9	145	66	0	92
walnut_01	1	-9	-9	64	0	94
walnut_02	1	106	142	68	1	92
walnut_02	1	106	148	64	0	94
walnut_03	2	110	145	-9	0	92
walnut_03	2	110	148	66	1	-9
```

### fastStructure file format:
- the file is `tab` or `space` delimited **but not both**
- no first row of loci names
- number of rows per sample = ploidy
    - e.g. if diploid, that sample would have 2 rows
- first data column is sample name
- second data column is population ID
- remaining columns are the genotype for that individual for that locus
- usually, first 6 colums are empty (but not necessary)
- **no** extra rows or columns.
### fastStructure file example:
```
chestnut_01	1	-9	145	66	0	92
chestnut_01	1	-9	-9	64	0	94
chestnut_02	1	106	142	68	1	92
chestnut_02	1	106	148	64	0	94
chestnut_03	2	110	145	-9	0	92
chestnut_03	2	110	148	66	1	-9
```

## Example
```
walnuts = structure("juglans_nigra.str", extracols = 0, extrarows = 0)
```
"""
function structure(infile::String; silent::Bool = false, extracols::Int = 0, extrarows::Int = 0, allow_monomorphic::Bool = false, missingval::String = "-9", faststructure::Bool = false)
    # find the delimiter
    delimcheck = join(open(readlines, `head -n $(2+extrarows) $(infile)`))

    if occursin("\t", delimcheck) & occursin(" ", delimcheck)
        error("$infile contains both tab and space delimiters. Please format the file so it uses either one or the other.")
    elseif occursin("\t", delimcheck)
        delim = "\t"
    elseif occursin(" ", delimcheck)
        delim = " "
    else
        error("Please format $infile to be either tab or space delimited")
    end
    
    if faststructure == true
        data_row = 1
        first_row = split(strip(open(readline, infile)), delim)
        n_loci = length(first_row)
        locinames = ["marker_$i" for i in 1:n_loci-2]
    else
        data_row = 2 + extrarows
        locinames = Symbol.(split(strip(open(readline, infile)), delim))
    end
    # read in the file as a table
    geno_parse = CSV.File(
        infile,
        delim = delim,
        header = false,
        datarow = data_row,
        missingstrings = [missingval],
        #type = type,
        ignorerepeated = true
        ) |> DataFrame
        # ignore any extra columns
        if !iszero(extracols)
            geno_parse = geno_parse[!, Not(collect(3:2+n))]
        end
        
    # determine marker from max genotype value
    maxval = map(i -> maximum(skipmissing(i)), eachcol(geno_parse[3:end])) |> maximum
    markertype = maxval < 100 ? Int8 : Int16
    if faststructure == true
        if markertype == Int8
            locinames = Symbol.(replace.(locinames, "marker" => "snp"))
        else
            locinames = Symbol.(replace.(locinames, "marker" => "msat"))
        end
    end
    
    rename!(geno_parse, append!([:name, :population], locinames))
    
    # fix names, just in case
    geno_parse.name .= replace.(geno_parse.name, "-" => "_")
    geno_parse.population = string.(geno_parse.population)
    by_sample = groupby(geno_parse, :name)
    
    
    # create new dataframe to add phased genotypes to
    loci_df = DataFrame(:locus => locinames)
    for (key, eachsample) in pairs(by_sample)
        insertcols!(
            loci_df,    
            Symbol(key.name) => map(i -> phase_structure(markertype, i...), eachcol(eachsample[3:end]))
        )
    end

    # hacky way of tranposing the dataframe
    loci_df = loci_df[!, Not(:locus)] |> Matrix |> permutedims |> DataFrame
    rename!(loci_df, locinames)

    # create the metadata from the original file info
    meta_df = DataFrames.combine(
        by_sample,
        :name => first => :name,
        :population => first => :population,
        :name => (i -> Int8(length(i))) => :ploidy  # ploidy derived from # of times a name appears
    )
    insertcols!(
        meta_df, 
        :longitude => Vector{Union{Missing, Float32}}(undef, length(meta_df.name)), 
        :latitude => Vector{Union{Missing, Float32}}(undef, length(meta_df.name))
    )
    
    if !silent
        @info "\n$(abspath(infile))\n$(length(meta_df.name)) samples from $(length(unique(meta_df.population))) populations detected\n$(length(locinames)) loci detected"
    end

    # create long-format DataFrame
    insertcols!(loci_df, 1, :name => meta_df.name, :population => meta_df.population)
    loci_df = DataFrames.stack(loci_df, DataFrames.Not([:name, :population]))
    rename!(loci_df, [:name, :population, :locus, :genotype])

    # convert columns to PooledArrays
    loci_df.name = PooledArray(loci_df.name)
    loci_df.population = PooledArray(loci_df.population)
    loci_df.locus = PooledArray(loci_df.locus)

    if allow_monomorphic 
        pd_out = PopData(meta_df, loci_df)
    else
        pd_out = drop_monomorphic!(PopData(meta_df, loci_df))
    end
    return pd_out
end

