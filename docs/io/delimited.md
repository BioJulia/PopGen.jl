## Import a delimited file as a `PopObj`

```julia
delimited(infile; delim = ",", digits = 3, marker = "snp", diploid = true)

# Example
julia> a = delimited("/data/cali_poppy.csv", digits = 2)
```

??? warning "Windows users"
    make sure to change your backslashes "\" to forward slashes "/" 

### Arguments

- `#!julia infile::String` : path to the input file, in quotes

### Keyword Arguments

- `#!julia delim::String` : delimiter of the file, as a string. must be a single delimiter. (default: `","`)

- `#!julia digits::Int64` : the number of digits used to denote an allele (default: `3`)
- `#!julila marker::String`  : "snp" (default) or "msat" for microsatellites
- `#!julila diploid::Bool` :  uses memory-optimized parsing for diploid samples (default: `true`)



## Formatting

- Loci names must be first row
- Individuals names must be first value in row
- Population ID's must be second value in row
- longitude (x) values third value in row, latitude (y) fourth value in row
  - fill with zeroes if no location data

### example
```bash
name,population,long,lat,Locus1,Locus2,Locus3   \n
sierra_01,mountain,11.11,-22.22,001001,002002,001001   \n
sierra_02,mountain,11.12,-22.21,001001,001001,001002   \n
snbarb_03,coast,0,0,001001,001001,001002 \n
snbarb_02,coast,11.14,-22.24,001001,001001,001001 \n
snbarb_03,coast,11.15,0,001002,001001,001001 \n
```

**Fun fact**

You can also use the command `csv()` synonymously with `delimited()`. 



## Acknowledgements

Thanks to the efforts of the [CSV.jl](https://github.com/JuliaData/CSV.jl) and [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) teams, we leverage those packages to do much of the heavy lifting within this parser. 