"""
    Read.genepop(infile::String; ploidy::Int64 = 2, popsep::Any = "POP", numpop::Int64)
Load a Genepop format file into memory as a PopObj object.
- `infile` : path to Genepop file
- `ploidy` : ploidy of the organism
- `popsep` : word that separates populations in `infile` (default: "POP")
- `numpop` : number of populations in `infile` (used for checking parser)
File must follow standard Genepop formatting:
- First line is a comment (and skipped)
- Loci are listed after comment as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default "POP") must delimit populations
- File is tab or space delimted

usage: `waspsNY = Read.genepop("wasp_hive.gen", ploidy = 2, popsep = "POP", numpops = 2);`

Genepop file example:  \n
---------------------
Wasp populations in New York \n
Locus1 \n
Locus2 \n
Locus3 \n
POP \n
Oneida_01,  250230 564568 110100  \n
Oneida_02,  252238 568558 100120  \n
Oneida_03,  254230 564558 090100  \n
POP \n
Newcomb_01,  254230 564558 080100 \n
Newcomb_02,  000230 564558 090080 \n
Newcomb_03,  254230 000000 090100 \n
Newcomb_04,  254230 564000 090120 \n
---------------------
"""
function genepop(infile::String; ploidy::Int64 = 2, popsep::Any = "POP", numpops::Int64)
    gpop = split(open(readlines,infile)[2:end], popsep)
    if length(gpop)-1 != numpops
        error("incorrect number of populations detected, see docstring for formatting
            expected : $numpops
            detected : $(length(gpop)-1) ")
    end
    if length(split(gpop[1], ",")) > 1
        locinames = split(gpop[1],",") |> Array{String,1}
    else
        locinames = gpop[1]
    end
    d = Dict(i => [] for i in locinames)
    popid = []
    indnames = []
    for i in 2:length(gpop)
        append!(popid, fill(i-1,length(gpop[i])))
        for j in 1:length(gpop[i])
            phasedloci = []
            push!(indnames, split( strip(gpop[i][j]), r"\,|\t")[1])
            unphasedloci = split( strip(gpop[i][j]), r"\s|\t")[2:end] |> Array{String,1}
            # phase genotypes by ploidy
            for locus in unphasedloci
                phasedlocus = parse.(Int64,[join(i) for i in Iterators.partition(locus,length(locus)÷ploidy)])  |> Tuple
                push!(phasedloci, phasedlocus)
            end
            for (loc,geno) in zip(locinames, phasedloci)
                push!(d[loc], geno)
            end
        end
    end
    PopObj(indnames,popid,locinames,ploidy,d,[],[])
end


"""
    csv(infile::String; delim::Union{Char,String,Regex}, ploidy::Int64 = 2, location::Bool = false)
Load a CSV-type file into memory as a PopObj object
- `infile` : path to Genepop file
- `delim` values can be space (" "), comma (","), tab ("\\t"), etc.
- `ploidy` : ploidy of the organism
- `location` : decimal degrees longitude/latitude provided as values 3/4
File formatting:
- Loci names must be first row
- Individuals names must be first value in row
- Population ID's must be second value in row
- [Optional] longitude (x) values third value in row, latitude (y) fourth

example: `lizardsCA = Read.csv("CA_lizards.csv", delim = ",", ploidy = 2);`

Formatting example:  \n
---------------------  \n
Locus1,Locus2,Locus3   \n
sierra_01,1,001001,002002,001001   \n
sierra_02,1,001001,001001,001002   \n
snbarb_03,2,001001,001001,001002 \n
snbarb_02,2,001001,001001,001001 \n
snbarb_03,2,001002,001001,001001 \n
---------------------
"""
function csv(infile::String; delim::Union{Char,String,Regex} = ",", ploidy::Int64 = 2, location::Bool = false)
    gpop = open(readlines,infile)
    locinames = split(gpop[1], delim)
    d = Dict(i => [] for i in locinames)
    popid = []
    indnames = []
    locx = []
    locy = []
    if location == false
        for i in gpop[2:end]
            tmp = split(i, delim) |> Array{String,1}
            phasedloci = []
            # phase genotypes by ploidy
            for locus in tmp[3:end]
                phasedlocus = parse.(Int64,[join(i) for i in Iterators.partition(locus,length(locus)÷ploidy)])|> Tuple
                push!(phasedloci, phasedlocus)
            end
            for (loc,geno) in zip(locinames, phasedloci)
                push!(d[loc], geno)
            end
            push!(indnames, tmp[1])
            push!(popid, tmp[2])
        end
    else
        for i in gpop[2:end]
            tmp = split(i, delim) |> Array{String,1}
            # phase loci according to ploidy
            phasedloci = []
            for locus in tmp[5:end]
                phasedlocus = parse.(Int64,[join(i) for i in Iterators.partition(locus,length(locus)÷ploidy)])|> Tuple
                push!(phasedloci, phasedlocus)
            end
            for (loc,geno) in zip(locinames, phasedloci)
                push!(d[loc], geno)
            end
            push!(indnames, tmp[1])
            push!(popid, tmp[2])
            push!(locx, tmp[3])
            push!(locy, tmp[4])
        end
    end
    println("\n", "Input File : ", abspath(infile) )
    PopObj(indnames, popid, locinames, ploidy, d, locx, locy)
end
