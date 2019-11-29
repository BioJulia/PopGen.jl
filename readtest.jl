struct PopSample
    name::String
    population::Union{String, Integer}
    ploidy::Integer
    longitude::Union{Missing,Float64}
    latitude::Union{Missing,Float64}
    genotypes::Vector{Union{Missing,Tuple}}
end


function phase_loci(indiv::String, locnames::Vector{String} ndigits::Integer, gtype::DataType)
   phasedloci = []
   row = split(strip(indiv), r"\,|\s|\t") |>Array{String,1}
   indname = row[1]
   unphasedloci = row[3:end]
   replace!(unphasedloci, "-9" => "0") #just in case -9 = missing
      for locus in unphasedloci
        phasedlocus = parse.(gtype,
            [join(i) for i in Iterators.partition(locus, ndigits)],) |> sort |> Tuple
         miss_geno = fill(0, length(phasedlocus)) |> Tuple
         if phasedlocus == miss_geno
            phasedlocus = missing
         end
         replace!(phasedloci, "0" => missing)
         push!(phasedloci, phasedlocus)
      end
      ploidy = mode(length.(phasedloci[1:end÷5] |> skipmissing))
      return (indname, ploidy, Array{Union{Missing,Tuple}}(phasedloci))
end



    geno_type = Int16
    ndigits = 3
    popid = []
    indnames = []
    loci = []
#    gpop = split(open(readlines, infile)[2:end], popsep)
#    gpop = split(open(readlines, "C:/Users/pdime/Desktop/PopGen.jl/data/data/gulfsharks.gen")[2:end], "POP")
    gpop = split(open(readlines, "/home/pdimens/PopGen.jl/data/data/gulfsharks.gen")[2:end], "POP")
#    if length(gpop) - 1 != numpops
#        error("incorrect number of populations detected, see docstring for formatting
#            expected : $numpops
#            detected : $(length(gpop)-1) ")
#    end
    if length(gpop[1]) == 1     # loci horizontally stacked
        locinames = strip.(split(gpop[1] |> join, ",") |> Array{String,1})
        replace!(locinames, "." => "_")
        locinames = Symbol.(locinames)
    else                        # loci vertically stacked
        locinames = Symbol.(replace(gpop[1], "." => "_"))
    end

    for i = 2:length(gpop)
        #append!(popid, fill(i - 1, length(gpop[i])))
        for j in gpop[i]
            ind, ploidy, loc = phase_loci(j, ndigits,geno_type)
            popsample = PopSample(ind, i-1, ploidy,missing, missing,loc)
            #push!(indnames, ind)
            #push!(loci, loc)
        end
    end
#TODO
    # infer ploidy by taking the mode of # alleles per locus per indiv for
    # the first 20% of loci


    out_tbl = [(name = indnames,
               population = popid,
               ploidy = ploidy,
               longitude = fill(missing, length(indnames)),
               latitude = fill(missing, length(indnames))
               )]

    out_out = [(name = i, population = j, ploidy = k) for (i,j,k) in zip(indnames, popid, ploidy)]

genotest = Dict(i => [] for i in names(x.loci))

for (i,j) in zip(names(x.loci),eachcol(x.loci))
    genotest[i] = j
end

typedgeno = Dict()
for (k,v) in genotest
    typedgeno[k] = Vector{Union{Missing, Tuple}}(v)
end

tb = (; typedgeno...) |> table ;
tb





### MAKE TABLE AND THEN MAKE IT A TABLE.TABLE
tbl = [(a = 3, b = 2)]
push!(tbl, (a = 2, b = 3))



map(x -> mode(length.(x[1:end÷2] |> skipmissing)), loci)
