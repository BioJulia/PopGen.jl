"""
    PopObj(ind::Array{String,1}, popid::Array{Union{Int64, String},1}, loci::Array{String,1}, ploidy::Int64, genotypes::Dict, NamedTuple{(:x, :y),Tuple{Float64,Float64}})
Type "PopObj", which stores population genetics genotype data
- `ind` ::Array{String,1} of individual names
- `popid` ::Array{Union{Int64,String},1} of population names/numbers
- `loci` ::Array{String,1} of locus names in order of appearance in `genotypes`
- `ploidy` ::Int64 single integer of ploidy
- `genotypes` ::Dict of [`ind`] => genotypes::Array{String,1} ordered by `loci`
- `locdata` ::NamedTuple of decimal degrees longitude (:x) and latitude (:y)
"""
mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genotypes::Dict
    longitude::Array{Union{Int64,Float64},1}
    latitude::Array{Union{Int64,Float64},1}
    #PopObj() = new(ind, popid, loci, ploidy, genotypes, longitude, latitude)
end

function Base.show(io::IO , x::PopObj)
    println(io, "Object of type PopObj:")
    if length(x.latitude) ==0 && length(x.longitude) == 0
        println(io,"No location data provided")
    else
        println(io, "\nLongitude:")
        println(io, string.(x.longitude[1:3]) , " \u2026 ", string.(x.longitude[end-2:end]), "\n")
        println(io, "Latitude:")
        println(io, string.(x.latitude[1:3]) , " \u2026 ", string.(x.latitude[end-2:end]), "\n")
    end
    println(io, "\nNumber of individuals: $(length(x.ind))")
    println(io, x.ind[1:3] , " \u2026 ", x.ind[end-2:end], "\n")
    println(io, "Number of loci: $(length(x.loci))")
    println(io, x.loci[1:3], " \u2026 " , x.loci[end-2:end], "\n" )
    println(io, "Ploidy: $(x.ploidy)")
    println(io, "Number of populations: $(length(x.popid |> unique))","\n")
    println(io, "   #Inds | Pop","\n", "   --------------")
    popcounts = hcat([sum(x.popid .== i) for i in unique(x.popid)],unique(x.popid))
    for eachpop in 1:length(popcounts)÷2
        println(io, "\t ", popcounts[eachpop], "\t", " |", "\t", popcounts[eachpop,2])
    end
    println(io, "\nAvailable fields: ind, popid, loci, ploidy, genotypes, longitude, latitude")
end
