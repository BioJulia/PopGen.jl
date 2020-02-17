using BenchmarkTools, JuliaDBMeta, CategoricalArrays, StatsBase, JuliaDB, PopGen

x = gulfsharks();

abstract type PopObj end

mutable struct PopSample <: PopObj
    name::String
    population::String
    ploidy::Int8
    longitude::Union{Missing, Float32}
    latitude::Union{Missing, Float32}
end

struct PopData <: PopObj
    samples::Vector{PopSample}
    loci::T where T<:IndexedTable
end

a = Vector{PopSample}()

for i in eachrow(x.samples)
    ps = PopSample(i.name, i.population, i.ploidy, i.longitude, i.latitude)
    push!(a, ps)
end


function Base.show(io::IO, ::MIME"text/plain", popsample::PopSample3)
    println(io, "PopSample")
    printstyled(io, "  name: ", bold = true)
    println(io, popsample.name)

    printstyled(io, "  population: ", bold = true)
    println(io, popsample.population)

    printstyled(io, "  ploidy: ", bold = true)
    println(io, popsample.ploidy)

    printstyled(io, "  longitude: ", bold = true)
    println(io, popsample.longitude)

    printstyled(io, "  latitude: ", bold = true)
    println(io, popsample.latitude)
end

function Base.show(io::IO, ::MIME"text/plain", z::Vector{PopSample})
    print(io, "Vector{PopSample} with $(length(z)) sample records\n")
    show(sampletable(z))
end

function samplenames(data::Vector{PopSample})
    map(i -> i.name, data)
end

function samplepop(data::Vector{PopSample})
    map(i -> i.population, data)
end

function sampleploidy(data::Vector{PopSample})
    map(i -> i.ploidy, data)
end

function samplelong(data::Vector{PopSample})
    map(i -> i.longitude, data)
end

function samplelat(data::Vector{PopSample})
    map(i -> i.latitude, data)
end

function sampleloc(data::Vector{PopSample})
    hcat(samplelong(data), samplelat(data))
end

function sampletable(data::Vector{PopSample})
    table(
        samplenames(data),
        samplepop(data),
        sampleploidy(data),
        samplelat(data),
        samplelong(data),
        names = [:name, :population, :ploidy, :longitude, :latitude],
        )
end

const PopSample3 = NamedTuple{(:name, :population, :ploidy, :longitude, :latitude),Tuple{String,String,Int8,Union{Missing,Float32},Union{Missing,Float32}}}

function Base.show(io, sample::NamedTuple{(:name, :population, :ploidy, :longitude, :latitude),Tuple{String,String,Int8,Float32,Float32}})
    x = PopSample(sample...)
    show(x)
end
