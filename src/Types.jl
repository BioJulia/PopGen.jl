export PopObj, PopData, show, Genotype, GenoArray

"""
    AbstractType PopObj
Generic AbstractType for use in PopGen.jl
"""
abstract type PopObj end

"""
```
PopData
    meta::DataFrame
    loci::DataFrame
```
The data struct used for the PopGen population genetics ecosystem. You are
STRONGLY discouraged from manually creating tables to pass into a PopObj,
and instead should use the provided file importers.

- `meta` ::DataFrame individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names
    - `ploidy` ::Int8 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `loci` ::DataFrame Long-format table of sample genotype records
    - name` ::CategoricalString the individual/sample names
    - `population`::CategoricalString population name
    - `locus` ::CategoricalString of locus name
    - `genotype` Tuple of Int8 or Int16 depending on SNP or microsatellite
"""
struct PopData <: PopObj
    meta::DataFrame
    loci::DataFrame
end

#=
#TODO needs to be included in the struct definition
function PopData(meta::DataFrame, loci::DataFrame)
    tmp = issorted(loci, [:locus, :name], lt = natural) ? loci : sort(loci, [:locus, :name], lt = natural)
    sort(meta.name) != sort(tmp.name.pool) && throw(ArgumentError("meta and loci dataframes do not contain the same sample names"))
    sort(unique(meta.population)) != sort(tmp.population.pool) && throw(ArgumentError("meta and loci dataframes do not contain the same population names"))
    PopData(meta, tmp)
end
=#

"""
    Genotype::DataType
For convenience purposes, an alias for `NTuple{N, <:Integer} where N`, which is
the type describing individual genotypes in PopData.
"""
const Genotype = NTuple{N, <:Integer} where N


"""
    GenoArray::DataType
For convenience purposes, an alias for an `AbstractVector` of elements `Missing`
and `Genotype`, which itself is of type `NTuple{N, <:Integer} where N`.
The definition as an `AbstractVector` adds flexibility for `SubArray`
cases.
"""
const GenoArray = AbstractVector{S} where S<:Union{Missing,Genotype}


"""
    Base.show(io::IO, data::PopData)
Overloads `Base.show` for concise summary printing of a PopData object.
"""
function Base.show(io::IO, data::PopData)
    println(io, "PopData Object")
    if occursin("Int16", string(eltype(data.loci.genotype)))
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    n_loc = length(unique(data.loci.locus))
    printstyled(io, n_loc, " ", marker, " markers", "\n" , bold = true)
    if "ploidy" ∈ names(data.meta)
        ploidy = unique(data.meta.ploidy) |> sort
        if length(ploidy) == 1
            print(io, "  Ploidy: ") ; printstyled(io, ploidy |> join, "\n", bold = true)
        else
            print(io, "  Ploidy (varies): ")
            printstyled(io, ploidy[1], bold = true); [printstyled(io, ", $i", bold = true) for i in ploidy[2:end]]
            print(io, "\n")
        end
    else
        print(io, "  Ploidy:") ; printstyled(io, " absent\n", color = :yellow)
    end
    print(io, "  Samples: ") ; printstyled(io, length(data.loci.name.pool), "\n", bold = true)
    print(io, "  Populations: ") ; printstyled(io, length(data.loci.population |> unique), "\n", bold = true)

    if "longitude" ∈ names(data.meta)
        if ismissing.(data.meta.longitude) |> all == true
            print(io, "  Coordinates:") ; printstyled(io, " absent\n", color = :yellow)
        elseif ismissing.(data.meta.longitude) |> all == false
            println(io, "  Coordinates: present")
        else
            println(io, "  Coordinates: present (", count(i -> i === missing, data.meta.longitude), " missing)")
        end
    else
        print(io, "  Coordinates:") ; printstyled(io, " absent\n", color = :yellow)
    end
    # add extracols section
end
