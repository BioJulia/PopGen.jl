#=
This file contains the helper functions necessary for file import/export
of various file formats.
=#
"""
    find_ploidy(genotypes::T where T<:SubArray)
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes
of a sample and return the ploidy of the first non-missing locus.
"""
@inline function find_ploidy(genotypes::T where T<:SubArray)
    @inbounds for i in genotypes
        i !== missing && return Int8(length(i))
    end
end

@inline function find_ploidy(genotypes::T where T<:AbstractArray)
    @inbounds for i in genotypes
        i !== missing && return Int8(length(i))
    end
end

"""
    phase(loc::String, type::DataType, digit::Int)
Takes a String of numbers returns a typed locus appropriate for PopGen.jl as used in the
`genepop` and `csv` file parsers. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

## Examples
```
ph_locus = phase("128114", Int16, 3)
map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
"""
@inline function phase(loc::T, type::DataType, digit::Int) where T<:AbstractString
    loc == "-9" || loc == "0"^length(loc) && return missing
    phased = map(i -> parse(type, join(i)), Iterators.partition(loc, digit))
    sort!(phased)
    tupled = Tuple(phased)
    return tupled
end

"""
    phase_dip(loc::String, type::DataType, digit::Int)
A diploid-optimized variant of `phase()` that uses integer division to split the alleles
of a locus into a tuple of `type` Type. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

## Examples
```
ph_locus = phase("128114", Int16, 3)
map(i -> phase_dip(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase_dip(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
"""
@inline function phase_dip(loc::T, type::DataType, digit::Int) where T<:Signed
    loc == -9 || iszero(loc) && return missing
    units = 10^digit
    allele1 = loc ÷ units |> type
    allele2 = loc % units |> type
    return [allele1, allele2] |> sort |> Tuple
end

@inline function phase_dip(loc::Missing, type::DataType, digit::Int)
    return missing
end

@inline function phase_dip(loc::String, type::DataType, digit::Int)
    loc == "-9" || loc == "0"^length(loc) && return missing
    newloc = parse(Int32, loc)
    units = 10^digit
    allele1 = newloc ÷ units |> type
    allele2 = newloc % units |> type
    return [allele1, allele2] |> sort |> Tuple
end
