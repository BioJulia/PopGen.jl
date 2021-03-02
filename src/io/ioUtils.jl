#=
This file contains the helper functions necessary for file import/export
of various file formats.
=#

"""
    determine_marker(geno_parse::T, digits::Int) where T<:AbstractDataFrame
Return either `Int8` or `Int16` depending on largest allelic value in all genotypes
in the first 10 samples of an input DataFrame (or all the samples if less than 10 samples).
If the largest allele is 20 or greater, the marker will be considered a Microsatellite
and coded in `PopData` as `Int16`, and the opposite is true for SNPs. There's no
specific reason 100 was chosen other than it being a reasonable buffer for edge
cases since SNP data <= 4, and haplotyped data could be a bit higher. Even if the
microsatellite markers are coded incorrectly, there will be zero impact to performance,
and considering how few microsatellite markers are used in typical studies, the
effect on in-memory size should be negligible (as compared to SNPs).
"""
function determine_marker(geno_parse::T, digits::Int) where T<:AbstractDataFrame
    # get the total # columns
    total_col = size(geno_parse,2)
    # find the 25% cutoff
    if total_col > 11
        num_test_loc = (total_col - 2) ÷ 8
    else
        num_test_loc = total_col - 2
    end
    # remove everything else
    test_df = @view geno_parse[:, 3:num_test_loc]

    # isolate the largest allele value
    max_allele = map(eachcol(test_df)) do i
        phase.(i, Int16, digits)  |>
        skipmissing |> Base.Iterators.flatten |> maximum
    end |> maximum

    if max_allele <= 100
        return Int8
    else
        return Int16
    end
end

"""
    find_ploidy(genotypes::T where T<:SubArray)
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes
of a sample and return the ploidy of the first non-missing locus.
"""
@inline function find_ploidy(genotypes::T) where T<:GenoArray
    return Int8(length(first(skipmissing(genotypes))))
end


"""
    phase(loc::Union{String, Int}, type::DataType, digit::Int)
Takes a String of numbers or Integers and returns a typed locus appropriate for PopGen.jl as used in the
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
    loc == "-9" || iszero(parse(Int, loc)) && return missing
    phased = map(i -> parse(type, join(i)), Iterators.partition(loc, digit))
    sort!(phased)
    Tuple(phased)
end

phase(loc::Missing, type::DataType, digit::Int) = missing

@inline function phase(loc::T, type::DataType, digit::Int) where T<:Integer
    loc == -9 || iszero(loc) && return missing
    units = 10^digit
    allele1 = loc ÷ units |> type
    allele2 = loc % units |> type
    return [allele1, allele2] |> sort |> Tuple
end

"""
    unphase(geno::T; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) where T <: Genotype
Takes a `Genotype` (e.g. `(131, 94)`) and returns a string of concatenated
alleles padded with *n* number of zeroes, where *n* is given by `digits = `.
`missing` values are returned as either a string of 'digits × ploidy' zeroes (`miss = 0`)
or `"-9"` (`miss = -9`). The `ploidy` flag is only relevant for unphasing `missing` genotypes
and not used otherwise.

### Example
```
unphase((1,2,3,4), digits = 3)
"001002003004"

unphase(missing, digits = 2, ploidy = 2, miss = -9)
"-9"

unphase(missing, digits = 2, ploidy = 2, miss = 0)
"0000"
```
"""
function unphase(geno::T; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) where T <: Genotype
    join(map(i -> lpad(i, digits, "0"), geno))
end


function unphase(geno::Missing; digits::Int = 3, ploidy::Int, miss::Int = 0)
    if miss == 0
        return "0"^(digits * ploidy)
    else
        return "$miss"
    end
end