#TODO add to docs
"""
    convert_coord(coordinate::string)
Takes a decimal minute format as a `String` and returns it as a decimal degree
`Float32`. Can be broadcasted over an array of decimal-minute strings to convert them.
## Formatting requirements
- Decimal Minutes: `"-11 43.11"` (must use space and be a `String`)
- **Must** use negative sign `-` instead of cardinal directions
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)

### Example
```
julia> convert_coord("-41 31.52")
-41.5253f0

julia> convert_coord.(["-41 31.52", "25 11.54"])
2-element Array{Float32,1}:
 -41.5253
  25.1923
```
"""
function convert_coord(coordinate::String)
    lowercase(coordinate) == "missing" && return missing
    coord_strip = replace(uppercase(coordinate), r"[NSEW]" => "")
    split_coord = split(coord_strip, " ")
    coord_degree = parse(Float32, split_coord[1])
    coord_minute = round(parse(Float32, split_coord[2])/60.0, digits = 4)
    # N + E are positive | S + W are negative
    if coord_degree < 0 || occursin(r"[SW]", coordinate)
        # if negative, subtract
        return coord_degree - coord_minute
    else
        # if positive, add
        return coord_degree +  coord_minute
    end
end

"""
    nonmissing(vec::T) where T<:AbstractArray
Convenience function to count the number of non-`missing` values
in a vector.
"""
function nonmissing(vec::T) where T<:AbstractArray
    count(!ismissing, vec)
end

"""
    reciprocal(num::T) where T <: Signed
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.
"""
function reciprocal(num::T) where T <: Real
    iszero(num) ? 1.0/Float64(num) : 0
end


"""
    multitest_missing(pvals::Array{Float64,1}, correction::String)
Modification to `MultipleTesting.adjust` to include `missing` values in the
returned array. Missing values are first removed from the array, the appropriate
correction made, then missing values are re-added to the array at their original
positions. See MultipleTesting.jl docs for full more detailed information.
#### example
`multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")`

### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forward stop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment
"""
@inline function multitest_missing(pvals::Vector{T}, correction::String) where T <: Union{Missing, <:AbstractFloat}
    # get indices of where original missing are
    miss_idx = findall(i -> i === missing, pvals)

    # make seperate array for non-missing P vals
    p_no_miss = skipmissing(pvals) |> collect

    # make a dict of all possible tests and their respective functions
    d = Dict(
        "bonferroni" => Bonferroni(),
        "holm" => Holm(),
        "hochberg" => Hochberg(),
        "bh" => BenjaminiHochberg(),
        "by" => BenjaminiYekutieli(),
        "bl" => BenjaminiLiu(),
        "hommel" => Hommel(),
        "sidak" => Sidak(),
        "forward stop" => ForwardStop(),
        "fs" => ForwardStop(),
        "bc" => BarberCandes(),
    )

    correct = adjust(p_no_miss, d[lowercase(correction)]) |> Vector{Union{Missing, Float64}}

    # re-add missing to original positions
    @inbounds for i in miss_idx
        @inbounds insert!(correct, i, missing)
    end
    return correct
end
