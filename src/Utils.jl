# TODO add to API docs
"""
    nonmissing(vec::T) where T<:AbstractArray
Convenience function to count the number of non-`missing` values
in a vector.
"""
function nonmissing(vec::T) where T<:AbstractArray
    count(!ismissing, vec)
end

# TODO add to API docs
"""
    reciprocal(num::T) where T <: Signed
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.
"""
function reciprocal(num::T) where T <: Real
    num != 0 ? 1/num : 0
end


# TODO replace location in docs
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
