export pairwise_fst

# generate a wrapper struct so we can nicely print the results
struct PairwiseFST
    results::DataFrame
    method::String
end

# pretty-printing of FST results
function Base.show(io::IO, data::PairwiseFST)
    show(
        io,
        data.results,
        show_row_number=false,
        rowlabel = Symbol(" "),
        eltypes = false,
        row_names = names(data.results),
        title = "Pairwise FST: " * data.method
    )
end


"""
    pairwise_fst(data::PopData; method::String, iterations::Int64)
Calculate pairwise FST between populations in a `PopData` object. Set `iterations` 
to a value greater than `0` to perform a single-tailed permutation test to obtain
P-values of statistical significance.
#### Methods:
- `"Hudson92"`: Hudson et al. (1992) method (only for biallelic data)
- `"Nei87"`: Nei (1987) method
- `"WC84"` : Weir & Cockerham (1984) method (default)

**Examples**
```julia
data = @nancycats
wc = pairwise_fst(data, method = "WC84")
wc_sig = pairwise_fst(data, iterations = 1000)
```
"""
function pairwise_fst(data::PopData; method::String = "WC84", iterations::Int64 = 0)
    if occursin("nei", lowercase(method))
        iterations > 0 && return _permuted_Nei(data, iterations)
        iterations == 0 && return _pairwise_Nei(data)
    elseif occursin("wc", lowercase(method))
        iterations > 0 && return _permuted_WeirCockerham(data, iterations)
        iterations == 0 && return _pairwise_WeirCockerham(data)
    elseif occursin("hudson", lowercase(method))
        iterations > 0 && return _permuted_Hudson(data, iterations)
        iterations == 0 && return _pairwise_Hudson(data)
    else
        throw(ArgumentError("please use one of \"WC84\" or \"Nei87\" methods"))
    end
end