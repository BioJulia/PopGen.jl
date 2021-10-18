export pairwisefst

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
    pairwisefst(data::PopData; method::String, by::String = "global", iterations::Int64)
Calculate pairwise FST between populations in a `PopData` object. Set `iterations` 
to a value greater than `0` to perform a single-tailed permutation test to obtain
P-values of statistical significance. Use `by = "locus"` to perform a locus-by-locus FST for
population pairs (iterations and significance testing ignored). 
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
function pairwisefst(data::PopData; method::String = "WC84", by::String = "global", iterations::Int64 = 0)
    if occursin("nei", lowercase(method))
        if by == "locus" 
            return _pairwise_Nei_lxl(data)
        else
            iterations > 0 && return _permuted_Nei(data, iterations)
            iterations == 0 && return _pairwise_Nei(data)
        end
    elseif occursin("wc", lowercase(method))
        if by == "locus" 
            throw(ArgumentError("locus-by-locus pairwise FST is not yet implemented for the WeirCockerham method." * feature_req()))
            #return _pairwise_WeirCockerham_lxl(data)
        else
            iterations > 0 && return _permuted_WeirCockerham(data, iterations)
            iterations == 0 && return _pairwise_WeirCockerham(data)
        end
    elseif occursin("hudson", lowercase(method))
        if by == "locus" 
            return _pairwise_Hudson_lxl(data)
        else
            iterations > 0 && return _permuted_Hudson(data, iterations)
            iterations == 0 && return _pairwise_Hudson(data)
        end
    else
        throw(ArgumentError("please use one of \"WC84\", \"Nei87\", or \"Hudson92\" methods"))
    end
end