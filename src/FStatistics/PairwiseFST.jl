# generate a wrapper struct so we can nicely print the results
"""
```julia
struct PairwiseFST
    results::DataFrame
    method::String
```
A convenience data structure which stores the `results` and `method` of a `pairwisefst` analysis.
The object pretty-prints a bit more information to the console, especially when doing a global pairwise FST. 
"""
struct PairwiseFST
    results::DataFrame
    method::String
end

# pretty-printing of FST results
function Base.show(io::IO, data::PairwiseFST)
    issymmetrical = size(data.results, 1) == size(data.results, 2)
    rwnames = issymmetrical ? names(data.results) : nothing
    show(
        io,
        round.(data.results, digits = 4),
        show_row_number = false,
        rowlabel = Symbol(" "),
        eltypes = false,
        row_names =  rwnames,
        title = "Pairwise FST: " * data.method
    )
end


#TODO update the docs with this implementation
"""
    pairwisefst(data::PopData; method::Function, by::String = "global", iterations::Int64)
Calculate pairwise FST between populations in a `PopData` object. Set `iterations` 
to a value greater than `0` to perform a single-tailed permutation test to obtain
P-values of statistical significance. Use `by = "locus"` to perform a locus-by-locus FST for
population pairs (iterations and significance testing ignored). Returns a `PairwiseFST` object,
stores a `DataFrame` of the `results`, along with the `method` used to obtain the estimates. 
#### Methods:
- `AMOVA`: AMOVA-based (Excoffier, 1992) FST as described by Bird et al. (2011) (default)
- `Hudson`: Hudson et al. (1992) method (only for biallelic data)
- `Nei`: Nei (1987) method
- `WeirCockerham` : Weir & Cockerham (1984) method

**Examples**
```julia
data = @nancycats
wc = pairwise_fst(data, method = WeirCockerham)
wc_sig = pairwise_fst(data, iterations = 1000)
```
"""
function pairwisefst(data::PopData; method::Function = AMOVA, by::String = "global", iterations::Int64 = 0)
    # sanity checks
    mth = Symbol(method)
    by = lowercase(by)
    if mth ∉ [:Hudson, :Nei, :WeirCockerham, :AMOVA]
        throw(ArgumentError("The \`method\` keyword argument ($method) is not recognized and must be one of Hudson, Nei, AMOVA, or WeirCockerham. See ?pairwisefst for usage information"))
    elseif mth == :Hudson
        isbiallelic(data) || throw(ArgumentError("Data must be biallelic to use the Hudson estimator"))
    end

    if by == "locus"
        if mth ∈ [:Hudson, :Nei]
            return _fst_lxl(data, Val(mth))
        else
            throw(ArgumentError("Method $mth is not (yet) supported for by-locus calculations."))
        end
    elseif by == "global"
        if iterations == 0
            _pairwisefst(data, Val(mth))
        else
            _fst_permutation(data, method, iterations)
        end
    else
        throw(ArgumentError("The keyword argument \`by\` must be either \"global\" or \"locus\""))
    end
end