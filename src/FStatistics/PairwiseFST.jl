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


#TODO add methods and complete docstring
"""
    pairwise_fst(data::PopData; method::String)
Calculate pairwise FST between populations in a `PopData` object.
#### Methods:
- `"Nei87"`: Nei (1987) method
- `"WC84"` : Weir-Cockerham (1984) method (default)
"""
function pairwise_fst(data::PopData; method::String = "WC")
    occursin("nei", lowercase(method)) && return _pairwise_Nei(data)
    occursin("wc", lowercase(method)) && return _pairwise_WeirCockerham(data)
end