#export

## experimental and not exported or documented!
function adjacency_matrix(data::PopData)
    data_loci = groupby(data.genodata, :locus)
    out_vec = Vector{Array{Int8,2}}(undef, length(data_loci))
    for (j,i) in enumerate(data_loci)
        uniq = unique(skipmissing(i.genotype))
        adj_mat = fill(Int8(0), length(samples(data)), length(uniq))
        for (j,k) in zip(i.genotype, eachrow(adj_mat))
            k .= Ref(j) .=== uniq 
        end
        out_vec[j] = adj_mat
    end
    return out_vec
end


"""
    multitest_missing(pvals::Vector{T}, method::String) where T <: Union{Missing, <:AbstractFloat}
Modification to `MultipleTesting.adjust` to include `missing` values in the
returned array. See MultipleTesting.jl docs for full more detailed information.

**Example**
```
julia> multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")

```
### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forwardstop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment
"""
@inline function multitest_missing(pvals::Vector{T}, method::String) where T <: Union{Missing, <:AbstractFloat}
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
        "forwardstop" => ForwardStop(),
        "fs" => ForwardStop(),
        "bc" => BarberCandes(),
    )
    p_copy = copy(pvals)
    p_copy[.!ismissing.(p_copy)] .= adjust(p_copy[.!ismissing.(p_copy)] |> Vector{Float64}, d[lowercase(method)])
    return p_copy
end