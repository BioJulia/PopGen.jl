export size, drop_monomorphic, drop_monomorphic!

## experimental and not exported or documented!
function adjacency_matrix(data::PopData)
    data_loci = groupby(data.loci, :locus)
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

#TODO change location in API docs and rename allele_pool?
#TODO replace alleles with universal Symbol type?
"""
    alleles(locus::T) where T<:GenoArray
Return an array of all the non-missing alleles of a locus.
"""
@inline function alleles(locus::T) where T<:GenoArray
    if all(ismissing.(locus))
        int_type = eltype(typeof(locus)) |> nonmissingtype |> eltype
        return Vector{Union{Missing, int_type}}(undef, length(locus))
    end
    alle_out = Base.Iterators.flatten(skipmissing(locus)) |> collect
end


"""
    allele_count(locus::T) where T<:GenoArray
Return the number of unique alleles present at a locus.
"""
@inline function allele_count(locus::T) where T<:GenoArray
    unique(locus) |> skipmissing |> Base.Iterators.flatten |> unique |> length
end


"""
    alleles(locus::T, miss::Bool = false) where T<:GenoArray
Return an array of all the non-missing alleles of a locus. Use the second positional
argument as `true` to include missing values.
"""
@inline function alleles(locus::T, miss::Bool) where T<:GenoArray
    int_type = eltype(typeof(locus)) |> nonmissingtype |> eltype
    if all(ismissing.(locus))
        return Vector{Union{Missing, int_type}}(undef, length(locus))
    end
    alle_out = Vector{Union{Missing, int_type}}(Base.Iterators.flatten(skipmissing(locus)) |> collect)
    if miss == true
        append!(alle_out, locus[locus .=== missing])
    end
    return alle_out
end


function Base.size(data::PopData)
    return (samples = size(data.meta)[1], loci = length(loci(data)))
end

"""
    unique_alleles(locus::T) where T<:GenoArray
Return an array of all the unique non-missing alleles of a locus.
"""
@inline function unique_alleles(locus::GenoArray)
    unique(alleles(locus))
end


"""
    convert_coord(coordinate::String)
Takes non-decimal-degree format as a `String` and returns it as a decimal degree
`Float32`. Can be broadcasted over an array of coordinate strings to convert them.
## Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
### Example
```
julia> convert_coord("-41 31.52")
-41.5253f0

julia> convert_coord.(["-41 31.52", "25 11:54S"])
2-element Array{Float32,1}:
-41.5253
-25.1983
```
"""
function convert_coord(coordinate::String)
    lowercase(coordinate) == "missing" && return missing
    coord_strip = replace(uppercase(coordinate), r"[NSEW]" => "")
    split_coord = parse.(Float32, split(coord_strip, r"\s|:"))
    split_coord[2] /= 60.0
    if length(split_coord) == 3
        split_coord[3] /= 3600.0
    end
    conv = mapreduce(abs, +, split_coord)
    # N + E are positive | S + W are negative
    if split_coord[1] < 0 || occursin(r"[SW]", uppercase(coordinate))
        # negative
        return round(conv * -1, digits = 4)
    else
        # positive
        return round(conv, digits = 4)
    end
end

function Base.copy(data::PopData)
    PopData(copy.([data.meta,data.loci])...)
end

"""
    count_nonzeros(x::AbstractVector{T}) where T<:Real
Return the number of non-zero values in a vector
"""
function count_nonzeros(x::AbstractVector{T}) where T<:Real
    mapreduce(!iszero, +, x)
end


"""
    drop_monomorphic(data::PopData)
Return a `PopData` object omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function drop_monomorphic(data::PopData)
    rm_loci = Vector{String}()
    for (loc, loc_sdf) in pairs(groupby(data.loci, :locus))
        length(unique(skipmissing(loc_sdf[:, :genotype]))) == 1 && push!(rm_loci, loc.locus)
    end

    if length(rm_loci) == 0
        return data
    elseif length(rm_loci) == 1
        @info "Dropping monomorphic locus " * rm_loci[1]
    else
        @info "Dropping $(length(rm_loci)) monomorphic loci" * "\n $rm_loci"
    end
    exclude(data, locus = rm_loci)
end


"""
    drop_monomorphic!(data::PopData)
Edit a `PopData` object in place by omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function drop_monomorphic!(data::PopData)
    rm_loci = Vector{String}()
    for (loc, loc_sdf) in pairs(groupby(data.loci, :locus))
        length(unique(skipmissing(loc_sdf[:, :genotype]))) == 1 && push!(rm_loci, loc.locus)
    end
    if length(rm_loci) == 0
        return data
    elseif length(rm_loci) == 1
        @info "Dropping monomorphic locus " * rm_loci[1]
    else
        @info "Dropping $(length(rm_loci)) monomorphic loci" * "\n $rm_loci"
    end
    exclude!(data, locus = rm_loci)
end

"""
    loci_dataframe(data::PopData)
Return a wide `DataFrame` of samples as columns, ommitting population information.

**Example**
```
julia> loci_dataframe(@nancycats)
9×237 DataFrame. Omitted printing of 232 columns
│ Row │ N215       │ N216       │ N217       │ N218       │ N219       │
│     │ Tuple…?    │ Tuple…?    │ Tuple…?    │ Tuple…?    │ Tuple…?    │
├─────┼────────────┼────────────┼────────────┼────────────┼────────────┤
│ 1   │ missing    │ missing    │ (135, 143) │ (133, 135) │ (133, 135) │
│ 2   │ (136, 146) │ (146, 146) │ (136, 146) │ (138, 138) │ (140, 146) │
│ 3   │ (139, 139) │ (139, 145) │ (141, 141) │ (139, 141) │ (141, 145) │
│ 4   │ (116, 120) │ (120, 126) │ (116, 116) │ (116, 126) │ (126, 126) │
│ 5   │ (156, 156) │ (156, 156) │ (152, 156) │ (150, 150) │ (152, 152) │
│ 6   │ (142, 148) │ (142, 148) │ (142, 142) │ (142, 148) │ (142, 148) │
│ 7   │ (199, 199) │ (185, 199) │ (197, 197) │ (199, 199) │ (193, 199) │
│ 8   │ (113, 113) │ (113, 113) │ (113, 113) │ (91, 105)  │ (113, 113) │
│ 9   │ (208, 208) │ (208, 208) │ (210, 210) │ (208, 208) │ (208, 208) │
```
"""
function loci_dataframe(data::PopData)
    unstack(select(data.loci, Not(:population)), :name, :genotype)[:, Not(:locus)]
end

#TODO make a SMatrix instead?
"""
    loci_matrix(data::PopData)
Return a matrix of genotypes with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> loci_matrix(@nancycats)
237×9 Array{Union{Missing, Tuple{Int16,Int16}},2}:
 missing     (136, 146)  (139, 139)  …  (199, 199)  (113, 113)  (208, 208)
 missing     (146, 146)  (139, 145)     (185, 199)  (113, 113)  (208, 208)
 (135, 143)  (136, 146)  (141, 141)     (197, 197)  (113, 113)  (210, 210)
 (133, 135)  (138, 138)  (139, 141)     (199, 199)  (91, 105)   (208, 208)
 (133, 135)  (140, 146)  (141, 145)     (193, 199)  (113, 113)  (208, 208)
 (135, 143)  (136, 146)  (145, 149)  …  (193, 195)  (91, 113)   (208, 208)
 (135, 135)  (136, 146)  (139, 145)     (199, 199)  (105, 113)  (208, 208)
 (135, 143)  (136, 146)  (135, 149)     (193, 197)  (91, 91)    (208, 212)
 (137, 143)  (136, 146)  (139, 139)     (197, 197)  (105, 113)  (208, 212)
 (135, 135)  (132, 132)  (141, 145)     (197, 197)  (91, 105)   (208, 208)
 (137, 141)  (130, 136)  (137, 145)  …  (193, 199)  (91, 91)    (182, 182)
 (129, 133)  (130, 136)  (135, 145)     (193, 199)  (91, 113)   (182, 208)
 ⋮                                   ⋱                          
 (133, 135)  (136, 136)  (135, 139)  …  (199, 199)  (113, 113)  (182, 182)
 (133, 141)  (136, 136)  (135, 139)     (197, 197)  (113, 113)  (182, 208)
 (133, 141)  (130, 146)  (141, 141)     (191, 199)  missing     (208, 208)
 (123, 133)  (138, 138)  (141, 145)     (191, 197)  missing     (208, 208)
 (123, 133)  (138, 138)  (139, 139)     (197, 199)  missing     (208, 208)
 (133, 141)  (136, 146)  (139, 139)  …  (197, 197)  missing     (208, 208)
 (133, 141)  (130, 136)  (139, 145)     (191, 199)  missing     (208, 208)
 (133, 141)  (136, 146)  (139, 145)     (199, 199)  missing     (208, 220)
 (133, 143)  (130, 130)  (135, 145)     (197, 197)  missing     (208, 208)
 (135, 141)  (136, 144)  (143, 143)     (191, 197)  (113, 117)  (208, 208)
 (137, 143)  (130, 136)  (135, 145)  …  (193, 199)  (113, 117)  (208, 208)
 (135, 141)  (130, 146)  (135, 139)     (197, 197)  missing     (208, 208)
 ```
"""
function loci_matrix(data::PopData)
    dims = size(data)
    sort_df = issorted(data.loci, [:name, :locus]) ? sort(data.loci, [:name, :locus]) : data.loci
    reshape(sort_df.genotype, (dims.samples, dims.loci)) |> collect
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


"""
    nonmissing(vec::T) where T<:AbstractArray
Convenience function to count the number of non-`missing` values
in a vector.
"""
@inline function nonmissing(vec::T) where T<:AbstractArray
    mapreduce(!ismissing, +, vec)
end


"""
    nonmissing(data::PopData, locus::String)
Convenience function to count the number of non-`missing` samples
at a locus.
"""
@inline function nonmissing(data::PopData, locus::String)
    data.loci[data.loci.locus .== locus, :genotype] |> nonmissing
end

"""
    nonmissings(vec1::AbstractVector, vec2::AbstractVector)
Return a vector of indices where neither input vectors have a `missing` value, i.e. an
intersection of the indices of their non-missing elements.
"""
@inline function nonmissings(vec1::T, vec2::T) where T <: AbstractVector
    mapreduce(i -> findall(!ismissing, i), intersect, (vec1, vec2))
end


"""
    pairwise_pairs(smp_names::Vector{T}) where T
Given a vector of some iterable, returns a vector of tuples of unique all x 
all combinations of pairs, excluding self-comparisons.

**Example**
```
julia> samps = ["red_1", "red_2", "blue_1", "blue_2"] ;

julia> pairwise_pairs(samps)
6-element Array{Tuple{String,String},1}:
 ("red_1", "red_2")
 ("red_1", "blue_1")
 ("red_1", "blue_2")
 ("red_2", "blue_1")
 ("red_2", "blue_2")
 ("blue_1", "blue_2")
```
"""
@inline function pairwise_pairs(smp_names::AbstractVector{T}) where T
    len = length(smp_names)
    [tuple(smp_names[i], smp_names[j]) for i in 1:len-1 for j in i+1:len]
end


"""
    partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
Like Base.Iterators.Partition, except you can apply arbitrary sizes to
partition the array by. The `steps` must add up to the total row length
of the array.

**Example**
```
julia> partitionmatrix(rand(20,5), [10,3,4,3]) .|> size
((10, 5), (3, 5), (4, 5), (3, 5))
```
"""
# solution brilliantly provided by @stevengj and @mcabbott on Slack and Discourse (https://discourse.julialang.org/t/is-there-a-simple-intuitive-way-to-partition-a-matrix-by-arbitrary-strides-like-i/55863)
function partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
    v = axes(array,1)
    v == 1:sum(steps) || error("Steps provided do not sum to length of the first dimension")
    i = firstindex(v)
    tmp = (view(v, i:(i+=s)-1) for s in steps)
    [view(array,r,:) for r in tmp]
end


"""
    phase(data::PopData)
Return a `Vector` of length `ploidy` composed of allele matrices with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> mtx = phase(@nancycats)
2-element Array{Array{Union{Missing, Int16},2},1}:
 [missing 136 … 113 208; missing 146 … 113 208; … ; 137 130 … 113 208; 135 130 … missing 208]
 [missing 146 … 113 208; missing 146 … 113 208; … ; 143 136 … 117 208; 141 146 … missing 208]

julia> mtx[1]
237×9 Array{Union{Missing, Int16},2}:
    missing  136  139  116         156  142  199  113         208
    missing  146  139  120         156  142  185  113         208
 135         136  141  116         152  142  197  113         210
 133         138  139  116         150  142  199   91         208
 133         140  141  126         152  142  193  113         208
 135         136  145  120         150  148  193   91         208
 135         136  139  116         152  142  199  105         208
 135         136  135  120         154  142  193   91         208
 137         136  139  116         150  142  197  105         208
 135         132  141  120         150  148  197   91         208
 137         130  137  128         152  142  193   91         182
 129         130  135  126         144  140  193   91         182
   ⋮                                      ⋮                   
 133         136  135     missing  146  142  199  113         182
 133         136  135     missing  150  142  197  113         182
 133         130  141     missing  148  142  191     missing  208
 123         138  141     missing  148  142  191     missing  208
 123         138  139     missing  150  142  197     missing  208
 133         136  139     missing  150  142  197     missing  208
 133         130  139     missing  152  142  191     missing  208
 133         136  139     missing  150  142  199     missing  208
 133         130  135     missing  148  142  197     missing  208
 135         136  143     missing  144  142  191  113         208
 137         130  135     missing  150  142  193  113         208
 135         130  135     missing  150  142  197     missing  208
```
"""
function phase(data::PopData)
    dims = size(data)
    ploidy = unique(data.meta.ploidy)
    ploidy = length(ploidy) != 1 ? error("Phasing will not work on mixed-ploidy samples") : ploidy[1]
    sort_df = issorted(data.loci, [:name, :locus]) ? sort(data.loci, [:name, :locus]) : data.loci
    matrices = map(j -> map(i -> ismissing(i) ? missing : i[j] , sort_df.genotype), 1:ploidy)
    map(i -> collect(reshape(i, (dims.samples, dims.loci))), matrices)
end



"""
    reciprocal(num::T) where T <: Signed
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.
"""
function reciprocal(num::T) where T <: Real
    !iszero(num) ? 1.0/float(num) : 0.0
end


"""
    reciprocal_sum(x::AbstractVector{T}) where T<:Real
Return the sum of the reciprocal values of `x`, skipping the `Inf` values
resulting from divide-by-zero errors.
"""
function reciprocal_sum(x::AbstractVector{T}) where T<:Real
    mapreduce(reciprocal, +, x)
end


# TODO add skip____ to docs/API/Utils.jl
"""
    skipnan(itr)
Return an iterator over the elements in `itr` skipping `Inf` and `-Inf` values. The returned
object can be indexed using indices of itr if the latter is indexable. Indices
corresponding to `Inf` values are not valid: they are skipped by keys and eachindex,   
and a MissingException is thrown when trying to use them. This is effectively `skipmissing`
for `Inf` and `-Inf` values.

Use collect to obtain an `Array` containing the non-`Inf` values in `itr`. Note that even  
if `itr` is a multidimensional array, the result will always be a `Vector` since it is not   
possible to remove `Inf`s while preserving dimensions of the input.
"""
function skipinf(itr)
    Iterators.filter(isfinite, itr)
end


"""
    skipnan(itr)
Return an iterator over the elements in `itr` skipping `NaN` values. The returned
object can be indexed using indices of itr if the latter is indexable. Indices
corresponding to `NaN` values are not valid: they are skipped by keys and eachindex,   
and a MissingException is thrown when trying to use them. This is effectively `skipmissing`
for `NaN` values.

Use collect to obtain an `Array` containing the non-`NaN` values in `itr`. Note that even  
if `itr` is a multidimensional array, the result will always be a `Vector` since it is not   
possible to remove `NaN`s while preserving dimensions of the input.
"""
function skipnan(itr) 
    Iterators.filter(!isnan, itr)
end

"""
    skipinfnan(itr)
Return an iterator over the elements in `itr` skipping `NaN`, `Inf` and `-Inf` values.
See the docstrings of `skipinf` and `skipnan` more details.
"""
function skipinfnan(itr)
    Iterators.filter(x -> (isfinite(x) & !isnan(x)), itr)
end


"""
    sim_pairs(data::Vector{String})
Takes a Vector of sample names and returns a Tuple of sample pairs, grouped by simulation
number. This is an internal function used for isolating sibship pairs from simulated shipship
pairs (via `PopGenSims.jl`) to perform `relatedness` estimates only on those pairs.

**Example**
```julia
julia> a = ["sim1_off1", "sim1_off2", "sim2_off1", "sim2_off2"] ;

julia> sim_pairs(a)
("sim1_off1", "sim1_off2")
("sim2_off1", "sim2_off2")
```
"""
function sim_pairs(data::Vector{String})
    n = length(data)
    isodd(n) && error("Expected an even number of samples, but got $n")
    Tuple.(Base.Iterators.partition(sort(data), 2))
end

"""
    Base.sort(x::NTuple{N,T}) where N where T <: Signed 
Sort the integers within a Tuple and return the sorted Tuple.
"""
function Base.sort(x::NTuple{N,T}) where N where T <: Signed 
    Tuple(sort(SVector(x)))
end


"""
    strict_shuffle(x::T) where T <: AbstractArray
Shuffle only the non-missing values of a Vector and return a copy of the vector,
keeping the `missing` values at their original locations.
Use `strict_shuffle!` to edit in-place instead of returning a copy.
"""
@inline function strict_shuffle(x::T) where T <: AbstractArray
    # get indices of where original missing are
    miss_idx = findall(i -> i === missing, x)
    out_vec = shuffle(Xoroshiro128Star(), x[.!ismissing.(x)])

    insert!.(Ref(out_vec), miss_idx, missing)
    return out_vec
end

"""
    strict_shuffle!(x::T)! where T <: AbstractArray
Shuffle only the non-missing values of a Vector, keeping the
`missing` values at their original locations. Use `strict_shuffle`
to return a copy instead of editing in-place.
"""
function strict_shuffle!(x::T) where T <: AbstractArray
    @inbounds shuffle!(Xoroshiro128Star(), @view x[.!ismissing.(x)])
    return x
end


"""
    generate_meta(data::DataFrame)
Given a genotype DataFrame formatted like `PopData.loci`, generates a corresponding
`meta` DataFrame. In other words, it creates the `.meta` part of `PopData` from the `.loci` part.

**Example**
```
julia> cats = @nancycats ;

julia> cats_nometa = cats.loci ;

julia> cats_meta = generate_meta(cats_nometa)
237×5 DataFrame
 Row │ name    population  ploidy  longitude  latitude 
     │ String  String      Int8    Float32?   Float32? 
─────┼─────────────────────────────────────────────────
   1 │ N215    1                2   missing   missing  
   2 │ N216    1                2   missing   missing  
   3 │ N217    1                2   missing   missing  
   4 │ N218    1                2   missing   missing  
   5 │ N219    1                2   missing   missing  
   6 │ N220    1                2   missing   missing  
   7 │ N221    1                2   missing   missing  
  ⋮  │   ⋮         ⋮         ⋮         ⋮         ⋮
 232 │ N295    17               2   missing   missing  
 233 │ N296    17               2   missing   missing  
 234 │ N297    17               2   missing   missing  
 235 │ N281    17               2   missing   missing  
 236 │ N289    17               2   missing   missing  
 237 │ N290    17               2   missing   missing  
                                       224 rows omitted
```
"""
function generate_meta(data::DataFrame)
    grp = groupby(data, :name)
    nms = map(z -> z.name, keys(grp))
    pops = map(z -> first(z.population), grp)
    ploids = map(z -> find_ploidy(z.genotype), grp)
    DataFrame(
        :name => nms,
        :population => pops,
        :ploidy => ploids,
        :longitude => Vector{Union{Missing, Float32}}(undef, (length(nms))),
        :latitude => Vector{Union{Missing, Float32}}(undef, (length(nms)))
    )
end