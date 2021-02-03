export quickstart, size, drop_monomorphic, drop_monomorphic!

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
        return Vector{Union{Missing, eltype(locus).b.types[1]}}(undef, length(locus))
    end
    alle_out = Base.Iterators.flatten(skipmissing(locus)) |> collect
end

"""
    alleles(locus::T, miss::Bool = false) where T<:GenoArray
Return an array of all the non-missing alleles of a locus. Use the second positional
argument as `true` to include missing values.
"""
@inline function alleles(locus::T, miss::Bool) where T<:GenoArray
    int_type = eltype(locus).b.types[1]
    if all(ismissing.(locus))
        return Vector{Union{Missing, int_type}}(undef, length(locus))
    end
    alle_out = Vector{Union{Missing, int_type}}(Base.Iterators.flatten(skipmissing(locus)) |> collect)
    if miss == true
        append!(alle_out, locus[locus .=== missing])
    end
    return alle_out
end


# TODO add to docs: Utils.jl API
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
    split_coord[2] /=60.0
    if length(split_coord) == 3
        split_coord[3] /=3600.0
    end
    conv = sum(abs.(split_coord))
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

##TODO add to docs API
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

#TODO add to docs API
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

````

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
    count(!ismissing, vec)
end

#TODO add to API docs
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
    intersect(map(i -> findall(!ismissing, i), (vec1, vec2))...)
end


#TODO add to docs API
"""
    pairwise_pairs(smp_names::Vector{String})
Given a vector of sample names, returns a vector of tuples of unique all x 
all combinations of sample pairs, excluding self-comparisons.

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
@inline function pairwise_pairs(smp_names::AbstractVector{String})
    [tuple(smp_names[i], smp_names[j]) for i in 1:length(smp_names)-1 for j in i+1:length(smp_names)]
end

#TODO add to docs API
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
    quickstart()
Prints helpful text of how to get started using PopGen.
"""
function quickstart()
    printstyled("\n        Quickstart for PopGen\n\n", bold = true)
    println("Documentation: https://pdimens.github.io/PopGen.jl/")
    println("Motivational(?) quote: ", motivational_quote())
    println("\nA few things things you can do to get started:")

    printstyled("\nLoad in data\n\n", color = :magenta)
    println("- read_from(filename; kwargs...)")
    println("- genepop(infile; kwargs...)  or similar file-specific importer")
    println("- use available @gulfsharks or @nancycats datasets")

    printstyled("\nExplore PopData\n\n", color = :blue)
    println("- populations(PopData) to view population information")
    println("- loci(PopData) to view locus names")
    println("- samples(PopData) to view sample names")
    println("- missing(PopData, by = ...) to view missing information")

    printstyled("\nManipulate PopData\n\n", color = :light_red)
    println("- populations!(PopData, ...) to rename populations")
    println("- locations!(PopData, ...) to add geographical coordinates")
    println("- exclude!(PopData, kwargs...) to selectively remove data")

    printstyled("\nAnalyses\n\n", color = :green)
    println("- richness(PopData) to calculate allelic richness")
    println("- allele_avg(PopData) to calculate average # of alleles")
    println("- summary(PopData) to calculate F-statistics, heterozygosity, etc.")
    println("- hwe_test(PopData) to test for Hardy-Weinberg Equilibrium")

    return
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
    sim_pairs(data::Vector{String})
Takes a Vector of sample names and returns a Tuple of sample pairs, grouped by simulation
number. This is an internal function used for isolating sibship pairs from simulated shipship
pairs (via `PopGenSims.jl`) to perform `relatedness` estimates only on those pairs.

**Example**
julia> a = ["sim1_off1", "sim1_off2", "sim2_off1", "sim2_off2"] ;

julia> sim_pairs(a)
("sim1_off1", "sim1_off2")
("sim2_off1", "sim2_off2")
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

#TODO add to API docs
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
    motivational_quote()
Randomly returns one of ~50 _somewhat_ motivational quotes.
"""
function motivational_quote()
# quotes retrieved on 06/11/2020 from https://www.coburgbanks.co.uk/blog/friday-funnies/50-funny-motivational-quotes/
    quotes =  [
    "\"The elevator to success is out of order. You’ll have to use the stairs, one step at a time.\" Joe Girard",
    "\"People often say that motivation doesn’t last. Well, neither does bathing – that’s why we recommend it daily.\" Zig Ziglar",
    "\"I always wanted to be somebody, but now I realise I should have been more specific.\" Lily Tomlin",
    "\"I am so clever that sometimes I don’t understand a single word of what I am saying.\" Oscar Wilde",
    "\"People say nothing is impossible, but I do nothing every day.\" Winnie the Pooh",
    "\"Life is like a sewer… what you get out of it depends on what you put into it.\" Tom Lehrer",
    "\"You can’t have everything. Where would you put it?\" Steven Wright",
    "\"Work until your bank account looks like a phone number.\" Unknown ",
    "\"Change is not a four letter word… but often your reaction to it is!\" Jeffrey Gitomer",
    "\"If you think you are too small to make a difference, try sleeping with a mosquito.\" Dalai Lama",
    "\"Bad decisions make good stories.\" Ellis Vidler",
    "\"I’ll probably never fully become what I wanted to be when I grew up, but that’s probably because I wanted to be a ninja princess.\" Cassandra Duffy",
    "\"When life gives you lemons, squirt someone in the eye.\" Cathy Guisewite",
    "\"A clear conscience is a sure sign of a bad memory.\" Mark Twain",
    "\"Well-behaved women seldom make history.\" Laurel Thatcher Ulrich",
    "\"I didn’t fail the test. I just found 100 ways to do it wrong.\" Benjamin Franklin",
    "\"I used to think I was indecisive, but now I’m not so sure.\" Unknown",
    "\"Don’t worry about the world coming to an end today. It’s already tomorrow in Australia.\" Charles Schulz",
    "\"Think like a proton. Always positive.\" Unknown",
    "\"Be happy – it drives people crazy.\" Unknown",
    "\"Optimist: someone who figures that taking a step backward after taking a step forward is not a disaster, it’s more like a cha-cha.\" Robert Brault",
    "\"The question isn’t who is going to let me, it’s who is going to stop me.\" Ayn Rand.",
    "\"You’re only given a little spark of madness. You mustn’t lose it.\" Robin Williams",
    "\"I am an early bird and a night owl… so I am wise and I have worms\" Michael Scott",
    "\"If you let your head get too big, it’ll break your neck.\" Elvis Presley",
    "\"The road to success is dotted with many tempting parking spaces.\" Will Rogers",
    "\"Live each day like it’s your second to the last. That way you can fall asleep at night.\" Jason Love",
    "\"Even a stopped clock is right twice every day. After some years, it can boast of a long series of successes.\" Marie Von Ebner-Eschenbach",
    "\"Honest criticism is hard to take, particularly from a relative, a friend, an acquaintance, or a stranger.\" Franklin P. Jones",
    "\"Opportunity is missed by most people because it is dressed in overalls and looks like work.\" Thomas Eddison",
    "\"A diamond is merely a lump of coal that did well under pressure.\" Unknown",
    "\"Nothing is impossible, the word itself says “I’m possible!\" Audrey Hepburn",
    "\"Friendship is like peeing on yourself: everyone can see it, but only you get the warm feeling that it brings.\" Robert Bloch",
    "\"Women who seek to be equal with men lack ambition.\" Marilyn Monroe",
    "\"By working faithfully eight hours a day you may eventually get to be boss and work twelve hours a day.\" Robert Frost",
    "\"The trouble with having an open mind, of course, is that people will insist on coming along and trying to put things in it.\" Terry Pratchett",
    "\"Age is of no importance unless you’re a cheese.\" Billie Burke",
    "\"When tempted to fight fire with fire, remember that the Fire Department usually uses water.\" Unknown",
    "\"Trying is the first step toward failure.\" Homer Simpson",
    "\"Happiness is just sadness that hasn’t happened yet.\" Unknown",
    "\"The best things in life are actually really expensive.\" Unknown",
    "\"A few harmless flakes working together can unleash an avalanche of destruction.\" Justin Sewell",
    "\"It could be that your purpose in life is to serve as a warning to others.\"  Ashleigh Brilliant",
    "\"If the world didn’t suck we’d all fly into space.\" Unknown",
    "\"Always remember that you are unique – just like everybody else.\" Unknown"
    ]
    return quotes[rand(1:length(quotes))]
end

#TODO add to docs API
"""
    generate_meta(data::DataFrame)
Given a genotype DataFrame formatted like `PopData.loci`, generates a corresponding
`meta` DataFrame. In other words, it creates the `.meta` part of `PopData` from the `.loci` part.

**Example**:
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