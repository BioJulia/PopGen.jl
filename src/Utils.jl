export quickstart, size, drop_monomorphic, drop_monomorphic!

#TODO change location in API docs
"""
    alleles(locus::T; miss::Bool = false) where T<:GenoArray
Return an array of all the non-missing alleles of a locus. Use
`miss = true` to include missing values.
"""
@inline function alleles(locus::T; miss::Bool = false) where T<:GenoArray
    if miss==true
        Base.Iterators.flatten(locus) |> collect
    else
        @inbounds Base.Iterators.flatten(locus[.!ismissing.(locus)]) |> collect
    end
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
        length(unique(@view loc_sdf[:, :genotype])) == 1 && push!(rm_loci, loc.locus)
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
        length(unique(@view loc_sdf[:, :genotype])) == 1 && push!(rm_loci, loc.locus)
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
function nonmissing(vec::T) where T<:AbstractArray
    count(!ismissing, vec)
end

#TODO add to API docs
"""
    nonmissing(data::PopData, locus::String)
Convenience function to count the number of non-`missing` samples
at a locus.
"""
function nonmissing(data::PopData, locus::String)
    data.loci[data.loci.locus .== locus, :genotype] |> nonmissing
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
function pairwise_pairs(smp_names::AbstractVector{String})
    [tuple(smp_names[i], smp_names[j]) for i in 1:length(smp_names)-1 for j in i+1:length(smp_names)]
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
    println("- use available gulfsharks() or nancycats() datasets")

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
    out_vec = shuffle(x[.!ismissing.(x)])

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
    @inbounds shuffle!(@view x[.!ismissing.(x)])
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