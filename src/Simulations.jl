 
"""
    simulate(data::PopData; n::Int = 100)

Simulate `n` number of individuals (default: `100`) per population using per-population
allele frequencies derived from a `PopData` object. Returns a new `PopData` object.

**Example**
```julia
cats = nancycats();

julia> sims = simulate(x , n = 100)
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 1700
  Loci: 9
  Populations: 17
  Coordinates: absent

  julia> sims.meta
  1700×5 DataFrame
  │ Row  │ name     │ population │ ploidy │ longitude │ latitude │
  │      │ String   │ String     │ Int64  │ Missing   │ Missing  │
  ├──────┼──────────┼────────────┼────────┼───────────┼──────────┤
  │ 1    │ sim_1    │ 1          │ 2      │ missing   │ missing  │
  │ 2    │ sim_2    │ 1          │ 2      │ missing   │ missing  │
  │ 3    │ sim_3    │ 1          │ 2      │ missing   │ missing  │
  │ 4    │ sim_4    │ 1          │ 2      │ missing   │ missing  │
  ⋮
  │ 1696 │ sim_1696 │ 17         │ 2      │ missing   │ missing  │
  │ 1697 │ sim_1697 │ 17         │ 2      │ missing   │ missing  │
  │ 1698 │ sim_1698 │ 17         │ 2      │ missing   │ missing  │
  │ 1699 │ sim_1699 │ 17         │ 2      │ missing   │ missing  │
  │ 1700 │ sim_1700 │ 17         │ 2      │ missing   │ missing  │  

  julia> sims.loci
  15300×4 DataFrame
  │ Row   │ name     │ population │ locus  │ genotype   │
  │       │ String   │ String     │ String │ Tuple…?    │
  ├───────┼──────────┼────────────┼────────┼────────────┤
  │ 1     │ sim_1    │ 1          │ fca8   │ (135, 135) │
  │ 2     │ sim_1    │ 1          │ fca23  │ (132, 140) │
  │ 3     │ sim_1    │ 1          │ fca43  │ (139, 139) │
  │ 4     │ sim_1    │ 1          │ fca45  │ (126, 126) │
  ⋮
  │ 15297 │ sim_1700 │ 17         │ fca78  │ (142, 142) │
  │ 15298 │ sim_1700 │ 17         │ fca90  │ (199, 199) │
  │ 15299 │ sim_1700 │ 17         │ fca96  │ (113, 113) │
  │ 15300 │ sim_1700 │ 17         │ fca37  │ (208, 208) │
```
"""
function simulate(data::PopData; n::Int = 100)
    length(unique(data.meta.ploidy)) != 1 && error("Simulations do not work on mixed-ploidy data (yet)")
    ploidy = first(unique(data.meta.ploidy))
    pops = unique(data.meta.population)
    npops = length(pops)
    nloci = length(data.loci.locus.pool)

    # instantiate output df
    simnames = reduce(vcat, fill.(["sim_" * "$i" for i in 1:(n*npops)], nloci))
    popnames = reduce(vcat, fill.(pops, (nloci * n)))
    locinames = reduce(vcat, fill.(Ref(unique(data.loci.locus)), (n * npops)))
    geno_out = DataFrame(:name => simnames, :population => popnames, :locus => locinames, :genotype => similar(data.loci.genotype, length(locinames)))

    # generate allele freqs per population
    gdf = groupby(data.loci, [:population, :locus])
    freqs = DataFrames.combine(
        gdf,
        :genotype => allele_freq => :frq
    )
    # create new genotypes
    transform!(freqs, :frq => (i -> sample_locus.(i,n,ploidy)) => :frq)
    
    # populate out df
    out_gdf = groupby(geno_out, :population)
    geno_gdf = groupby(freqs, :population)
    for pop in pops
        out_gdf[(population = pop,)].genotype .= reduce(hcat, geno_gdf[(population = pop,)].frq) |> permutedims |> vec
    end

    # regenerate meta info
    meta_df = DataFrames.combine(
        groupby(geno_out, :name),
        :name => first => :name,
        :population => first => :population
    )
    meta_df[!, :ploidy] .= 2
    meta_df[!, :longitude] .= missing
    meta_df[!, :latitude] .= missing

    PopData(meta_df, geno_out)
end


"""
    sample_locus(locus::Dict, n::Int, ploidy::Signed)

Internal function used by `simulate` to take a `Dict` of alleles => frequencies of a locus and return
`n` number of genotypes (n_alleles = `ploidy`) by using weighted sampling of the
allele-frequency pairs. 

**Example**
```julia
d = Dict(
  133 => 0.125,
  135 => 0.5625,
  143 => 0.25,
  137 => 0.0625
  )

julia> sample_locus(d, 3, 2)
5-element Array{Tuple{Int16,Int16},1}:
 (133, 135)
 (135, 135)
 (143, 137)

julia> sample_locus(d, 3, 3)
5-element Array{Tuple{Int16,Int16,Int16},1}:
 (135, 135, 133)
 (143, 135, 133)
 (137, 135, 135)
```
"""
function sample_locus(locus::Dict, n::Int, ploidy::Signed)
    isempty(locus) && return fill(missing, n)
    k,v = collect(keys(locus)), collect(values(locus))
    alleles = [sample(k, Weights(v), n) for i in 1:ploidy]
    Tuple.(sort.((eachrow(hcat(alleles...)))))
end