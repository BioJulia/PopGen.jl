"""
    richness(data::PopObj)
Calculates various allelic richness and returns vector of per locus allelic richness. To be called internally by
function calculating overall or per population richness either rarefied or not.

"""
function richness(data::PopObj)
    rich = []
    @inbounds for locus in eachcol(data.loci, false)
        unique_genos = unique(locus)
        unique_genos = collect(skipmissing(unique_genos))

        unique_alleles = vec(reduce(hcat, getindex.(unique_genos,i) for i in eachindex(unique_genos[1]))) |> unique

        rich_locus = length(unique_alleles)
        push!(rich, rich_locus)
    end
    return rich
end

nancycats() |> richness
