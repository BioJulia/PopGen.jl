"""
    AbstractType PopObj
Generic AbstractType for use in PopGen.jl
"""
abstract type PopObj end

"""
    PopData(meta::IndexedTable, loci::IndexedTable)
The data struct used for the PopGen population genetics ecosystem. You are
STRONGLY discouraged from manually creating tables to pass into a PopObj,
and instead should use the provided genepop, csv, or vcf file importers.

- `meta` ::IndexedTable individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names
    - `ploidy` ::Int8 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `loci` ::IndexedTable Long-format table of sample genotype records
    - name` ::CategoricalString the individual/sample names
    - `population`::CategoricalString population name
    - `locus` ::CategoricalString of locus name
    - `genotype` Tuple of Int8 or Int16 depending on SNP or microsatellite
"""
struct PopData <: PopObj
    meta::IndexedTable
    loci::IndexedTable
end


function Base.show(io::IO, catst::T) where T<:CategoricalString
    print(io, "\"$(catst)\"")
end
