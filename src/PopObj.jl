"""
    PopObj(samples::DataFrame, loci::DataFrame)
The data struct used for the PopGen population genetics ecosystem. You are
STRONGLY discouraged from manually creating dataframes to pass into a PopObj,
and instead should use the provided genepop, csv, or vcf file importers.

- `samples` ::DataFrame individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names/numbers
    - `ploidy` ::Int8 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `loci` ::DataFrame loci and their genotypes
    - columns are named by loci
    - genotypes are Tuples of Int16 or Int16, arraged in order of `.samples.name`
"""
struct PopObj
    samples::DataFrame
    loci::DataFrame
    function PopObj(x::DataFrame, y::DataFrame)
        if sort(names(x)) != [:latitude, :longitude, :name, :ploidy, :population]
            error("Incorrect column names in samples dataframe.
            Columns should be: name, population, ploidy, longitude, latitude")
        end
        size(x,1) != size(y,1) && error("length mismatch of dataframes. samples: $(size(x,1)) | loci: $(size(y, 1))")
        typeof(x.name) != Array{String,1} && @error ":name values must be of type String"
        new(x,y)
    end
end

"""
    PopOpt(samples::DataFrame, loci::DataFrame)
An immutable version of a `PopObj`, used for under-the-hood processes in PopGen.jl.
"""
struct PopOpt
    samples::DataFrame
    loci::DataFrame
    PopOpt(x::PopObj) = new(x.samples, x.loci)
end
