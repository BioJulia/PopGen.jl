using .GZip

# if GZip is loaded in, overwrite openvcf to this method
function openvcf(infile::String)
    if endswith(infile, ".vcf") || endswith(infile, ".bcf")
        return open(infile, "r")
    elseif endswith(infile, ".vcf.gz") || endswith(infile, ".bcf.gz")
        return GZip.open(infile, "r")
    else
        throw(ArgumentError("The filename must end with .vcf/.bcf or .vcf.gz/.bcf.gz"))
    end
end
