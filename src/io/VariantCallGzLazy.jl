using .GZip

# if GZip is loaded in, overwrite openvcf to this method
function openvcf(infile::String)
    if endswith(infile, ".gz")
        return GZip.open(infile, "r")
    else
        return open(infile, "r")
    end
end