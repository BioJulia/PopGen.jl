## test data available for use in PopGen

"""
    nancycats()
Returns a `PopObj` of corresponding "nancycats" dataset as featured in
the R package `adegenet`. This is microsatellite data corresponding to 4
populations.

Example:

nancy = nancycats()
"""
function nancycats()
    println("Downloading nancycats data")
    download("https://github.com/pdimens/PopGen.jl/raw/master/test/nancycats.gen","testdata")
    x = genepop("testdata", popsep = "pop", numpops = 4)
    rm("testdata")
    return x
end


"""
    gulfsharks()
Returns a `PopObj` of corresponding the Blacknose shark dataset as used in
Dimens et al. 2019. This is a mid-sized SNP dataset of 2213 SNPs across
212 individuals, within 7 populations.

Example:

sharks = gulfsharks()
"""
function gulfsharks()
    println("Downloading Blacknose shark data from Dimens et al. 2019")
    download("https://github.com/pdimens/PopGen.jl/raw/master/test/testdata.gen","testdata")
    x = genepop("testdata", numpops = 7)
    rm("testdata")
    return x
end
