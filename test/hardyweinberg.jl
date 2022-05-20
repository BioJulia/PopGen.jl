module  TestHardyWeinberg

using PopGen
using DataFrames
using Test

cats = @nancycats;
testarray = cats.genodata.genotype[1:10]

@testset "Hardy Weinberg" begin
    @test PopGen._chisqlocus(testarray) == (5.135802469135802, 6, 0.5265173055755079)
    tmp = hwetest(cats)
    @test tmp isa DataFrame
    @test size(tmp) == (9,4)
    tmp = hwetest(cats, correction = "bh")
    @test size(tmp) == (9,5)
    tmp = hwetest(cats, by = "population", correction = "bh")
    @test size(tmp) == (153,6)
end

end #module