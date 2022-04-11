module  TestAlleleMatrices

using PopGen
using Test

cats = @nancycats ;

@testset "AlleleMatrices.jl" begin
    @testset "_allelematrix" begin
        @test size(PopGen._allelematrix(cats, by = "count", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, by = "frequency", missings = "missing", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, by = "frequency", missings = "zero", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, by = "frequency", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, scale = true)) == (237, 108)
        @test size(PopGen._allelematrix(cats, center = true)) == (237, 108)
        @test PopGen._allelematrix(cats, scale = true) != PopGen._allelematrix(cats, center = true)
        @test PopGen._allelematrix(cats, by = "count") != PopGen._allelematrix(cats)
        @test PopGen._allelematrix(cats, missings = "missing") != PopGen._allelematrix(cats, missings = "zero")
        @test PopGen._allelematrix(cats, missings = "zero") != PopGen._allelematrix(cats)
    end
end

end # module