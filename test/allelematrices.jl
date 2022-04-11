module  TestAlleleMatrices

using PopGen
using Test

cats = @nancycats ;

@testset "AlleleMatrices.jl" begin
    @testset "_allelematrix" begin
        @test size(_allelematrix(cats, by = "count", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(_allelematrix(cats, by = "frequency", missings = "missing", scale = false, center = false)) == (237, 108)
        @test size(_allelematrix(cats, by = "frequency", missings = "zero", scale = false, center = false)) == (237, 108)
        @test size(_allelematrix(cats, by = "frequency", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(_allelematrix(cats, scale = true)) == (237, 108)
        @test size(_allelematrix(cats, center = true)) == (237, 108)
        @test _allelematrix(cats, scale = true) != _allelematrix(cats, center = true)
        @test _allelematrix(cats, by = "count") != _allelematrix(cats)
        @test _allelematrix(cats, missings = "missing") != _allelematrix(cats, missings = zero)
        @test _allelematrix(cats, missings = zero) != _allelematrix(cats)
    end
end

end # module