module  TestClustering

using PopGen
using Test
using Clustering
using DataFrames
using Distances

cats = @nancycats;

@testset "Clustering.jl" begin
    @testset "k-means" begin
        @test kmeans(cats, k = 2) isa Clustering.KmeansResult
        @test kmeans(cats, k = 2, matrixtype = :freq) isa Clustering.KmeansResult
        @test kmeans(cats, k = 2, iterations = 30) isa Clustering.KmeansResult
    end

    @testset "k-medoids" begin
        @test kmedoids(cats, k = 2) isa Clustering.KmedoidsResult
        @test kmedoids(cats, k = 2, distance = sqeuclidean) isa Clustering.KmedoidsResult
        @test kmedoids(cats, k = 2, distance = sqeuclidean, matrixtype = :freq) isa Clustering.KmedoidsResult
    end

    @testset "h-clust" begin
        hcout = hclust(cats)
        @test hcout isa Hclust
        @test hclust(cats, matrixtype = :freq) isa Hclust
        @test hclust(cats, distance = sqeuclidean) isa Hclust
        @test cutree(cats, hcout, krange = 2:5) isa DataFrame
    end

    @testset "fuzzy-c" begin
        @test fuzzycmeans(cats, c = 2) isa FuzzyCMeansResult
        @test fuzzycmeans(cats, c = 2, fuzziness = 3) isa FuzzyCMeansResult
        @test fuzzycmeans(cats, c = 2, iterations = 30) isa FuzzyCMeansResult
        @test fuzzycmeans(cats, c = 2, matrixtype = :freq) isa FuzzyCMeansResult
    end

    @testset "dbscan" begin
        @test dbscan(cats, radius = 0.5, minpoints = 2) isa DbscanResult
        @test dbscan(cats, radius = 0.5, minpoints = 2, distance = sqeuclidean) isa DbscanResult
        @test dbscan(cats, radius = 0.5, matrixtype = :freq) isa DbscanResult
    end
end

end #module