using Test,
      GeneticBitArrays

@testset "Construction" begin
  @test sum(DNASeq("AAAA").data, dims=2)[:] == [4;0;0;0]
  @test sum(DNASeq("CCCC").data, dims=2)[:] == [0;4;0;0]
  @test sum(DNASeq("GGGG").data, dims=2)[:] == [0;0;4;0]
  @test sum(DNASeq("TTTT").data, dims=2)[:] == [0;0;0;4]
  @test DNASeq("TTTT").data == RNASeq("UUUU").data
  @test length(DNASeq("ACGT")[2]) == length(DNASeq("C"))
end

@testset "Random" begin
  a = rand(RNASeq, Weights(fill(0.25, 4)), 1000)
  @test typeof(a) == RNASeq
  @test sum(a.data) == 1000
  @test all(sum(a.data, dims=1) .== 1)
end

@testset "Indexing" begin
  @test DNASeq("ACGT")[2] == DNASeq("C")
  b =  RNASeq("ACGU")
  b[1] = 'U'
  @test b[1] == RNASeq("U")
end

@testset "Errors" begin
  @test_throws ErrorException DNASeq("AAAU")
  @test_throws ErrorException RNASeq("AAAT")
  @test_throws ErrorException rand(RNASeq, Weights(fill(0.25, 3)), 1000)
end
