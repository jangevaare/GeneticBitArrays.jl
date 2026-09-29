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
  b = rand(DNASeq, fill(Weights(fill(0.25, 4)), 1000))
  @test typeof(b) == DNASeq
  @test sum(b.data) == 1000
  @test all(sum(b.data, dims=1) .== 1)
end

@testset "Indexing" begin
  for T in (DNASeq, RNASeq)
    first = T('A')
    second = T('A')
    first[1] = 'C'
    @test second == T("A")
    @test T('A') == T("A")
  end
  @test DNASeq("ACGT")[2] == DNASeq("C")
  c =  RNASeq("ACGU")
  c[1] = 'U'
  @test c[1] == RNASeq("U")
end

@testset "Errors" begin
  @test_throws ErrorException DNASeq("AAAU")
  @test_throws ErrorException RNASeq("AAAT")
  @test_throws ErrorException rand(RNASeq, Weights(fill(0.25, 3)), 1000)
end

@testset "String conversion" begin
  for (T, alphabet) in ((DNASeq, "NVHMDRWABSYCKGT-"),
                        (RNASeq, "NVHMDRWABSYCKGU-"))
    for text in ("", alphabet, repeat(alphabet, 100))
      sequence = T(text)
      @test String(sequence) == text
      @test convert(String, sequence) == text
      @test GeneticBitArrays.convert(String, sequence) == text
    end
    @test occursin(alphabet, sprint(show, T(alphabet)))
  end
end
