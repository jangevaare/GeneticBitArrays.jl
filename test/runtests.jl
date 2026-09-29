using Random
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

@testset "Explicit random generators" begin
  generators = isdefined(Random, :Xoshiro) ? (MersenneTwister, Random.Xoshiro) : (MersenneTwister,)
  w = Weights([1, 2, 3, 4])
  for R in generators, T in (DNASeq, RNASeq), args in ((w, 100), (fill(w, 100),))
    rng = R(42)
    reference = R(42)
    a = rand(rng, T, args...)
    @test a == rand(reference, T, args...)
    @test rand(rng) == rand(reference)
    @test rand(R(43), T, args...) != a
    @test all(sum(a.data, dims=1) .== 1)
    Random.seed!(123)
    expected = rand()
    Random.seed!(123)
    rand(R(42), T, args...)
    @test rand() == expected
    Random.seed!(123)
    default = rand(T, args...)
    Random.seed!(123)
    @test rand(T, args...) == default
  end
  for T in (DNASeq, RNASeq)
    @test length(rand(MersenneTwister(42), T, w, 0)) == 0
    @test length(rand(MersenneTwister(42), T, typeof(w)[])) == 0
    @test rand(MersenneTwister(42), T, Weights([0, 1, 0, 0]), 4) == T("CCCC")
    @test rand(MersenneTwister(42), T, [Weights([1, 0, 0, 0]), Weights([0, 0, 1, 0])]) == T("AG")
    @test_throws ErrorException rand(MersenneTwister(42), T, Weights([1, 1, 1]), 4)
    @test_throws ErrorException rand(MersenneTwister(42), T, [w, Weights([1, 1, 1])])
  end
end
