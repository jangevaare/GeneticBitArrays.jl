module GeneticBitArrays

import Base.show,
       Base.length,
       Base.getindex,
       Base.setindex!,
       Base.==,
       Base.rand

import StatsBase.Weights,
       StatsBase.sample

abstract type GeneticSeq end

include("DNA.jl")
include("RNA.jl")
include("Common.jl")

export DNASeq, RNASeq, Weights

end # module
