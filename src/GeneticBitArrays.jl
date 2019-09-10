module GeneticBitArrays

import Base.show,
       Base.length

include("DNA.jl")
include("RNA.jl")

export DNASeq, RNASeq

end # module
