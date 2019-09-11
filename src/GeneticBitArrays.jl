module GeneticBitArrays

import Base.show,
       Base.length,
       Base.getindex,
       Base.setindex,
       Base.==

abstract type GeneticSeq end

include("DNA.jl")
include("RNA.jl")
include("Common.jl")

export DNASeq, RNASeq

end # module
