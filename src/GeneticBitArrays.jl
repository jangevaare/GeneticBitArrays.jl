module GeneticBitArrays

import Base.getindex,
       Base.setindex,
       Base.show,
       Base.length

function onehotDNA(x::Char)
  return ['A', 'C', 'G', 'T'] .== Ref(x)
end


function onehotDNA(x::BitArray{1})
  return ['A', 'C', 'G', 'T'][x]
end


struct DNABitArray
  data::BitArray{2}

  function DNABitArray(seq::Vector{Char})
     x = BitArray(undef, (4, length(seq)))
    for i in 1:length(seq)
      x[onehotDNA(i), i] += 1
    end
    return new(x)
  end
end


function DNABitArray(x::String)
  return DNABitArray(Vector{Char}(x))
end


function length(x::DNABitArray)
  return size(x.data, 1)
end


function getindex(x::DNABitArray, i)
  return DNABitArray(x.data[i, :])
end


function show(io::IO, x::DNABitArray)
  if length(x) <= 26
    print(io, prod([onehotDNA(x[i]) for i=1:length(object)]))
  else
    print(io, prod([onehotDNA(x[i]) for i=1:13]) *
              "..." *
              prod([onehotDNA(x[i]) for i=length(x)-13:length(x)]))
  end
end

export DNABitArray

end # module
