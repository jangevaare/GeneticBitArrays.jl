struct DNASeq <: GeneticSeq
  data::BitArray{2}
end


function onehot(::Type{DNASeq}, x::Char)
  return ['A', 'C', 'G', 'T'] .== Ref(x)
end


function onehot(::Type{DNASeq}, x::BitArray{1})
  return ['A', 'C', 'G', 'T'][x][1]
end


function DNASeq(seq::Vector{Char})
  x = BitArray(undef, (4, length(seq)))
  for i in 1:length(seq)
    x[:, i] = onehot(DNASeq, seq[i])
  end
  return DNASeq(x)
end


function DNASeq(x::String)
  return DNASeq(Vector{Char}(x))
end


function DNASeq(x::BitArray{1})
  return DNASeq(reshape(x, (4, 1)))
end
