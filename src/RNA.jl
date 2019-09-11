struct RNASeq <: GeneticSeq
  data::BitArray{2}
end


function onehot(::Type{RNASeq}, x::Char)
  return ['A', 'C', 'G', 'U'] .== Ref(x)
end


function onehot(::Type{RNASeq}, x::BitArray{1})
  return ['A', 'C', 'G', 'U'][x][1]
end


function RNASeq(seq::Vector{Char})
  x = BitArray(undef, (4, length(seq)))
  for i in 1:length(seq)
    x[:, i] = onehot(RNASeq, seq[i])
  end
  return RNASeq(x)
end


function RNASeq(x::String)
  return RNASeq(Vector{Char}(x))
end


function RNASeq(x::BitArray{1})
  return RNASeq(reshape(x, (4, 1)))
end
