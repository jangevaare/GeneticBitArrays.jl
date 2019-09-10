function onehotRNA(x::Char)
  return ['A', 'C', 'G', 'U'] .== Ref(x)
end


function onehotRNA(x::BitArray{1})
  return ['A', 'C', 'G', 'U'][x][1]
end

struct RNASeq
  data::BitArray{2}
end


function RNASeq(seq::Vector{Char})
  x = BitArray(undef, (4, length(seq)))
  for i in 1:length(seq)
    x[:, i] = onehotRNA(seq[i])
  end
  return RNASeq(x)
end


function RNASeq(x::String)
  return RNASeq(Vector{Char}(x))
end


function length(x::RNASeq)
  return size(x.data, 2)
end


function show(io::IO, x::RNASeq)
  len = length(x)
  println(io, "$(len)nt RNA sequence")
  if len <= 26
    print(io, prod([onehotRNA(x.data[:, i]) for i=1:len]))
  else
    print(io, prod([onehotRNA(x.data[:, i]) for i=1:13]) *
              "..." *
              prod([onehotRNA(x.data[:, i]) for i=len-13:len]))
  end
end
