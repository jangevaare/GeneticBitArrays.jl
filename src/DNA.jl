function onehotDNA(x::Char)
  return ['A', 'C', 'G', 'T'] .== Ref(x)
end


function onehotDNA(x::BitArray{1})
  return ['A', 'C', 'G', 'T'][x][1]
end

struct DNASeq
  data::BitArray{2}
end


function DNASeq(seq::Vector{Char})
  x = BitArray(undef, (4, length(seq)))
  for i in 1:length(seq)
    x[:, i] = onehotDNA(seq[i])
  end
  return DNASeq(x)
end


function DNASeq(x::String)
  return DNASeq(Vector{Char}(x))
end


function length(x::DNASeq)
  return size(x.data, 2)
end


function show(io::IO, x::DNASeq)
  len = length(x)
  println(io, "$(len)nt DNA sequence")
  if len <= 26
    print(io, prod([onehotDNA(x.data[:, i]) for i=1:len]))
  else
    print(io, prod([onehotDNA(x.data[:, i]) for i=1:13]) *
              "..." *
              prod([onehotDNA(x.data[:, i]) for i=len-13:len]))
  end
end
