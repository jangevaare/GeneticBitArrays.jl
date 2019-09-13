struct RNASeq <: GeneticSeq
  data::BitArray{2}

  function RNASeq(x; checkinput::Bool=true)
    if checkinput
      if _checkinput(x)
        return new(x)
      else
        throw(ErrorException("Invalid input"))
      end
    else
      return new(x)
    end
  end
end


function onehot(::Type{RNASeq}, x::Char)
  if x     == 'A'
    return BitArray{1}([1, 0, 0, 0])
  elseif x == 'C'
    return BitArray{1}([0, 1, 0, 0])
  elseif x == 'G'
    return BitArray{1}([0, 0, 1, 0])
  elseif x == 'U'
    return BitArray{1}([0, 0, 0, 1])
  else
    throw(ErrorException("Invalid `Char`"))
  end
end


function onehotinv(::Type{RNASeq}, x::BitArray{1})
  return ['A', 'C', 'G', 'U'][x][1]
end


function RNASeq(seq::Vector{Char})
  x = BitArray{2}(fill(0, (4, length(seq))))
  for i in 1:length(seq)
    if seq[i]     == 'A'
      x[1, i] = 1
    elseif seq[i] == 'C'
      x[2, i] = 1
    elseif seq[i] == 'G'
      x[3, i] = 1
    elseif seq[i] == 'U'
      x[4, i] = 1
    else
      throw(ErrorException("Invalid `Char` at index $i"))
    end
  end
  return RNASeq(x, checkinput=false)
end


function RNASeq(x::String)
  return RNASeq(Vector{Char}(x))
end


function RNASeq(x::BitArray{1})
  return RNASeq(reshape(x, (4, 1)))
end
