function length(x::T) where {T <: GeneticSeq}
  return size(x.data, 2)
end


function show(io::IO, x::T) where {T <: GeneticSeq}
  len = length(x)
  println(io, "$(len)nt DNA sequence")
  if len <= 26
    print(io, prod([onehot(T, x.data[:, i]) for i=1:len]))
  else
    print(io, prod([onehot(T, x.data[:, i]) for i=1:13]) *
              "..." *
              prod([onehot(T, x.data[:, i]) for i=len-13:len]))
  end
end


function getindex(x::T, i::Int64) where {T <: GeneticSeq}
  return T(x.data[:, i])
end


function setindex(x::T, a::BitArray{1}, i) where {T <: GeneticSeq}
  return x.data[:, i] = a
end


function setindex(x::T, a::Char, i) where {T <: GeneticSeq}
  return x.data[:, i] = onehot(T, a)
end


function ==(x::T, y::T) where {T <: GeneticSeq}
  return x.data == y.data
end
