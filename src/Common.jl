function rand(::Type{T}, w::W, n::Int64) where {T <: GeneticSeq, W<: Weights}
  x = BitArray(fill(0, (4, n)))
  @simd for i = 1:n
    @inbounds x[sample(1:4, w), i] = 1
  end
  return T(x)
end


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


function setindex!(x::T, a::BitArray{1}, i) where {T <: GeneticSeq}
  return x.data[:, i] = a
end


function setindex!(x::T, a::Char, i) where {T <: GeneticSeq}
  return x.data[:, i] = onehot(T, a)
end


function ==(x::T, y::T) where {T <: GeneticSeq}
  return x.data == y.data
end
