function _checkinput(x::BitArray{2})
  return size(x, 1) == 4 && all(sum(x, dims=1) .== 1)
end


function rand(::Type{T}, w::W, n::Int64; checkinput::Bool=true) where {T <: GeneticSeq, W<: Weights}
  if checkinput && length(w) != 4
    throw(ErrorException("Invalid sampling weights for $T generation"))
  end
  x = BitArray(fill(0, (4, n)))
  @simd for i = 1:n
    @inbounds x[sample(1:4, w), i] = 1
  end
  return T(x, checkinput=false)
end


function rand(::Type{T}, wv::Vector{W}; checkinput::Bool=true) where {T <: GeneticSeq, W<: Weights}
  len = length(wv)
  x = BitArray(fill(0, (4, len)))
  @simd for i = 1:len
    if checkinput && length(wv[i]) != 4
      throw(ErrorException("Invalid sampling weights for $(i)th nt in $T generation"))
    end
    @inbounds x[sample(1:4, wv[i]), i] = 1
  end
  return T(x, checkinput=false)
end


function length(x::T) where {T <: GeneticSeq}
  return size(x.data, 2)
end


function show(io::IO, x::T) where {T <: GeneticSeq}
  len = length(x)
  println(io, "$(len)nt DNA sequence")
  if len <= 26
    print(io, prod([onehotinv(T, x.data[:, i]) for i=1:len]))
  else
    print(io, prod([onehotinv(T, x.data[:, i]) for i=1:13]) *
              "..." *
              prod([onehotinv(T, x.data[:, i]) for i=len-13:len]))
  end
end


function getindex(x::T, i::Int64) where {T <: GeneticSeq}
  return T(x.data[:, i])
end


function setindex!(x::T, a::BitArray{1}, i::Int64) where {T <: GeneticSeq}
  if length(a) == 4 && sum(a) == 1
    return x.data[:, i] = a
  else
    throw(ErrorException("Invalid input"))
  end
end


function setindex!(x::T, a::Char, i::Int64) where {T <: GeneticSeq}
  return x.data[:, i] = onehot(T, a)
end


function ==(x::T, y::T) where {T <: GeneticSeq}
  return x.data == y.data
end
