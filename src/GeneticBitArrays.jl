module GeneticBitArrays

  import Base.show,
         Base.length,
         Base.getindex,
         Base.setindex!,
         Base.==,
         Base.rand,
         StatsBase.Weights,
         StatsBase.sample,
         StaticArrays.SVector

  abstract type GeneticSeq end

  struct DNASeq <: GeneticSeq
    data::BitArray{2}

    function DNASeq(x::BitArray{2}; checkinput::Bool=true)
      if checkinput
        if size(x, 1) == 4
          return new(x)
        else
          throw(ErrorException("Invalid input, must be `BitArray{4,n}`"))
        end
      else
        return new(x)
      end
    end
  end

  struct RNASeq <: GeneticSeq
    data::BitArray{2}

    function RNASeq(x::BitArray{2}; checkinput::Bool=true)
      if checkinput
        if size(x, 1) == 4
          return new(x)
        else
          throw(ErrorException("Invalid input, must be `BitArray{4,n}`"))
        end
      else
        return new(x)
      end
    end
  end

  function length(x::GeneticSeq)
    return size(x.data, 2)
  end

  function getindex(x::T, i) where {T <: GeneticSeq}
    return T(x.data[:, i])
  end

  const _bitslookup = SVector{16, BitArray{1}}([1;1;1;1], [1;1;1;0],
                                               [1;1;0;1], [1;1;0;0],
                                               [1;0;1;1], [1;0;1;0],
                                               [1;0;0;1], [1;0;0;0],
                                               [0;1;1;1], [0;1;1;0],
                                               [0;1;0;1], [0;1;0;0],
                                               [0;0;1;1], [0;0;1;0],
                                               [0;0;0;1], [0;0;0;0])

  const _dnacharlookup = SVector{16, Char}('N', 'V', 'H', 'M',
                                           'D', 'R', 'W', 'A',
                                           'B', 'S', 'Y', 'C',
                                           'K', 'G', 'T', '-')

  const _rnacharlookup = SVector{16, Char}('N', 'V', 'H', 'M',
                                           'D', 'R', 'W', 'A',
                                           'B', 'S', 'Y', 'C',
                                           'K', 'G', 'U', '-')

  const _lookup(::Type{DNASeq}) = _dnacharlookup
  const _lookup(::Type{RNASeq}) = _rnacharlookup
  const _seq(::Type{DNASeq}) = "DNA"
  const _seq(::Type{RNASeq}) = "RNA"

  function _bitarray(::Type{T}, x::Union{String, Vector{Char}}) where {T <: GeneticSeq}
    l = _lookup(T)
    s = BitArray{2}(undef, (4, length(x)))
    for i = eachindex(x)
      ind = findfirst(l .== x[i])
      if isnothing(ind)
        throw(ErrorException("Unrecognized $(_seq(T)) `Char` $(x[i]) at index $i"))
      else
        s[:, i] = _bitslookup[ind]
      end
    end
    return s
  end

  function _bitarray(::Type{T}, x::Char) where {T <: GeneticSeq}
    ind = findfirst(_lookup(T) .== x)
    if isnothing(ind)
      @error "Unrecognized $(_seq(T)) `Char` $x"
    end
    return _bitslookup[ind]
  end

  function _bitarray(::Type{T}, x::BitArray{1}) where {T <: GeneticSeq}
    if length(x) != 4
      @error "Invalid input, must have a length of 4"
    end
    return reshape(x, (4,1))
  end

  DNASeq(x) = DNASeq(_bitarray(DNASeq, x), checkinput=false)
  RNASeq(x) = RNASeq(_bitarray(RNASeq, x), checkinput=false)

  function setindex!(x::T, a, i) where {T <: GeneticSeq}
    x.data[:,i] = _bitarray(T, a)
    return x
  end

  function convert(::Type{String}, x::T) where {T <: GeneticSeq}
    s = ""
    for i in 1:length(x)
      if x.data[1,i]
        if x.data[2,i]
          if x.data[3,i]
            if x.data[4,i]
              # A C G T/U
              s *= "N"
            else
              # A C G
              s *= "V"
            end
          else
            if x.data[4,i]
              # A C T/U
              s *= "H"
            else
              # A C
              s *= "M"
            end
          end
        else
          if x.data[3,i]
            if x.data[4,i]
              # A G T/U
              s *= "D"
            else
              # A G
              s *= "R"
            end
          else
            if x.data[4,i]
              # A T/U
              s *= "W"
            else
              # A
              s *= "A"
            end
          end
        end
      else
        if x.data[2,i]
          if x.data[3,i]
            if x.data[4,i]
              # C G T/U
              s *= "B"
            else
              # C G
              s *= "S"
            end
          else
            if x.data[4,i]
              # C T/U
              s *= "Y"
            else
              # C
              s *= "C"
            end
          end
        else
          if x.data[3,i]
            if x.data[4,i]
              # G T/U
              s *= "K"
            else
              # G
              s *= "G"
            end
          else
            if x.data[4,i]
              if T == DNASeq
                # T
                s *= "T"
              else
                # U
                s *= "U"
              end
            else
              # gap
              s *= "-"
            end
          end
        end
      end
    end
    return s
  end

  function show(io::IO, x::T) where {T <: GeneticSeq}
    len = length(x)
    println(io, "$(len)nt $(_seq(T)) sequence")
    if len <= 26
      print(io, convert(String, x))
    else
      print(io, convert(String, x[1:13]) *
                "..." *
                convert(String, x[len-13:len]))
    end
  end

  function ==(x::T, y::T) where {T <: GeneticSeq}
    return x.data == y.data
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

  export GeneticSeq, DNASeq, RNASeq, Weights

end # module
