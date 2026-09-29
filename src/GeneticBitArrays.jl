module GeneticBitArrays

  import Random

  import Base.convert,
         Base.show,
         Base.length,
         Base.getindex,
         Base.setindex!,
         Base.==,
         Random.rand,
         Random.AbstractRNG,
         StatsBase.Weights,
         StatsBase.sample

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

  const _dnacharlookup = ('N', 'V', 'H', 'M', 'D', 'R', 'W', 'A',
                          'B', 'S', 'Y', 'C', 'K', 'G', 'T', '-')
  const _rnacharlookup = ('N', 'V', 'H', 'M', 'D', 'R', 'W', 'A',
                          'B', 'S', 'Y', 'C', 'K', 'G', 'U', '-')

  # Four-bit masks encode membership in A, C, G, T/U, in that order.
  # 0xff distinguishes invalid input from the valid gap mask (zero).
  function _masklookup(chars)
    masks = fill(UInt8(0xff), 128)
    for (i, c) in enumerate(chars)
      masks[Int(c) + 1] = UInt8(16 - i)
    end
    return Tuple(masks)
  end

  const _dnamasklookup = _masklookup(_dnacharlookup)
  const _rnamasklookup = _masklookup(_rnacharlookup)

  _lookup(::Type{DNASeq}) = _dnacharlookup
  _lookup(::Type{RNASeq}) = _rnacharlookup
  _masklookup(::Type{DNASeq}) = _dnamasklookup
  _masklookup(::Type{RNASeq}) = _rnamasklookup
  _seq(::Type{DNASeq}) = "DNA"
  _seq(::Type{RNASeq}) = "RNA"

  function _mask(::Type{T}, c::Char) where {T <: GeneticSeq}
    code = UInt32(c)
    return code < 128 ? _masklookup(T)[Int(code) + 1] : UInt8(0xff)
  end

  function _setmask!(bits::BitArray{2}, mask::UInt8, i::Int)
    bits[1, i] = (mask & 0x08) != 0
    bits[2, i] = (mask & 0x04) != 0
    bits[3, i] = (mask & 0x02) != 0
    bits[4, i] = (mask & 0x01) != 0
    return bits
  end

  function _bitarray(::Type{T}, x::Union{String, Vector{Char}}) where {T <: GeneticSeq}
    bits = BitArray{2}(undef, (4, length(x)))
    for (i, c) in enumerate(x)
      mask = _mask(T, c)
      if mask == 0xff
        throw(ErrorException("Unrecognized $(_seq(T)) `Char` $c at index $i"))
      end
      _setmask!(bits, mask, i)
    end
    return bits
  end

  function _bitarray(::Type{T}, x::Char) where {T <: GeneticSeq}
    mask = _mask(T, x)
    if mask == 0xff
      throw(ErrorException("Unrecognized $(_seq(T)) `Char` $x"))
    end
    return _setmask!(BitArray{2}(undef, (4, 1)), mask, 1)
  end

  function _bitarray(::Type{T}, x::BitArray{1}) where {T <: GeneticSeq}
    if length(x) != 4
      throw(ErrorException("Invalid input, must have a length of 4"))
    end
    return reshape(x, (4,1))
  end

  DNASeq(x) = DNASeq(_bitarray(DNASeq, x), checkinput=false)
  RNASeq(x) = RNASeq(_bitarray(RNASeq, x), checkinput=false)

  function setindex!(x::T, a, i) where {T <: GeneticSeq}
    x.data[:,i] = _bitarray(T, a)
    return x
  end

  """
      String(sequence::GeneticSeq)
      convert(String, sequence::GeneticSeq)

  Return the sequence as an IUPAC nucleotide string, including ambiguity codes
  and `-` for gaps.
  """
  function Base.String(x::T) where {T <: GeneticSeq}
    chars = _lookup(T)
    bytes = Vector{UInt8}(undef, length(x))
    for i in 1:length(x)
      # The lookup lists the four-bit patterns from 1111 through 0000.
      mask = 8 * x.data[1,i] + 4 * x.data[2,i] +
             2 * x.data[3,i] + x.data[4,i]
      bytes[i] = UInt8(chars[16 - mask])
    end
    return String(bytes)
  end

  convert(::Type{String}, x::GeneticSeq) = String(x)

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

  @static if isdefined(Random, :default_rng)
    _default_rng() = Random.default_rng()
  else
    _default_rng() = Random.GLOBAL_RNG
  end

  """
      rand([rng::AbstractRNG], T, weights::Weights, n::Integer; checkinput=true)
      rand([rng::AbstractRNG], T, weights::Vector{<:Weights}; checkinput=true)

  Generate a `DNASeq` or `RNASeq` with weights for A, C, G, and T/U.
  Supply one set of weights for all `n` sites, or one set per site.
  An explicit RNG controls every draw; omitted RNGs use Julia's default RNG.
  Exact random streams may change across Julia or dependency versions.

  # Example
  ```julia
  using GeneticBitArrays, Random
  rand(MersenneTwister(42), DNASeq, Weights([1, 1, 1, 1]), 100)
  ```
  """
  function rand(rng::AbstractRNG, ::Type{T}, w::W, n::Integer; checkinput::Bool=true) where {T <: GeneticSeq, W <: Weights}
    if checkinput && length(w) != 4
      throw(ErrorException("Invalid sampling weights for $T generation"))
    end
    x = falses(4, n)
    for i in 1:n
      x[sample(rng, 1:4, w), i] = true
    end
    return T(x, checkinput=false)
  end

  function rand(rng::AbstractRNG, ::Type{T}, wv::Vector{W}; checkinput::Bool=true) where {T <: GeneticSeq, W <: Weights}
    len = length(wv)
    x = falses(4, len)
    for i in 1:len
      if checkinput && length(wv[i]) != 4
        throw(ErrorException("Invalid sampling weights for $(i)th nt in $T generation"))
      end
      x[sample(rng, 1:4, wv[i]), i] = true
    end
    return T(x, checkinput=false)
  end

  rand(::Type{T}, w::Weights, n::Integer; checkinput::Bool=true) where {T <: GeneticSeq} =
    rand(_default_rng(), T, w, n; checkinput=checkinput)

  rand(::Type{T}, wv::Vector{<:Weights}; checkinput::Bool=true) where {T <: GeneticSeq} =
    rand(_default_rng(), T, wv; checkinput=checkinput)

  export GeneticSeq, DNASeq, RNASeq, Weights

end # module
