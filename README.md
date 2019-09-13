# GeneticBitArrays.jl
[![Build Status](https://travis-ci.com/jangevaare/GeneticBitArrays.jl.svg?branch=master)](https://travis-ci.com/jangevaare/GeneticBitArrays.jl)

Minimal representations of genetic sequences.

## Use
Input using `String`, `Vector{Char}`, or use an appropriate `BitArray`.

*e.g.*
```
julia> DNASeq("AAAAGCCT")

8 nt DNA sequence
AAAAGCCT
```

## Representation
A 4 x n `BitArray` is used to represent sequences internally, where n is sequence length.


```
julia> x = DNASeq("AAAAGCCT")

julia> x.data

4×8 BitArray{2}:
 1  1  1  1  0  0  0  0
 0  0  0  0  0  1  1  0
 0  0  0  0  1  0  0  0
 0  0  0  0  0  0  0  1
```

## Notes
* Ambiguities are not currently supported.
* For a full featured package for genetic sequences see [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl).
