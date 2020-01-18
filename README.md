# GeneticBitArrays.jl
[![Build Status](https://travis-ci.org/jangevaare/GeneticBitArrays.jl.svg?branch=master)](https://travis-ci.org/jangevaare/GeneticBitArrays.jl)

Minimal representations of DNA and RNA genetic sequences.

## Use
Input using `String`, `Vector{Char}` with nucleotides indicated by their [IUPAC code](https://www.bioinformatics.org/sms/iupac.html). You may also construct a sequence with a 4 x n `BitArray`, which is how sequences are represented internally with this package.

## Example
```
julia> x = DNASeq("NVHMDRWABSYCKGT-")

15nt DNA sequence
NVHMDRWABSYCKG-

julia> x.data

4×15 BitArray{2}:
 1  1  1  1  1  1  1  1  0  0  0  0  0  0  0
 1  1  1  1  0  0  0  0  1  1  1  1  0  0  0
 1  1  0  0  1  1  0  0  1  1  0  0  1  1  0
 1  0  1  0  1  0  1  0  1  0  1  0  1  0  0
```

## Notes
* For a full featured package for genetic sequences see [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl).
