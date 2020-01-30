# GeneticBitArrays.jl
[![Latest Release](https://img.shields.io/github/release/jangevaare/GeneticBitArrays.jl.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/releases/latest)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/jangevaare/GeneticBitArrays.jl.svg?branch=master)](https://travis-ci.org/jangevaare/GeneticBitArrays.jl)

Minimal representations of DNA and RNA genetic sequences using `BitArray`s in [Julia](https://julialang.org).

## Use
Input using `String`, `Vector{Char}` with nucleotides indicated by their [IUPAC code](https://www.bioinformatics.org/sms/iupac.html). You may also construct a sequence with a 4 x n `BitArray` - the same way sequences are represented internally with this package.

## Example
```
julia> x = DNASeq("NVHMDRWABSYCKGT-")

16nt DNA sequence
NVHMDRWABSYCKGT-

julia> x.data

4×15 BitArray{2}:
 1  1  1  1  1  1  1  1  0  0  0  0  0  0  0
 1  1  1  1  0  0  0  0  1  1  1  1  0  0  0
 1  1  0  0  1  1  0  0  1  1  0  0  1  1  0
 1  0  1  0  1  0  1  0  1  0  1  0  1  0  0
```

## Notes
* For a full featured package for biological sequences see [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl).
