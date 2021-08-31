# GeneticBitArrays.jl
[![DOI](https://zenodo.org/badge/207670056.svg)](https://zenodo.org/badge/latestdoi/207670056)
[![Latest Release](https://img.shields.io/github/release/jangevaare/GeneticBitArrays.jl.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/releases/latest)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/blob/master/LICENSE)

[![test-lts](https://github.com/jangevaare/GeneticBitArrays.jl/actions/workflows/test-lts.yml/badge.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/actions/workflows/test-lts.yml)
[![test-stable](https://github.com/jangevaare/GeneticBitArrays.jl/actions/workflows/test-stable.yml/badge.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/actions/workflows/test-stable.yml)
[![test-nightly](https://github.com/jangevaare/GeneticBitArrays.jl/actions/workflows/test-nightly.yml/badge.svg)](https://github.com/jangevaare/GeneticBitArrays.jl/actions/workflows/test-nightly.yml)
[![codecov.io](http://codecov.io/github/jangevaare/GeneticBitArrays.jl/coverage.svg?branch=master)](http://codecov.io/github/jangevaare/GeneticBitArrays.jl?branch=master)

## Description
Minimal representations of DNA and RNA genetic sequences using `BitArray`s in [Julia](https://julialang.org).

## Installation
The current release can be installed from the Julia REPL with:

```julia
pkg> add GeneticBitArrays
```

The development version (master branch) can be installed with:

```julia
pkg> add GeneticBitArrays#master
```

## Use
Input using `String`, `Vector{Char}` with nucleotides indicated by their [IUPAC code](https://www.bioinformatics.org/sms/iupac.html). You may also construct a sequence with a 4 x n `BitArray` - the same way sequences are represented internally with this package.

## Example
```
julia> x = DNASeq("NVHMDRWABSYCKGT-")

16nt DNA sequence
NVHMDRWABSYCKGT-

julia> x.data

4×16 BitArray{2}:
 1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0
 1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0
 1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0
 1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0
```

## Notes
* For a full featured package for biological sequences see [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl).
