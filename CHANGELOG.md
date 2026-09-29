# Changelog

This history was reconstructed from the commits and source changes between
successive tags. Changes after the latest tag are listed as unreleased, even
where the package version has already been updated.

## Unreleased (version set to 0.4.0)

- Add optional explicit random number generators to weighted DNA and RNA
  sequence sampling, while keeping calls without an RNG available.
- Make string conversion and sequence construction faster for long sequences.
  Character lookup now uses immutable masks and tuples; StaticArrays is no
  longer required.
- Fix construction of single-nucleotide sequences so each has its own
  correctly shaped data.
- Allow StatsBase 0.34 and move testing to refreshed GitHub Actions workflows.

## [0.3.2](https://github.com/jangevaare/GeneticBitArrays.jl/compare/v0.3.1...v0.3.2)

- Allow StaticArrays 1 alongside 0.12.
- Refresh CI and release automation, add citation metadata, and remove the
  checked-in dependency manifest. No sequence behavior changed in this tag.

## [0.3.1](https://github.com/jangevaare/GeneticBitArrays.jl/compare/v0.3.0...v0.3.1)

- Allow StatsBase 0.33 and test against Julia 1.4.
- Add code coverage reporting, automated release tagging, and clearer README
  installation instructions. No sequence code changed in this tag.

## [0.3.0](https://github.com/jangevaare/GeneticBitArrays.jl/compare/v0.2.1...v0.3.0)

- Represent IUPAC ambiguity codes and gaps in DNA and RNA sequences, in
  addition to the four unambiguous nucleotides. String conversion and display
  now preserve those symbols.
- Consolidate DNA and RNA handling into one implementation and improve
  construction and lookup performance.
- Keep compatibility with Julia 1.0 and expand the sequence examples.

## [0.2.1](https://github.com/jangevaare/GeneticBitArrays.jl/compare/v0.2.0...v0.2.1)

- Generate a sequence from a vector of weights, allowing different nucleotide
  probabilities at each site.
- Export the common `GeneticSeq` type.

## [0.2.0](https://github.com/jangevaare/GeneticBitArrays.jl/compare/v0.1.0...v0.2.0)

- Add random DNA and RNA sequence generation from weighted nucleotide
  probabilities.
- Validate sequence input and sampling weights, with tests for invalid values.

## [0.1.0](https://github.com/jangevaare/GeneticBitArrays.jl/tree/v0.1.0)

- Initial release of bit-array-backed `DNASeq` and `RNASeq` types, with
  construction from nucleotide data, indexing, display, equality, and basic
  tests.
