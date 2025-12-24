# MUMmer-rs: Rust Implementation of MUMmer

This is a Rust implementation of MUMmer, a versatile alignment tool for DNA and protein sequences. The original MUMmer is written in C/C++ and is widely used for genome alignment tasks.

## Features

- **Suffix Array Construction**: Efficient construction of suffix arrays for fast pattern matching
- **Multiple Algorithms**: Support for MUM (Maximal Unique Match), MAM (Maximal Almost-Unique Match), and MEM (Maximal Exact Match) algorithms
- **Command-Line Interface**: Compatible with the original MUMmer command-line interface
- **FASTA Support**: Can process FASTA format sequence files

## Algorithms

- **MUM (Maximal Unique Match)**: Finds matches that are unique in both sequences
- **MAM (Maximal Almost-Unique Match)**: Finds matches that are unique in the reference sequence (default)
- **MEM (Maximal Exact Match)**: Finds all maximal matches regardless of uniqueness

## Usage

```bash
# Basic usage
cargo run -- [options] <reference-file> <query-file>

# Examples:
# Find matches with minimum length 20 (default)
cargo run -- reference.fa query.fa

# Find matches with minimum length 10
cargo run -- -l 10 reference.fa query.fa

# Use MUM algorithm
cargo run -- -mum -l 10 reference.fa query.fa

# Use MEM algorithm (finds all matches)
cargo run -- -maxmatch -l 10 reference.fa query.fa
```

## Options

- `-mum`: Compute maximal matches that are unique in both sequences
- `-mumreference`: Compute maximal matches that are unique in the reference sequence (default)
- `-maxmatch`: Compute all maximal matches regardless of their uniqueness
- `-l <n>`: Set the minimum length of a match (default: 20)

## Architecture

The implementation consists of several key modules:

- `sequence.rs`: DNA sequence handling and manipulation
- `suffix_array.rs`: Suffix array construction and search algorithms
- `algorithms.rs`: Core MUMmer algorithms (MUM, MAM, MEM)
- `main.rs`: Command-line interface

## Performance

The Rust implementation maintains the algorithmic efficiency of the original while providing memory safety and modern language features. The core algorithms use suffix arrays with LCP (Longest Common Prefix) arrays for efficient pattern matching.

## Testing

Run the tests with:
```bash
cargo test
```

## Building

To build the project:
```bash
cargo build --release
```

The resulting binary will be in `target/release/mummer-rs`.

## License

This project is licensed under the same terms as the original MUMmer software.