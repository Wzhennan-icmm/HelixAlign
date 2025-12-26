# HelixAlign

English | [中文](README.zh.md)

## Overview

HelixAlign is an efficient sequence alignment tool implemented in Rust, providing functionality similar to nucmer with support for multiple matching algorithms and output formats. This project aims to offer a high-performance, memory-friendly alternative for genomic sequence alignment analysis.

## Main Features

### 1. Sequence Alignment Algorithms

- **MUM (Maximal Unique Matches)**: Compute maximal matches that are unique in both sequences
- **MAM (Maximal Unique Reference Matches)**: Compute maximal matches that are unique in the reference sequence (default algorithm)
- **MEM (Maximal Exact Matches)**: Compute all maximal matches regardless of their uniqueness

### 2. Core Components

#### Suffix Array Implementation
- Implements sparse suffix arrays for efficient pattern matching
- Provides `find_matches` method that returns match position values
- Supports LCP (Longest Common Prefix) array calculation

#### Sequence Processing
- Implements `DnaSequence` struct for DNA sequence handling
- Provides reverse complement calculation functionality

#### Genomic Statistics
- Calculates genomic statistics such as N50/N90
- Supports multi-sequence statistics calculation
- Implements `Display` trait for easy output formatting

### 3. Nucmer Parameter Support

The project implements all nucmer parameters to ensure feature parity:

#### Basic Matching Parameters
- `-mum`: Compute maximal matches that are unique in both sequences
- `-mumreference`/`-mumcand`: Compute maximal matches that are unique in the reference sequence (default)
- `-maxmatch`: Compute all maximal matches regardless of their uniqueness
- `-l`/`--minmatch`: Set the minimum length of a single exact match (default: 20)

#### Clustering and Extension Parameters
- `-b`/`--breaklen`: Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default: 200)
- `-c`/`--mincluster`: Sets the minimum length of a cluster of matches (default: 65)
- `-D`/`--diagdiff`: Set the maximum diagonal difference between two adjacent anchors in a cluster (default: 5)
- `-d`/`--diagfactor`: Set the maximum diagonal difference between two adjacent anchors in a cluster as a differential fraction of the gap length (default: 0.12)
- `-g`/`--maxgap`: Set the maximum gap between two adjacent matches in a cluster (default: 90)
- `-L`/`--minalign`: Minimum length of an alignment (default: 0)

#### Processing Options
- `-noextend`: Do not perform cluster extension step
- `-nooptimize`: No alignment score optimization
- `-nosimplify`: Don't simplify alignments by removing shadowed clusters
- `-f`/`--forward`: Use only the forward strand of the Query sequences
- `-r`/`--reverse`: Use only the reverse complement of the Query sequences

#### Output and File Options
- `-p`/`--prefix`: Write output to PREFIX.delta (default: out)
- `--delta`: Output delta file to specified path
- `--sam-short`: Output SAM file, short format
- `--sam-long`: Output SAM file, long format
- `--save`: Save suffix array to files
- `--load`: Load suffix array from file

#### Advanced Options
- `-banded`: Enforce absolute banding of dynamic programming matrix based on diagdiff parameter
- `-large`: Force the use of large offsets
- `-G`/`--genome`: Map genome to genome (long query sequences)
- `-M`/`--max-chunk`: Set maximum chunk size
- `-t`/`--threads`: Set number of threads to use
- `-batch`: Proceed by batch of chunks from the reference
- `-format`: Specify output format (default, delta, paf, sam)
- `-stats`: Show reference and query sequence statistics (N50, N90, etc.)

### 4. Output Formats
- **Default**: Default format
- **Delta**: nucmer-compatible delta format
- **PAF**: Pairwise mApping Format
- **SAM**: Sequence Alignment/Map format

### 5. Multi-threading Support
- Parallel processing using Rayon
- Support for custom thread count
- Progress bar implementation

## Installation

### Compile from Source

```bash
git clone https://github.com/yourusername/helixalign.git
cd helixalign
cargo build --release
```

The compiled executable will be located at `target/release/helixalign`.

## Usage

### Basic Usage

```bash
helixalign reference.fa query.fa
```

### Using Specific Algorithm and Parameters

```bash
helixalign -maxmatch -l 20 -t 4 -format sam reference.fa query.fa
```

### Display Statistics

```bash
helixalign -stats reference.fa query.fa
```

### Complete Example with All Parameters

```bash
helixalign -maxmatch -b 200 -c 65 -D 5 -d 0.12 -noextend -f -g 90 -l 20 -L 0 -nooptimize -r -nosimplify -p results --delta results.delta --sam-short results.sam --batch 10000 -banded -large -G -M 50000 -t 8 -format paf -stats reference.fa query.fa
```

## Technical Highlights

### 1. Efficient Algorithm Implementation
- Sparse suffix arrays for efficient pattern matching
- Support for multiple matching algorithms (MUM, MAM, MEM)
- Implementation of maximal match finding and extension

### 2. Memory Optimization
- Sparse suffix arrays reduce memory usage
- Batch processing mechanism for large file handling
- Optional large offset support

### 3. User-Friendly
- Detailed command-line help information
- Progress bar for processing status
- Multiple output format options

## Project Structure

```
helixalign/
├── src/
│   ├── main.rs              # Main program entry and command-line parsing
│   ├── lib.rs               # Library file, exports all modules
│   ├── sequence.rs          # DNA sequence processing
│   ├── suffix_array.rs      # Suffix array implementation
│   ├── algorithms.rs        # Matching algorithm implementation
│   ├── nucmer.rs            # Nucmer alignment implementation
│   ├── genomic_stats.rs     # Genomic statistics calculation
│   └── output_format.rs     # Output format handling
├── Cargo.toml               # Project configuration and dependencies
├── README.md                # Project documentation
└── README.zh.md             # Chinese documentation
```

## Dependencies

- `rayon`: Parallel computing support
- `indicatif`: Progress bar display
- `clap`: Command-line argument parsing (planned)

## Performance

While maintaining feature parity with nucmer, HelixAlign provides the following performance advantages:

- Faster suffix array construction
- Lower memory usage
- Better multi-threading scalability

## Contributing

Bug reports and feature requests are welcome! If you'd like to contribute code, please:

1. Fork this repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This project is inspired by MUMmer4 and aims to provide similar functionality
- Thanks to the Rust community for excellent tools and libraries

## Contact

For questions or suggestions, please contact:

- Submit Issues: [GitHub Issues](https://github.com/yourusername/helixalign/issues)
- Email: your.email@example.com
