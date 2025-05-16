# Singletrack

This repository contains an example implementation of the Singletrack technique for gap-affine and dual gap-affine sequence alignment algorithms presented [here](TODO). This technique enables efficient backtracking in sequence alignment by utilizing only the scores from the substitution matrix (M), i.e., without access to the indel matrices (I1, D1, I2, D2). As a result, only the M matrix and a minimal scope of the indel matrices need to be stored during the alignment step, reducing memory consumption by 3x for gap-affine and 5x for dual gap-affine alignments and improving memory hierarchy utilization on modern hardware, often translating into improved performance.

The repository includes two dynamic programming-based sequence aligners: one implementing the traditional backtracking approach (`src/dp_aligner_base.*`), and another utilizing the Singletrack technique (`src/dp_aligner_singletrack.*`). The Singletrack-based aligner incorporates further optimizations in the alignment phase to fully exploit the benefits of the technique.

## Getting Started

### Prerequisites

```
CMake (for building the project)
A C++ compiler with full C++23 support (e.g., GCC 12 or later)
Python 3 (for synthetic dataset generation)
```

### Installation

To build the project, execute the following commands:

```
mkdir build
cd build
cmake ..
make
```

If you need to specify a particular C++23-compatible compiler (for example, GCC 13), you can set it explicitly when running CMake:

```
cmake -DCMAKE_CXX_COMPILER=g++-13 ..
```

### Generating Example Datasets

A utility script is provided to generate synthetic DNA sequence datasets:

```
python ./tools/generate_dataset.py -n NUMBER_OF_SEQUENCE_PAIRS -l LENGTH_OF_THE_SEQUENCES
```

Replace NUMBER_OF_SEQUENCE_PAIRS and LENGTH_OF_THE_SEQUENCES with your desired values.

### Benchmarking

A benchmarking tool is included to evaluate and compare the traditional and Singletrack-based aligners in terms of execution time and memory usage:

```
./build/benchmark --help

Usage: benchmark [OPTIONS]
Options:
  -h, --help            Show this help message
  -v, --verbose         Verbose output, i.e., print alignments and scores
  -d, --dataset <file>  Dataset file
  --match <int>         The match score [default=0]
  --mismatch <int>      The mismatch score [default=1]
  --gapo <int>          The gap open penalty (Optional)
  --gape <int>          The gap extension penalty [default=1]
  --gapo2 <int>         The second gap open penalty for dual gap-affine (Optional)
  --gape2 <int>         The second gap extension penalty for dual gap-affine (Optional)
```

The benchmark tool reports execution time and memory consumption for both aligners, allowing direct performance comparison:

```
./build/benchmark --gapo 2 --gapo2 4 --gape2 6 -d ./dataset.txt

*** Using dual gap-affine ***
  Match: 0
  Mismatch: 1
  Gap open 1: 2
  Gap extension 1: 1
  Gap open 2: 4
  Gap extension 2: 6

Traditional DP Time: 1.01182 s
Singletrack DP Time: 0.613686 s

Traditional DP memory usage: 19.685 MB
Singletrack DP memory usage: 3.94464 MB
```

### Dataset Format

The dataset file should contain one sequence per line, organized as pairs. The first sequence in each pair (target) starts with >, and the second (query) starts with <. For example:

```
>AA
<AG
>ACAC
<CTG
```

## Cite Us

TODO