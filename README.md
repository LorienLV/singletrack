# Singletrack

This repository contains implementations of the Singletrack algorithm presented [here](TODO). Singletrack enables efficient backtracking in sequence alignment by utilizing only the scores from the substitution matrix (M), i.e., without access to the indel matrices (I1, D1, I2, D2). As a result, only the M matrix and a minimal scope of the indel matrices need to be stored during the alignment step, reducing memory consumption by 3x for gap-affine and 5x for dual gap-affine alignments and improving memory hierarchy utilization on modern hardware, often translating into improved performance.

The repository contains three different implementations of Singletrack, each with its own README:
- basic: A basic C++ implementation of the classical dynamic programming algorithm with and without using Singletrack.
- WFA2-lib: Singletrack integrated into WFA2-lib (commit 2ec28919af3a1b545acfb38e9bdefd160a87f266).
- KSW2: TODO.

Additionally, we include a tool to generate simulated datasets under `tools/python`.

## Dataset Format

All aligners in this repository expect a dataset file containing one sequence per line, organized as pairs. The first sequence in each pair (target) starts with >, and the second (query) starts with <. For example:

```
>AA
<AG
>ACAC
<CTG
```

## Cite Us

TODO