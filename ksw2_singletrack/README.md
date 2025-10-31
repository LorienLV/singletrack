# KSW2 Singletrack
This repository includes the implementation of the difference recurrence relations proposed by the Suzuki–Kasahara algorithm. Specifically, it uses KSW2, the implementation used by the widely adopted read-mapping tool minimap2.
The folder contains both the baseline version (without SingleTrack) and the SingleTrack-enhanced version of the algorithm. 
## Prerequisites 
To compile and use this application, you need a C/C++ compiler (e.g., gcc, clang).
KSW2 relies on SIMD instructions for computing alignments, and therefore, this implementation requires SIMD support.
You must have either SSE2/SSE4.2 (for x86 architectures) or NEON (for ARM architectures) enabled to compile and run the code.
## Instalation 
To compile the library, run the following commands:
```console
git clone https://github.com/criniguez/ksw2_singletrack.git
cd ksw2_singletrack
make clean all
```
## Benchmarking
We provide a simple alignment tool to align files containing pairs of sequences.
 Bear in mind that the Singletrack is not used by default, neither is used gap-affine nor dual gap-affine. Also if u use --only-score Singletrack is not used as is not needed.
### Usage
```console
Usage: benchmark_singletrack [OPTIONS]
Options:
  -h, --help            Show this help message
  -d, --dataset <file>  Dataset file
  -o, --output <file>   Output file
  --match <int>         The match score [default=0]
  --mismatch <int>      The mismatch score [default=1]
  --gapo <int>          The gap open penalty (Optional)
  --gape <int>          The gap extension penalty [default=1]
  --gapo2 <int>         The second gap open penalty for dual gap-affine (Optional)
  --gape2 <int>         The second gap extension penalty for dual gap-affine (Optional)
  --single-track        Apply SingleTrack strategy
  --only-score          Compute only alignment score (no CIGAR)
```
### Example
```console
./benchmark_singletrack -d sequences.seq --match 0 --mismatch 4 --gapo 6 --gape 2 -o cigar.out
Running AFFINE
Time Alignment: 2.26605s
```
The benchmark tool reports the execution time spent during the alignment for the selected algorithm.
