# Singletrack

This repository contains implementations of the Singletrack algorithm presented [here](https://www.biorxiv.org/content/10.1101/2025.10.31.685625v1). Singletrack enables efficient backtracking in sequence alignment by utilizing only the scores from the substitution matrix (M), i.e., without access to the indel matrices (I1, D1, I2, D2). As a result, only the M matrix and a minimal scope of the indel matrices need to be stored during the alignment step, reducing memory consumption by 3x for gap-affine and 5x for dual gap-affine alignments and improving memory hierarchy utilization on modern hardware, often translating into improved performance.

The repository contains three different implementations of Singletrack, each with its own README:
- basic: A basic C++ implementation of the classical dynamic programming algorithm with and without using Singletrack.
- WFA2-lib: Singletrack integrated into WFA2-lib (commit 2ec28919af3a1b545acfb38e9bdefd160a87f266).
- KSW2: Singletrack integrated into KSW2 (commit 289609bd9e5381a13b16239d0a7703f1ff03f9ca).

Additionally, we include a tool to generate simulated datasets under `tools/python`.

## Datasets

The datasets used in the article are available in [Zenodo](https://zenodo.org/records/17525721) 

All aligners in this repository expect a dataset file containing one sequence per line, organized as pairs. The first sequence in each pair (target) starts with >, and the second (query) starts with <. For example:

```
>AA
<AG
>ACAC
<CTG
```

## Cite Us

> **Lorién López-Villellas, Cristian Iñiguez, Albert Jiménez-Blanco, Quim Aguado-Puig, Miquel Moretó, Jesús Alastruey-Benedé, Pablo Ibáñez, Santiago Marco-Sola.** *Singletrack: An Algorithm for Improving Memory Consumption and Performance of Gap-Affine Sequence Alignment*, bioRXiv, 2025.

```
@article{LpezVillellas2025,
  title = {Singletrack: An Algorithm for Improving Memory Consumption and Performance of Gap-Affine Sequence Alignment},
  url = {http://dx.doi.org/10.1101/2025.10.31.685625},
  DOI = {10.1101/2025.10.31.685625},
  publisher = {Cold Spring Harbor Laboratory},
  author = {López-Villellas,  Lorién and Iñiguez,  Cristian and Jiménez-Blanco,  Albert and Aguado-Puig,  Quim and Moretó,  Miquel and Alastruey-Benedé,  Jesús and Ibáñez,  Pablo and Marco-Sola,  Santiago},
  year = {2025},
  month = nov 
}
```
