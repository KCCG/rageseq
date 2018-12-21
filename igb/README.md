<p align="center"><img src="images/RAGE_SEQ_LOGO" alt="RAGE-seq" width="25%"></p>

## RAGE-seq

Computational pipeline for the RAGE-seq workflow to accurately characterize gene-regulatory events of complex and heterogeneous cellular populations.

### Dependencies

RAGE-seq requires the following software dependencies (supported version in parentheses):

* cellranger (2.1.1)
* Albacore ()
*


### Overview

The pipeline is designed for execution on massively parallel computing infrastructure, with the following operations:

### General pipeline

1. Demultiplex sequencing data
2. _in silico_ and _in vitro_ capture
3. UMI clustering
4. Assembly
5. fast5_fetcher
6. Consensus polishing

### Immune receptor specific pipeline

* IgBLAST
* V(D)J indel correction
* IgBLAST
* Constant region calling


## Required inputs

### Cellranger

* BCLs
* HDF5 of molecule_info
* List of barcodes

### Nanopore

* FASTQs
* FAST5s
* Sequencing summary

## Example

```bash
# Cellranger

# ONT

# Pull out barcodes

# Filter with baits

# Fuzzy matching

# Assembly

# fast5_fetcher

# Polishing

# IgBLAST

# Indel correction

# IgBLAST

# V(D)J rearrangements



```
