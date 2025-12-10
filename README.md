# [Fuse: Multiple Network Alignment via Data Fusion](https://doi.org/10.1093/bioinformatics/btv731)

A novel global multiple network aligner that fuses sequence similarities and network wiring patterns using non-negative matrix tri-factorization (NMTF).

## Authors

- Vladimir Gligorijevic
- Noel Malod-Dognin
- Natasa Przulj

**Corresponding Author:** Prof. Natasa Przulj  
📧 natasa.przulj@mbzuai.ac.ae

## Overview

Fuse is a novel global multiple network aligner that operates in two steps:

1. **Data Fusion**: Computes novel similarity scores between proteins by fusing sequence similarities and network wiring patterns over all proteins in all PPI networks being aligned, using non-negative matrix tri-factorization (NMTF)
2. **Alignment**: Constructs a one-to-one global multiple network alignment using an approximate maximum weight k-partite matching solver

## Resources

- **[Supplementary Material](http://www0.cs.ucl.ac.uk/staff/natasa/FUSE/FUSE_supplement.pdf)** (PDF)
- **[Dataset for Fuse](http://www0.cs.ucl.ac.uk/staff/natasa/FUSE/FUSE_dataset.7zt)**
- **[Predicted Clusters Table](https://github.com/przulj-lab/fuse/blob/main/predicted_clusters_official_gene_symbol.xls)** (XLS)

*Note: If download links don't work when clicked, copy and paste the URLs directly into your browser's address bar.*

## Requirements

### Data-Fusion Step
- MATLAB R2013a or higher
- List of all PPI network files in edgelist format
- Input file containing sequence similarity between proteins (e.g., 1-eval)
- Number of NMTF iterations
- List of rank parameters
- Gamma parameter

### Alignment Step
- PPI network files in LEDA format
- Network list file containing full path to each network (one per line)
- File containing sequence similarity between proteins
- File containing NMTF-based similarity between proteins

## Installation

```bash
git clone https://github.com/przulj-lab/fuse.git
cd fuse
```

## Usage

### Step 1: Data-Fusion

Run in MATLAB console:

```matlab
run_nmtf({'./Nets/CElegans.edgelist', './Nets/DMelanogaster.edgelist', './Nets/HSapiens.edgelist', './Nets/MMusculus.edgelist', './Nets/SCerevisiae.edgelist'}, 'sequence_scores.lst', [80, 90, 80, 70, 50], 1000, 0.01, 'nmft_scores.lst')
```

**Parameters:**
- List of PPI network file paths
- Sequence similarity file
- Rank parameters array
- Number of iterations (1000)
- Gamma parameter (0.01)
- Output file name

This step produces `nmft_scores.lst`, which serves as input for the alignment step.

### Step 2: Alignment

Example network list file (`network_list.txt`):
```
./Nets/CElegans.gw
./Nets/DMelanogaster.gw
./Nets/HSapiens.gw
./Nets/MMusculus.gw
./Nets/SCerevisiae.gw
```

Run the alignment:

```bash
./FUSE.exe -n network_list.txt -s sequence_scores.lst -t nmft_scores.lst -o output_alignment
```

**Parameters:**
- `-n` : Network list file
- `-s` : Sequence similarity scores file
- `-t` : NMTF-based similarity scores file
- `-o` : Output alignment file

For complete list of parameters:
```bash
./FUSE.exe -h
```

## Other Multiple Network Alignment Algorithms

- **[NetCoffee](https://code.google.com/archive/p/netcoffee/)** 
- **[Smetana](https://biomlsp.com/)**
- **[IsoRankN](https://cb.csail.mit.edu//mna/)**

## Citation

If you use this code or data in your research, please cite:

```
Vladimir Gligorijević, Noël Malod-Dognin, Nataša Pržulj,
Fuse: multiple network alignment via data fusion,
Bioinformatics, Volume 32, Issue 8, April 2016,
Pages 1195–1203, https://doi.org/10.1093/bioinformatics/btv731
```

## Contact

For questions or issues, please contact Prof. Natasa Przulj at natasa.przulj@mbzuai.ac.ae
