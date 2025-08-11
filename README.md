## PhIP-Seq Analysis with edgeR in R
This repository contains an R script for analyzing PhIP-Seq data using the `edgeR` package. It performs:
- Data preprocessing
- Normalization
- Dispersion estimation
- MDS visualization
- Differential expression analysis (exact test)
- Export of results and diagnostic plots

## Requirements

- R >= 4.0
- Packages:
  - edgeR
  - ggplot2
  - gtools
  ## Usage
  Be sure to define:
<pre>  data_file = "data/all_counts.tsv" 
  metadata_file = "metadata/metadata_file.txt" 
  group = "Group"
  baseLineGroup = "seronegative"
  n_replicat = 2
  number_of_peptide = 385968
  output = "output/"  </pre>
  ## Output
  - Normalized counts
  - MDS plot
  - Differential expression results
  - MD plots
  - Saved edgeR objects

