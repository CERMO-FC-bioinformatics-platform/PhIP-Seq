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
<pre>  
  # Path to the raw read count table (must include a column named "ID" for peptide identifiers)
  data_file = "data/all_counts.tsv" 
  
  # Path to the metadata file describing the samples (e.g., group assignment, conditions)
  metadata_file = "metadata/metadata_file.txt" 
  
  # Column name in metadata that indicates the experimental groups
  group = "Group"
  
  # The baseline/control group name in the metadata, used for comparisons
  baseLineGroup = "seronegative"
  
  # Number of replicates per group
  n_replicat = 2
  
  # Total number of peptides (rows) in the count table
  number_of_peptide = 385968
  
  # Output directory where plots, results, and R objects will be saved
  output = "output/"  </pre>
  ## Output
  - Normalized counts
  - MDS plot
  - Differential expression results
  - MD plots
  - Saved edgeR objects

