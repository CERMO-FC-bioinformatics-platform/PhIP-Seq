# Load necessary packages
library(edgeR)      # For differential expression analysis
library(ggplot2)    # For plotting
library(gtools)     # For sorting (e.g., mixedsort)

# Main function to run edgeR analysis on PhIP-Seq data
run_edgeR_analysis = function(data_file, metadata_file, out_dir, baseLineGroup, group, n_replicat, number_of_peptide) {
  
  # Define and create output directories
  outPath_Robj = file.path(out_dir, "R_objects")
  outPath_plots = file.path(out_dir, "plots")
  outPath_tables = file.path(out_dir, "tables")
  
  dir.create(outPath_Robj, showWarnings = TRUE, recursive = TRUE)
  dir.create(outPath_plots, showWarnings = TRUE, recursive = TRUE)
  dir.create(outPath_tables, showWarnings = TRUE, recursive = TRUE)
  
  # Load the count table
  count_table = read.table(data_file, sep = "\t", header = TRUE)
  
  # Extract numeric columns (sample counts) and the "ID" column
  numeric_cols = sapply(count_table, is.numeric)
  numeric_and_id_cols = c("ID", names(count_table)[numeric_cols])
  count_table_new = count_table[numeric_and_id_cols]
  
  # Remove "INDEX" column if present
  if ("INDEX" %in% colnames(count_table_new)) {
    count_table_new = count_table_new[, -which(colnames(count_table_new) == "INDEX")]
  }
  
  # Store annotation (non-numeric) columns
  count_table_new_cols = colnames(count_table_new)
  annotation_df = count_table[c("ID", setdiff(colnames(count_table), count_table_new_cols))]
  
  # Simplify column names by removing suffixes after underscores
  new_colnames = sub("_.*$", "", colnames(count_table_new))
  colnames(count_table_new) = new_colnames
  
  # Convert count table to a matrix and set rownames
  count_table_new_mat = as.matrix(count_table_new[, -1])
  rownames(count_table_new_mat) = count_table_new[, 1]
  
  # Sort column names in alphanumeric order
  sorted_names = mixedsort(colnames(count_table_new_mat))
  count_table_new_mat = count_table_new_mat[, sorted_names]
  
  # Load metadata (sample info)
  sample_metaData = read.table(metadata_file, sep = "\t", header = TRUE, row.names = 1)
  
  # Check that count matrix columns match metadata row names
  if (!identical(colnames(count_table_new_mat), rownames(sample_metaData))) {
    stop("The order of columns in the count matrix does not match the rows in the metadata file.")
  } else {
    cat("The order of columns matches between count matrix and metadata.\n")
  }
  
  # Create DGEList object
  edgeR_obj = DGEList(counts = count_table_new_mat, 
                      group = sample_metaData[[group]], 
                      genes = annotation_df)
  
  # Get group levels
  level_group = levels(edgeR_obj$samples$group)
  
  # Plot library sizes
  plot_name = file.path(outPath_plots, "1_Librarysize.png")
  png(plot_name, width = 25, height = 8, units = "in", res = 300)
  barplot(edgeR_obj$samples$lib.size * 1e-6, 
          names = rownames(edgeR_obj$samples), 
          las = 2, 
          ylab = "Library size (millions)", 
          col = "aquamarine3", 
          cex.axis = 0.8, 
          cex.names = 0.90)
  dev.off()
  
  # Normalize counts using TMM (edgeR does not transform counts)
  edgeR_obj2 = calcNormFactors(edgeR_obj)
  
  # Save normalized counts (CPM)
  cpm_table = cpm(edgeR_obj2, log = FALSE)
  cpm_file = file.path(outPath_tables, "edgeR_normalized_count.tsv")
  write.table(cpm_table, cpm_file, sep = "\t", quote = FALSE, row.names = TRUE)
  
  # Get number of groups
  n_group = length(level_group)
  
  # Plot MDS (multidimensional scaling) to visualize sample distances
  plot_name = file.path(outPath_plots, "2_plotMDS_afterNorm.png")
  png(plot_name, width = 25, height = 10, units = "in", res = 300)
  pltMDS = plotMDS(edgeR_obj2, col = rep(1:n_group, each = n_replicat))
  print(pltMDS)
  dev.off()
  
  # Save MDS coordinates to file
  mds_df = data.frame(Sample = rownames(pltMDS$distance.matrix.squared), 
                      X = pltMDS$x, 
                      Y = pltMDS$y)
  mds_file = file.path(outPath_tables, "MDS_XYcoordinates.tsv")
  write.table(mds_df, file = mds_file, row.names = FALSE, sep = "\t")
  
  # Prepare design matrix for DE analysis
  groups = as.factor(sample_metaData[[group]])
  design = model.matrix(~0 + groups)
  colnames(design) = levels(groups)
  
  # Reorder columns so that the baseline group is last
  baseline_col = design[, baseLineGroup]
  design = design[, !colnames(design) %in% baseLineGroup]
  design = cbind(design, baseline_col)
  colnames(design)[colnames(design) == "baseline_col"] = baseLineGroup
  
  # Estimate dispersion
  edgeR_obj3 = estimateDisp(edgeR_obj2, design, robust = TRUE)
  
  # Calculate biological coefficient of variation (BCV)
  bcv = sqrt(edgeR_obj3$common.dispersion)
  rounded_bcv = round(bcv, 2)
  print(paste0("BCV: ", rounded_bcv))
  
  # Save edgeR object
  edgeR_file = file.path(outPath_Robj, "edgeR_obj3.rds")
  saveRDS(edgeR_obj3, file = edgeR_file)
  
  # Differential expression analysis: compare each group to baseline
  pair_names = colnames(design)[colnames(design) != baseLineGroup]
  
  for (condition in pair_names) {
    # Define comparison pair
    comparison_pair = c(baseLineGroup, condition)
    
    # Run exact test
    res_edgeR = exactTest(edgeR_obj3, pair = comparison_pair)
    
    # Save MD plot
    plot_name = file.path(outPath_plots, paste0("3_plotMD_", condition, "_vs_", baseLineGroup, ".png"))
    png(plot_name, width = 8, height = 5, units = "in", res = 300)
    plotMD(res_edgeR)
    abline(h = c(-1, 1), col = "blue")
    dev.off()
    
    # Save top results
    topTags_table = topTags(res_edgeR, adjust.method = "fdr", n = number_of_peptide)
    topTags_file = file.path(outPath_tables, paste0("FDR_edgeR_DE_res_", condition, "_vs_", baseLineGroup, ".tsv"))
    write.table(topTags_table, topTags_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# =============================
# Example usage
# =============================

# Define input and output paths and parameters

# Path to raw count table (must include 'ID' column)
data_file = "../data/all_counts.tsv"

# Path to metadata file (must match sample order in count table)
metadata_file = "../metadata/metadata_Sep2023_F9F10_neg.txt"

# Column in metadata indicating groups
group = "Group"

# Baseline group name (e.g., "seronegative")
baseLineGroup = "seronegative"

# Total number of peptides (rows) in count table
number_of_peptide = 385968

# Number of replicates per group
n_replicat = 2

# Output directory
output = "../output_F9F10_test"

# Run the analysis
run_edgeR_analysis(data_file, metadata_file, output, baseLineGroup, group, n_replicat, number_of_peptide)
