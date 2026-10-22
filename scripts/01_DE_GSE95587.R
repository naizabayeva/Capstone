# 01_DE_GSE95587.R
#
# Differential expression (DE) analysis for GSE95587.
# Comparing Alzheimer's disease (AD) vs Control.
#
# This is the third step for GSE95587:
# 1. Load the two datasets prepared in preprocessing
# 2. Compare my gene list with Shah et al.'s results
# 3. Run limma-voom on both filtering branches
# 4. Identify differentially expressed genes (DEGs)
# 5. Compare my results with Shah et al.
# 6. Save tables, plots, and a short summary
#
# The two branches are:
# shah     = filtering method described by Shah et al.
# standard = edgeR's filterByExpr() method
#
# The published DE results are used ONLY for comparison.
# They are not used to fit my models.


# 1. Packages and file paths ----------------------------------------

library(edgeR)
library(limma)

acc <- "GSE95587"

out_dir <- file.path("results", acc)
proc_dir <- file.path("data", "processed", acc)
ref_dir <- file.path("data", "reference")

shah_file <- file.path(
  ref_dir, "shah_GSE95587_full_DE_results.csv"
)

# Create output folders if they don't exist.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)


# 2. Load the processed data ----------------------------------------

# The .rds file was created in the preprocessing script.
# It should contain two DGEList objects:
# branches$shah and branches$standard.

branches <- readRDS(
  file.path(proc_dir, "dge_branches.rds")
)

# Check that the list contains the correct objects.
stopifnot(all(c("shah", "standard") %in% names(branches)))

stopifnot(
  inherits(branches$shah, "DGEList"),
  inherits(branches$standard, "DGEList")
)

cat("Loaded processed datasets:\n")
cat("Shah-style:", dim(branches$shah), "\n")
cat("Standard:", dim(branches$standard), "\n")


# 3. Check normalization factors ------------------------------------

# TMM normalization was already done during preprocessing.
# I don't need to repeat it here.
#
# I only check that the saved normalization factors exist
# and contain valid values.
#
# This check does not prove how TMM was calculated.
# The preprocessing script documents that decision.

for (b in names(branches)) {
  
  factors <- branches[[b]]$samples$norm.factors
  
  cat("\nBranch:", b, "\n")
  cat("Samples:", ncol(branches[[b]]), "\n")
  cat("TMM factors:", length(factors), "\n")
  
  if (length(factors) != ncol(branches[[b]]) ||
      any(!is.finite(factors)) ||
      any(factors <= 0)) {
    stop("Invalid normalization factors in branch: ", b)
  }
}


# 4. Load Shah et al.'s published results ---------------------------

# Download the reference DE results only if the file is missing.
# I am not using their code to perform the DE analysis.

if (!file.exists(shah_file)) {
  
  download.file(
    paste0(
      "https://raw.githubusercontent.com/",
      "Code-krafter22/meta-analysis_AD_biomarker/",
      "main/results/differential_expression/",
      "GSE95587/full_DE_results.csv"
    ),
    shah_file,
    mode = "wb"
  )
}

# Check that the downloaded file exists.
if (!file.exists(shah_file)) {
  stop("Shah reference file was not found.")
}

# Load the results.
# row.names = 1 assumes the first column contains gene names.
# This needs to match the actual reference CSV format.

shah <- read.csv(
  shah_file,
  row.names = 1,
  check.names = FALSE
)

# These are the columns needed for comparison.
required_cols <- c("logFC", "P.Value", "adj.P.Val")

if (!all(required_cols %in% names(shah))) {
  stop(
    "The Shah reference file does not have the expected columns. ",
    "Check names(shah) and the downloaded CSV."
  )
}

# Check that the relevant values are numeric.
stopifnot(
  is.numeric(shah$logFC),
  is.numeric(shah$P.Value),
  is.numeric(shah$adj.P.Val)
)

cat("\nShah reference file loaded.\n")
cat("Genes in reference:", nrow(shah), "\n")


# 5. Compare gene lists ---------------------------------------------

# My Shah-style filter kept 15,904 genes.
# Shah's results file also reportedly contains 15,904 genes.
#
# But having the same number is not enough.
# I also need to check whether they are the SAME genes.

my_genes <- rownames(branches$shah)
shah_genes <- rownames(shah)

common_genes <- intersect(my_genes, shah_genes)

cat("\nGene list comparison:\n")
cat("Genes in my Shah-style branch:", length(my_genes), "\n")
cat("Genes in Shah's results:", length(shah_genes), "\n")
cat("Genes in both:", length(common_genes), "\n")
cat("Only in mine:",
    length(setdiff(my_genes, shah_genes)), "\n")
cat("Only in Shah's:",
    length(setdiff(shah_genes, my_genes)), "\n")

# Save the genes that are different between the lists.
write.csv(
  data.frame(gene = setdiff(my_genes, shah_genes)),
  file.path(out_dir, "genes_only_in_mine.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene = setdiff(shah_genes, my_genes)),
  file.path(out_dir, "genes_only_in_shah.csv"),
  row.names = FALSE
)

if (length(common_genes) == 0) {
  stop("No matching gene names. Check the reference gene IDs.")
}


# 6. Run limma-voom -------------------------------------------------

# limma-voom is used for differential expression analysis.
#
# The steps are:
# design = describe the groups being compared
# voom   = calculate log2 CPM and precision weights
# lmFit  = fit a linear model for each gene
# eBayes = improve variance estimates using information
#          across all genes
#
# Control is the reference group.
# Positive logFC means higher expression in AD.
#
# This model only includes diagnosis, with no other covariates.

run_de <- function(dge, plot_file) {
  
  # Check the sample groups.
  group <- factor(
    dge$samples$group,
    levels = c("Control", "AD")
  )
  
  if (any(is.na(group))) {
    stop("Missing or unexpected sample group labels.")
  }
  
  # Create the design matrix.
  design <- model.matrix(~ group)
  
  # voom also makes a mean-variance plot.
  png(plot_file, width = 800, height = 700)
  
  v <- voom(dge, design, plot = TRUE)
  
  dev.off()
  
  # Fit a model for each gene.
  fit <- lmFit(v, design)
  
  # Apply empirical Bayes moderation.
  fit <- eBayes(fit)
  
  # Return all genes sorted by P-value.
  tt <- topTable(
    fit,
    coef = "groupAD",
    number = Inf,
    sort.by = "P"
  )
  
  return(tt)
}


# Run the same analysis on both filtering branches.

results <- list(
  shah = run_de(
    branches$shah,
    file.path(out_dir, "voom_shah.png")
  ),
  
  standard = run_de(
    branches$standard,
    file.path(out_dir, "voom_standard.png")
  )
)

# Save the complete DE results for both branches.
for (b in names(results)) {
  
  write.csv(
    results[[b]],
    file.path(out_dir, paste0("DE_", b, "_full.csv"))
  )
}

cat("\nDifferential expression analysis complete.\n")


# 7. Identify differentially expressed genes -----------------------

# DEG = differentially expressed gene.
#
# Here I use the threshold stated in my replication plan:
# unadjusted P < 0.05 and absolute log2FC > 0.58.
#
# I also report the number of genes with FDR < 0.05.
# FDR corrects for testing many genes at once.
#
# IMPORTANT: These are two different significance criteria.

is_deg <- function(tt) {
  
  tt$P.Value < 0.05 & abs(tt$logFC) > 0.58
  
}

shah_deg <- rownames(shah)[is_deg(shah)]

cat("\nDEGs in Shah's reference file:",
    length(shah_deg), "\n")


# 8. Compare my results with Shah et al. ----------------------------

# Compare the genes tested, DEG counts, fold changes,
# directions, and overlap with Shah's results.
#
# Spearman correlation tells me whether genes have similar
# fold-change rankings in my results and Shah's results.
#
# Top-500 overlap tells me how many of the 500 genes with
# the smallest P-values are shared.

compare_results <- function(tt, branch_name) {
  
  my_deg <- rownames(tt)[is_deg(tt)]
  
  common <- intersect(rownames(tt), rownames(shah))
  
  if (length(common) < 2) {
    stop("Not enough shared genes for comparison.")
  }
  
  # Sort both tables by P-value before taking the top 500.
  top_mine <- head(
    rownames(tt)[order(tt$P.Value)],
    500
  )
  
  top_shah <- head(
    rownames(shah)[order(shah$P.Value)],
    500
  )
  
  # Make a table with one row for this branch.
  data.frame(
    branch = branch_name,
    
    genes_tested = nrow(tt),
    
    genes_shared_with_shah = length(common),
    
    DEGs = length(my_deg),
    
    DEGs_up = sum(tt[my_deg, "logFC"] > 0),
    
    DEGs_down = sum(tt[my_deg, "logFC"] < 0),
    
    DEGs_FDR05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
    
    shah_DEGs_found = length(intersect(my_deg, shah_deg)),
    
    pct_shah_DEGs_found = round(
      100 * length(intersect(my_deg, shah_deg)) /
        length(shah_deg),
      1
    ),
    
    logFC_spearman = round(
      cor(
        tt[common, "logFC"],
        shah[common, "logFC"],
        method = "spearman",
        use = "complete.obs"
      ),
      4
    ),
    
    max_abs_logFC_diff = signif(
      max(
        abs(tt[common, "logFC"] - shah[common, "logFC"]),
        na.rm = TRUE
      ),
      3
    ),
    
    pct_same_direction = round(
      100 * mean(
        sign(tt[common, "logFC"]) ==
          sign(shah[common, "logFC"]),
        na.rm = TRUE
      ),
      1
    ),
    
    top500_overlap = length(
      intersect(top_mine, top_shah)
    )
  )
}


# Create a reference row using Shah's published results.

shah_row <- data.frame(
  branch = "Shah et al. (reference)",
  genes_tested = nrow(shah),
  genes_shared_with_shah = nrow(shah),
  DEGs = length(shah_deg),
  DEGs_up = sum(shah[shah_deg, "logFC"] > 0),
  DEGs_down = sum(shah[shah_deg, "logFC"] < 0),
  DEGs_FDR05 = sum(shah$adj.P.Val < 0.05, na.rm = TRUE),
  shah_DEGs_found = length(shah_deg),
  pct_shah_DEGs_found = 100,
  logFC_spearman = 1,
  max_abs_logFC_diff = 0,
  pct_same_direction = 100,
  top500_overlap = min(500, nrow(shah))
)


# Combine Shah's results and my two branches.

comparison <- rbind(
  shah_row,
  compare_results(results$shah, "Shah-style filter"),
  compare_results(results$standard, "Standard filter")
)

cat("\nComparison with Shah et al.:\n")
print(comparison)

write.csv(
  comparison,
  file.path(out_dir, "DE_comparison_with_shah.csv"),
  row.names = FALSE
)

cat("\nShah et al. Table 1 reports",
    "459 up / 746 down for GSE95587.\n")


# 9. Fold-change comparison plots -----------------------------------

# Compare the log2 fold changes from my analysis and Shah's.
#
# Each point represents one gene.
# The red diagonal shows where points would fall if the
# fold changes matched exactly.

png(
  file.path(out_dir, "DE_logFC_vs_shah.png"),
  width = 1400,
  height = 700
)

par(mfrow = c(1, 2))

for (b in names(results)) {
  
  common <- intersect(
    rownames(results[[b]]),
    rownames(shah)
  )
  
  plot(
    shah[common, "logFC"],
    results[[b]][common, "logFC"],
    pch = 16,
    cex = 0.4,
    col = rgb(0, 0, 0, 0.3),
    xlab = "logFC - Shah et al.",
    ylab = "logFC - my analysis",
    main = paste("GSE95587 -", b)
  )
  
  abline(0, 1, col = "red")
}

dev.off()


# 10. Volcano plots -------------------------------------------------

# A volcano plot shows:
# x-axis = log2 fold change
# y-axis = statistical significance (-log10 P-value)
#
# Red genes pass the DEG threshold and are higher in AD.
# Blue genes pass the threshold and are lower in AD.
# Grey genes do not pass both thresholds.

png(
  file.path(out_dir, "DE_volcano.png"),
  width = 1400,
  height = 700
)

par(mfrow = c(1, 2))

for (b in names(results)) {
  
  tt <- results[[b]]
  
  point_col <- ifelse(
    is_deg(tt) & tt$logFC > 0,
    "tomato",
    ifelse(
      is_deg(tt) & tt$logFC < 0,
      "steelblue",
      "grey70"
    )
  )
  
  plot(
    tt$logFC,
    -log10(tt$P.Value),
    pch = 16,
    cex = 0.4,
    col = point_col,
    xlab = "logFC (AD vs Control)",
    ylab = "-log10(P)",
    main = paste("GSE95587 -", b)
  )
  
  # Dashed lines show the DEG cutoffs.
  abline(
    v = c(-0.58, 0.58),
    h = -log10(0.05),
    lty = 2
  )
}

dev.off()


# 11. Save summary and software versions ----------------------------

# Save the main settings and comparison results.
# This will help with my capstone process documentation.

summary_lines <- c(
  paste("Dataset:", acc),
  paste("Date:", Sys.time()),
  "Model: limma-voom",
  "Design: ~ group (Control reference)",
  "Covariates: none",
  "DEG threshold: P < 0.05 and |logFC| > 0.58",
  paste("edgeR version:",
        as.character(packageVersion("edgeR"))),
  paste("limma version:",
        as.character(packageVersion("limma"))),
  "",
  "Gene list comparison:",
  paste("Genes in my Shah-style branch:", length(my_genes)),
  paste("Genes in Shah reference:", length(shah_genes)),
  paste("Genes shared:", length(common_genes)),
  "",
  "Differential expression comparison:",
  capture.output(print(comparison))
)

writeLines(
  summary_lines,
  file.path(out_dir, "DE_summary.txt")
)

# Save R and package versions.
writeLines(
  capture.output(sessionInfo()),
  file.path("logs", "sessionInfo_01_DE_GSE95587.txt")
)

cat("\nDE analysis complete!\n")
cat("Results saved to:", out_dir, "\n")
