#01_preprocess_GSE95587.R

#Preprocess GSE95587 before DE analysis

#Preprocessing here is: 
#1. Give every gene a usable name 
#2. Normalize counts using TMM
#3. Filter low expression genes in 2 ways
#4. Compare how many genes each filter keeps
#5. Make MDS plots to check sample patterns
#6. Save both versions

#The goal is to compare Shah et al.'s described filtering method 
#with the standard edgeR approach

library(edgeR)

#1. Load the data from script 01
out_dir <- file.path("results", "GSE95587")
proc_dir <- file.path("data", "processed", "GSE95587")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

#2. Give every gene a usable name
#script 1 matched geneIDs to gene symbols but did not change them, 
#so I'll remove genes without symbols and make duplicate names unique

has_name <- !is.na(symbol) & symbol !=""
counts2 <- counts[has_name, , drop=FALSE]
rownames(counts2) <- make.unique(as.character(symbol[has_name]))

cat("Genes without symbols:", sum(!has_name), "\n")
cat("Genes remaining:", nrow(counts2), "\n")

#3. Set up the groups
#script 01 assigned AD and control labels, so I only need to 
#convert them into a factor for edgeR
#"Control" is 1st since it will be the reference group, 
#when comparing AD against control later

group2 <- factor(group, levels = c("Control", "AD"))

if(any(is.na(group2))) {
  stop("Some samples have missing or unexpected diagnosis labels")
}

stopifnot(ncol(counts2) == length(group2))

cat("Samples:", ncol(counts2), "\n")
print(table(group2))

#4. TMM Normalization
#create edgeR object with the counts and groupd levels
#DGElist keeps the count matrix and sample information together
dge <- DGEList(counts = counts2, group = group2)

dge <- calcNormFactors(dge, method = "TMM")

#5. Filter low expression genes
#5.1 Shah et al filtering
#shah et al. kept genes with Log2 CPM (counts per million) > 1 
#in at least as many samples as the smaller diagnosis group
#so min_n should be 33 (control) for this dataset

min_n <- min(table(group2))
logcpm <- cpm(dge, log=TRUE, prior.count = 1)

keep_shah <- rowSums(logcpm > 1) >= min_n

#B Standard edgeR filtering
#filterByExpr() is edgeR's recommend filtering function
#it uses the same groups and library sizes to decide which
#genes have enough expression for DE testing

keep_std <- filterByExpr(dge, group=group2)

#compare two filters
cat("\nGenes before filtering:", nrow(dge), "\n")
cat("Shah-style filter:", sum(keep_shah), "\n")
cat("Standard edgeR filter:", sum(keep_std), "\n")
cat("Kept by both:", sum(keep_shah & keep_std), "\n")

#6. Create and save both branches
#keep only genes that passed each filter
#recalculate library sizes and TMM factors after filtering
dge_shah <- dge[keep_shah, , keep.lib.sizes = FALSE]
dhe_shah <- calcNormFactors(dge_shah, method = "TMM")

dge_std <- dge[keep_std, , keep.lib.sizes = FALSE]
dhe_std <- calcNormFactors(dge_std, method = "TMM")

#save 2 versions in a names list
#similar to Python dictionary with 2 keys
branches <- list(
  shah = dge_shah,
  standard = dge_std
)

saveRDS(branches, file.path(proc_dir, "dge_branches.rds"))

#7. MDS plots
#MDS shows how similar or different the samples are overall
#Im looking for unusual samples or possible group patterns
#red - AD, blue - Control

for (b in names(branches)) {
  
  logcpm_plot <- cpm(branches[[b]], log = TRUE, prior.count = 1)
  
  print(dim(logcpm_plot))
  
  png(
    file.path(out_dir, paste0("qc_MDS_", b, ".png")),
    width = 900,
    height = 800
  )
  
  plotMDS(
    logcpm_plot,
    labels = group2,
    col = ifelse(group2 == "AD", "tomato", "steelblue"),
    main = paste("GSE95587 MDS -", b, "filter")
  )
  
  dev.off()
}

#8. Save preprocessing summary
#save results so I can document the decisions and compare #s
summary_lines <- c(
  paste("Dataset:", acc),
  paste("Date:", Sys.time()),
  paste("Genes before filtering:", nrow(dge)),
  paste("Shah-style filter:", sum(keep_shah)),
  paste("Standard filter:", sum(keep_std)),
  paste("Genes kept by both:", sum(keep_shah & keep_std)),
  paste("Normalization: TMM"),
  paste("edgeR version:", as.character(packageVersion("edgeR")))
)

writeLines(
  summary_lines,
  file.path(out_dir, "preprocess_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path("logs", "sessionInfo_02_preprocess_GSE95587.txt")
)

cat("\nPreprocessing complete!\n")
cat("Processed data saved to:", proc_dir, "\n")
