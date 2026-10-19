#1st Load - GSE95587
#Download and load GSE95587 data: 117 brain samples, loads it into R. 
#Checks that the samples match their labels, makes 2 quality checks plots
#and save short summary. 
#Nothing analytical yet. No genes are filtered, no statistics run
#Load and check step

#0. Setup
#Load the two packages GEOquery, edgeR
#We only use it for "DGEList" (a container for counts) and to
#compute log CMP values
library(GEOquery)
library(edgeR)

options(timeout = 3600)  # allow slow downloads
#give up on downloads after 60 secs by default. Allow up to 1 hour

#store the dataset ID once, so its not mistyped later
acc     <- "GSE95587"
raw_dir <- file.path("data", "raw", acc) #downloaded count files
meta_dir <- file.path("data", "metadata", acc) #sample information
out_dir <- file.path("results", acc) #QC plots and summary 
#create each folder if it doesn't exist yet
for (d in c(raw_dir, meta_dir, out_dir, "logs")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

#1. Download the counts and the gene annotation from GEO 
#These are the NCBI-generated counts (GEO page -> "Download RNA-seq counts"), 
#the same files Shah et al. used for this cohort.
counts_file <- file.path(raw_dir, paste0(acc, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
annot_file  <- file.path(raw_dir, "Human.GRCh38.p13.annot.tsv.gz")

#download only when the files are not already saved
if (!file.exists(counts_file)) {
  download.file(
    paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", acc,
           "&format=file&file=", acc, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
    counts_file, mode = "wb")
}
if (!file.exists(annot_file)) {
  download.file(
    "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz",
    annot_file, mode = "wb")
}

#2. Load the counts and gene annotation
#read compressed counts, rows are genes and columns are samples 
counts <- read.delim(gzfile(counts_file), check.names = FALSE) 
rownames(counts) <- counts$GeneID
counts <- as.matrix(counts[, -1])
cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

#match each GeneID to its gene symbol, check missing or repeated symbols
annot <- read.delim(gzfile(annot_file), check.names = FALSE, quote = "")
symbol <- annot$Symbol[match(rownames(counts), annot$GeneID)]
cat("Genes with no symbol:", sum(is.na(symbol) | symbol == ""), "\n")
cat("Gene symbols that appear more than once:",
    sum(duplicated(symbol[!is.na(symbol)])), "\n")
#decision for later: how to handle duplicate symbols. No genes are changed here

#3. Load the sample metadata (diagnosis, age, sex, Braak, etc)
#get diagnosis and other sample details from GEO, save a copy for future runs
meta_file <- file.path(meta_dir, paste0(acc, "_pData_full.csv"))
if (!file.exists(meta_file)) {
  gse <- getGEO(acc, GSEMatrix = TRUE, destdir = meta_dir)
  pdata <- Biobase::pData(gse[[1]])
  write.csv(pdata, meta_file, row.names = FALSE)
}
pdata <- read.csv(meta_file, check.names = FALSE)
cat("Metadata:", nrow(pdata), "samples x", ncol(pdata), "columns\n")
cat("Sample characteristics available:\n")
print(grep(":ch1$", names(pdata), value = TRUE))

dx_col <- grep("diagnosis|disease|condition", grep(":ch1$", names(pdata), value = TRUE),
               ignore.case = TRUE, value = TRUE)[1]
if (is.na(dx_col)) stop("Couldn't find a diagnosis column; see the list above.")
cat("\nUsing diagnosis column:", dx_col, "\n")
print(table(pdata[[dx_col]], useNA = "ifany"))

#4. Match samples to dignoses
#Put metadata in the sampe sample order as the counts. Stop if IDs don't match
pdata <- pdata[match(colnames(counts), pdata$geo_accession), ]
stopifnot(identical(colnames(counts), pdata$geo_accession))

#Assign AD or Control from the diagnosis text. Leave unknown labels as NA
group <- ifelse(grepl("alzheimer|^ad$", pdata[[dx_col]], ignore.case = TRUE), "AD",
                ifelse(grepl("control|normal", pdata[[dx_col]], ignore.case = TRUE), "Control", NA))
cat("\nGroups used:\n")
print(table(group, useNA = "ifany"))
cat("Shah et al. Table 1 reports 84 AD and 33 controls.\n")

#5. Basic quality checks
#store counts in edgeR and calculate the total reads/sample (in millions)
dge <- DGEList(counts = counts, group = factor(group, levels = c("Control", "AD")))
lib_millions <- dge$samples$lib.size / 1e6
cat("\nLibrary size (millions of reads): min", round(min(lib_millions), 1),
    "median", round(median(lib_millions), 1), "max", round(max(lib_millions), 1), "\n")

#Log CPM adjusts for library size and makes expression values easier to compare
#This is for QC only. Raw counts are unchanged
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
png(file.path(out_dir, "qc_logCPM_boxplot.png"), width = 1600, height = 600)
boxplot(logcpm, las = 2, cex.axis = 0.5, outline = FALSE,
        col = ifelse(group == "AD", "tomato", "steelblue"),
        main = "GSE95587 log2-CPM before filtering (red = AD, blue = control)")
dev.off()

#Plot 1 - compare the distribution of gene expression across samples
png(file.path(out_dir, "qc_library_sizes.png"), width = 1600, height = 600)
barplot(lib_millions, las = 2, cex.names = 0.5,
        col = ifelse(group == "AD", "tomato", "steelblue"),
        main = "GSE95587 library sizes (millions)")
dev.off()

# Plot 2: compare total read counts across samples.
png(file.path(out_dir, "qc_library_sizes.png"), width = 1600, height = 600)
barplot(lib_millions, las = 2, cex.names = 0.5,
        col = ifelse(group == "AD", "tomato", "steelblue"),
        main = "GSE95587 library sizes (millions)")
dev.off()

#6. Save summary
#record sample counts, diagnosis labels, gene symbols, and file checksum
summary_lines <- c(
  paste("Dataset:", acc),
  paste("Run on:", format(Sys.time())),
  paste("Genes:", nrow(counts), " Samples:", ncol(counts)),
  paste("AD:", sum(group == "AD", na.rm = TRUE),
        " Control:", sum(group == "Control", na.rm = TRUE),
        " Unlabelled:", sum(is.na(group))),
  paste("Diagnosis column:", dx_col),
  paste("Duplicate gene symbols:", sum(duplicated(symbol[!is.na(symbol)]))),
  paste("Library size median (M):", round(median(lib_millions), 1)),
  paste("Counts file MD5:", unname(tools::md5sum(counts_file)))
)
writeLines(summary_lines, file.path(out_dir, "load_summary.txt"))

#save R and package versions so the analysis can be reproduced
writeLines(capture.output(sessionInfo()), file.path("logs", "sessionInfo_01_load_GSE95587.txt"))

cat("\nDone. Look in results/GSE95587/ for the summary and two QC plots.\n")

