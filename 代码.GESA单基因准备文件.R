

# ---- 0. Working directory (edit to your local path) ----
setwd("D:\\")

# ---- 1. Define the target gene ----
geneName <- "ANTXR1"                    # set the gene name directly here
cat("Looking for gene:", geneName, "\n")

# Optional: read gene name from command line
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) > 0) {
#   geneName <- args[1]
#   cat("Gene from command line:", geneName, "\n")
# } else {
#   geneName <- "ANTXR1"
# }

# Uppercase for case-insensitive matching
geneName_upper <- toupper(trimws(geneName))

# ---- 2. Initialize variables ----
geneArr <- c()
indexs <- c()          # column indices of samples
hash <- list()
colNum <- 0

# ---- 3. Read the gene-level expression matrix ----
# Row 1 = sample names (starting at column 2), column 1 = gene symbols
data <- read.table("geneMatrix.txt", header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, comment.char = "")
cat("Data dimension:", dim(data), "\n")

# First pass: collect sample column indices and build sample hash
for (i in 2:ncol(data)) {
  sampleName <- trimws(data[1, i])
  indexs <- c(indexs, i)
  hash[[sampleName]] <- 1
  colNum <- colNum + 1
}
cat("Number of samples:", colNum, "\n")
cat("First 5 sample names:", names(hash)[1:min(5, length(hash))], "\n")

# ---- 4. Locate the target gene in the matrix ----
all_genes_upper <- toupper(trimws(data[2:nrow(data), 1]))
cat("Total genes in matrix:", length(all_genes_upper), "\n")

found <- FALSE
target_row <- NULL

if (geneName_upper %in% all_genes_upper) {
  matching_rows <- which(all_genes_upper == geneName_upper)
  original_gene_names <- data[matching_rows + 1, 1]
  cat("Exact match found:", original_gene_names, "at row(s)",
      matching_rows + 1, "\n")
  target_row <- matching_rows[1] + 1
  found <- TRUE
} else {
  cat("No exact match, trying fuzzy search...\n")
  matches <- grep(geneName_upper, all_genes_upper)
  if (length(matches) > 0) {
    target_row <- matches[1] + 1
    found <- TRUE
    cat("Using fuzzy match:", data[target_row, 1], "at row", target_row, "\n")
  }
}

# Stop if the gene cannot be found
if (!found) {
  cat("Error: gene", geneName, "not found in the matrix.\n")
  print(head(data[2:min(21, nrow(data)), 1], 20))
  stop("Gene lookup failed.")
}

# Extract the expression values of the target gene
geneArr <- as.numeric(data[target_row, indexs])
names(geneArr) <- trimws(data[1, indexs])
cat("Expression summary for", data[target_row, 1], ":\n")
print(summary(geneArr))

# Save expression values to CSV
result_df <- data.frame(Sample = names(geneArr), Expression = geneArr)
output_file <- paste0(geneName, "_expression_values.csv")
write.csv(result_df, output_file, row.names = FALSE)
cat("Expression values saved to:", output_file, "\n")

# Expose as a variable in the global environment
assign(paste0(geneName, "_expr"), geneArr, envir = .GlobalEnv)
cat("Variable", paste0(geneName, "_expr"), "created in global env.\n")

# ---- 5. Compute the median expression of the target gene ----
firstGeneVal <- geneArr[1]
geneMed <- median(geneArr)

# Total number of gene rows (exclude header row)
rowNum <- nrow(data) - 1

# ---- 6. Open output files ----
gct_file <- file(paste0(geneName, ".gct"), "w")
cls_file <- file(paste0(geneName, ".cls"), "w")

# ---- 7. Write GCT header ----
cat("#1.2\n", file = gct_file)
cat(paste(rowNum, colNum, sep = "\t"), "\n", file = gct_file)

# Column titles: all sample names
sampleNames <- sapply(2:ncol(data), function(i) trimws(data[1, i]))
cat(paste(sampleNames, collapse = "\t"), "\n", file = gct_file)

# ---- 8. Write CLS header ----
# [numSamples, numClasses, 1]
cat(paste(colNum, 2, 1, sep = "\t"), "\n", file = cls_file)
# Class label order is based on the target gene's first sample vs. median
if (firstGeneVal > geneMed) {
  cat("#\tl\th\n", file = cls_file)
} else {
  cat("#\th\tl\n", file = cls_file)
}

# ---- 9. Write expression rows and assign class labels ----
typeArr <- c()
for (i in 2:nrow(data)) {
  # Keep only the part before "|" in gene symbols
  symbolName <- gsub("\\|.*", "", data[i, 1])
  symbolName <- trimws(symbolName)

  # Write gene symbol and fixed "na" description field
  cat(symbolName, "\tna", file = gct_file)

  # Write expression values across all samples
  for (col in indexs) {
    cat("\t", data[i, col], file = gct_file)
  }
  cat("\n", file = gct_file)

  # For the target gene, classify each sample as high (h) or low (l)
  currentGene <- toupper(trimws(data[i, 1]))
  if (currentGene == geneName_upper) {
    for (col in indexs) {
      if (as.numeric(data[i, col]) > geneMed) {
        typeArr <- c(typeArr, "h")
      } else {
        typeArr <- c(typeArr, "l")
      }
    }
  }
}

# Write class labels in sample order
cat(paste(typeArr, collapse = "\t"), "\n", file = cls_file)

close(gct_file)
close(cls_file)
cat("GCT and CLS files written for gene:", geneName, "\n")

# ---- 10. Post-process all .gct files ----
# Add "NAME" and "DESCRIPTION" headers to row 3 of every .gct file
gct_files <- list.files(pattern = "\\.gct$")
if (length(gct_files) == 0) {
  stop("No .gct files found in the current directory.")
}

for (file in gct_files) {
  lines <- readLines(file)
  if (length(lines) < 3) {
    warning(paste("File", file, "has fewer than 3 lines, skipped."))
    next
  }
  lines[3] <- paste("NAME", "DESCRIPTION", lines[3], sep = "\t")
  writeLines(lines, file)
  cat("Updated file:", file, "\n")
}
cat("All done.\n")
