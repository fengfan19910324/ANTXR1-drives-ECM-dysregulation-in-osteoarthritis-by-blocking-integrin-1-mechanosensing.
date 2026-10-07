# 1. Set working directory
setwd("D:\\")

# 2. Set parameter: target row for gene annotation (1‑based), convert to 0‑based index
inputRow <- 11
targetColIdx <- inputRow - 1

# 3. Define input and output file paths
exprFilePath <- "GSE.txt"
platformFilePath <- "GPL.txt"
outputFilePath <- "geneMatrix.txt"

# 4. Read expression matrix (GSE.txt)
cat("Loading expression data file:", exprFilePath, "...\n")
exprData <- read.delim(
  exprFilePath,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  comment.char = "!"
)
colnames(exprData)[1] <- "ProbeID"
cat("Expression data loaded,", ncol(exprData) - 1, "samples detected.\n")

# 5. Read platform annotation file (GPL.txt)
cat("Loading platform annotation file:", platformFilePath, "...\n")
platformData <- read.delim(
  platformFilePath,
  header = FALSE,
  sep = "\t",
  quote = "\"",
  comment.char = "#",
  stringsAsFactors = FALSE
)

# 6. Build probe‑to‑gene symbol mapping with progress bar
cat("Processing platform annotation...\n")
totalRows <- nrow(platformData)
geneMapping <- list()
pbPlat <- txtProgressBar(min = 0, max = totalRows, style = 3)

for (row in 1:totalRows) {
  # Skip empty rows / header lines starting with "ID" or "!"
  if (all(platformData[row, ] == "") || grepl("^(ID|\\!)", platformData[row, 1])) {
    setTxtProgressBar(pbPlat, row)
    next
  }

  # Extract target annotation column if available
  if (ncol(platformData) >= (targetColIdx + 1)) {
    rawGene <- platformData[row, targetColIdx + 1]
    # Keep entries without whitespace
    if (rawGene != "" && !grepl(".+\\s+.+", rawGene)) {
      # Split "///" composite annotations, take the first entry
      cleanGene <- sub("(.+?)///(.+)", "\\1", rawGene)
      cleanGene <- gsub('"', '', cleanGene)
      geneMapping[[as.character(platformData[row, 1])]] <- cleanGene
    }
  }
  setTxtProgressBar(pbPlat, row)
}
close(pbPlat)
cat("\nPlatform processing finished, total valid probe‑gene pairs:", length(geneMapping), "\n")

# 7. Merge expression matrix with probe‑gene mapping
cat("Merging expression data with gene annotation...\n")
mappingDF <- data.frame(
  ProbeID = names(geneMapping),
  geneSymbol = unlist(geneMapping),
  stringsAsFactors = FALSE
)
mergedData <- merge(exprData, mappingDF, by = "ProbeID")
cat("Merge completed,", nrow(mergedData), "probes with gene annotation.\n")

# 8. Aggregate multiple probes per gene: compute mean expression
cat("Aggregating probes to gene‑level expression (mean value)...\n")
sampleCols <- setdiff(colnames(mergedData), c("ProbeID", "geneSymbol"))
geneGroups <- split(mergedData, mergedData$geneSymbol)
geneNamesUnique <- names(geneGroups)
totalGenes <- length(geneNamesUnique)
resultList <- vector("list", totalGenes)

pbAgg <- txtProgressBar(min = 0, max = totalGenes, style = 3)
for (i in seq_along(geneNamesUnique)) {
  geneName <- geneNamesUnique[i]
  groupData <- geneGroups[[i]]
  means <- colMeans(groupData[, sampleCols, drop = FALSE], na.rm = TRUE)
  resultList[[i]] <- c(geneSymbol = geneName, means)
  setTxtProgressBar(pbAgg, i)
}
close(pbAgg)

# Convert list output to data frame and cast numeric columns
aggData <- do.call(rbind, resultList)
aggData <- as.data.frame(aggData, stringsAsFactors = FALSE)
for (col in colnames(aggData)[-1]) {
  aggData[[col]] <- as.numeric(aggData[[col]])
}

# Sort output by gene symbol
aggData <- aggData[order(aggData$geneSymbol), ]
cat("Aggregation finished,", nrow(aggData), "unique genes obtained.\n")

# 9. Export gene‑level expression matrix
cat("Writing output to:", outputFilePath, "...\n")
write.table(
  aggData,
  file = outputFilePath,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("Output file successfully generated:", outputFilePath, "\n")
