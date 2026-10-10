# ------------------------ Load required packages -------------------------------
library(reshape2)
library(WGCNA)
library(limma)
library(pheatmap)
library(ggplot2)
library(grid)
library(RColorBrewer)

set.seed(1234)
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Working directory & input file ----------------------
# Please set your own working directory before running
# workDir <- "D:\\"
# setwd(workDir)

expFilePath <- "Sample_Type_Matrix_HighCV_Top30percent.csv"

if (!file.exists(expFilePath)) {
  stop("Input expression matrix file not found. Please check file path.")
}
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Data reading and preprocessing ---------------------------
tmp_head <- readLines(expFilePath, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")

rawData <- read.table(expFilePath, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(rawData) <- rawData[, 1]
exprData <- rawData[, -1]
exprDataCopy <- exprData

dimNames <- list(rownames(exprDataCopy), colnames(exprDataCopy))
numericData <- matrix(as.numeric(as.matrix(exprDataCopy)), nrow = nrow(exprDataCopy), dimnames = dimNames)

# log2‑transform if raw count data detected
if (max(numericData, na.rm = TRUE) > 1000) {
  numericData <- log2(numericData + 1)
}

normalizedData <- normalizeBetweenArrays(numericData)
filteredData <- normalizedData[apply(normalizedData, 1, sd, na.rm = TRUE) > 0.5, ]

# ------------------------ Auto‑detect sample groups & transpose ----------------------
sample_names <- colnames(filteredData)
group_info <- ifelse(grepl("_Con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_OA$", sample_names, ignore.case = TRUE), "Osteoarthritis", "Unknown"))

if (any(group_info == "Unknown")) stop("Unknown sample group detected. Check sample name suffix _Con/_OA.")

numControl <- sum(group_info == "Control")
numTreatment <- sum(group_info == "Osteoarthritis")

dataControl   <- filteredData[, which(group_info == "Control"), drop=FALSE]
dataTreatment <- filteredData[, which(group_info == "Osteoarthritis"), drop=FALSE]
combinedData  <- cbind(dataControl, dataTreatment)

exprMatrix    <- t(combinedData)

if (nrow(exprMatrix) < 10) {
  warning("Sample size is small; network robustness may be limited.")
}
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Sample quality check & clustering -------------------------
sampleCheck <- goodSamplesGenes(exprMatrix, verbose = 3)
if (!sampleCheck$allOK) {
  if (sum(!sampleCheck$goodGenes) > 0) {
    cat("Removed genes:", paste(names(exprMatrix)[!sampleCheck$goodGenes], collapse = ", "), "\n")
  }
  if (sum(!sampleCheck$goodSamples) > 0) {
    cat("Removed samples:", paste(rownames(exprMatrix)[!sampleCheck$goodSamples], collapse = ", "), "\n")
  }
  exprMatrix <- exprMatrix[sampleCheck$goodSamples, sampleCheck$goodGenes]
}

sampleDendro <- hclust(dist(exprMatrix), method = "average")
pdf("Sample_Clustering.pdf", width = 12, height = 9)
par(cex = 0.8)
par(mar = c(0, 4, 2, 0))
plot(sampleDendro, main = "Sample Clustering for Outlier Detection", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.2, cex.main = 2)
abline(h = 20000, col = "red", lwd = 2)
dev.off()

clusterCut <- cutreeStatic(sampleDendro, cutHeight = 20000, minSize = 10)
if (length(unique(clusterCut)) > 1) {
  keepSamples <- (clusterCut == 1)
  exprMatrix <- exprMatrix[keepSamples, ]
} else {
  cat("No obvious outlier samples detected; using all samples.\n")
}
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Trait data and sample matching ---------------------------
clinicalData <- data.frame(
  Normal = c(rep(1, numControl), rep(0, numTreatment)),
  Disease = c(rep(0, numControl), rep(1, numTreatment))
)
rownames(clinicalData) <- colnames(combinedData)

commonSamples <- intersect(rownames(exprMatrix), rownames(clinicalData))
exprMatrix <- exprMatrix[commonSamples, ]
clinicalData <- clinicalData[commonSamples, ]

if (nrow(exprMatrix) != nrow(clinicalData)) {
  stop("Mismatch between expression matrix and clinical trait sample count.")
}

sampleDendro2 <- hclust(dist(exprMatrix), method = "average")
traitColors <- numbers2colors(clinicalData, signed = FALSE)
pdf("Sample_Heatmap.pdf", width = 12, height = 12)
plotDendroAndColors(sampleDendro2, traitColors,
                    groupLabels = names(clinicalData),
                    main = "Sample Dendrogram and Trait Heatmap", dendroLabels = FALSE)
dev.off()
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Soft‑threshold selection & network construction ---------------------------
closeAllConnections()
enableWGCNAThreads()

powerVector <- 1:20
sftResult <- pickSoftThreshold(exprMatrix, powerVector = powerVector, verbose = 5)
optimalPower <- 14
cat("Soft‑threshold power set to:", optimalPower, "\n")

pdf("Scale_Independence_and_Mean_Connectivity.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sftResult$fitIndices[, 1], -sign(sftResult$fitIndices[, 3]) * sftResult$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale Independence Analysis", cex.axis = 1.2, cex.lab = 1.3)
text(sftResult$fitIndices[, 1], -sign(sftResult$fitIndices[, 3]) * sftResult$fitIndices[, 2],
     labels = powerVector, cex = cex1, col = "blue")
abline(h = 0.90, col = "red", lwd = 2)
points(optimalPower, -sign(sftResult$fitIndices[powerVector == optimalPower, 3]) *
         sftResult$fitIndices[powerVector == optimalPower, 2],
       col = "red", pch = 19, cex = 1.5)

plot(sftResult$fitIndices[, 1], sftResult$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Network Connectivity Analysis", cex.axis = 1.2, cex.lab = 1.3)
text(sftResult$fitIndices[, 1], sftResult$fitIndices[, 5],
     labels = powerVector, cex = cex1, col = "blue")
points(optimalPower, sftResult$fitIndices[powerVector == optimalPower, 5],
       col = "red", pch = 19, cex = 1.5)
dev.off()

adjacencyMatrix <- adjacency(exprMatrix, power = optimalPower)
TOMMatrix <- TOMsimilarity(adjacencyMatrix)
dissTOM <- 1 - TOMMatrix
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Gene clustering & module detection ---------------------------
geneDendro <- hclust(as.dist(dissTOM), method = "average")
pdf("Gene_Clustering.pdf", width = 12, height = 9)
plot(geneDendro, xlab = "", sub = "", main = "Gene Clustering Based on TOM", labels = FALSE, hang = 0.04)
dev.off()

minModuleSize <- 50
dynamicModuleLabels <- cutreeDynamic(dendro = geneDendro, distM = dissTOM,
                                     deepSplit = 2, pamRespectsDendro = FALSE,
                                     minClusterSize = minModuleSize)
moduleColors <- labels2colors(dynamicModuleLabels)

pdf("Dynamic_Tree_Modules.pdf", width = 8, height = 6)
plotDendroAndColors(geneDendro, moduleColors, "Dynamic Tree Cut Modules",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

moduleEigengenesList <- moduleEigengenes(exprMatrix, colors = moduleColors)
MEs <- moduleEigengenesList$eigengenes
moduleDiss <- 1 - cor(MEs)
moduleEigDendro <- hclust(as.dist(moduleDiss), method = "average")
mergeThreshold <- 0.25

pdf("Module_Clustering.pdf", width = 8, height = 6)
plot(moduleEigDendro, main = "Clustering of Module Eigengenes", xlab = "", sub = "")
abline(h = mergeThreshold, col = "red", lwd = 2)
dev.off()

mergeResult <- mergeCloseModules(exprMatrix, moduleColors, cutHeight = mergeThreshold, verbose = 3)
mergedModuleColors <- mergeResult$colors
mergedMEs <- mergeResult$newMEs

pdf("Merged_Modules_Comparison.pdf", width = 10, height = 6)
plotDendroAndColors(geneDendro, cbind(moduleColors, mergedModuleColors),
                    c("Original Modules", "Merged Modules"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Module Comparison Pre‑ and Post‑Merging")
dev.off()

moduleColors <- mergedModuleColors
moduleEigengenes <- mergedMEs
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Module‑trait relationship analysis ---------------------------
moduleTraitCor <- cor(moduleEigengenes, clinicalData, use = "p")
moduleTraitPvalues <- corPvalueStudent(moduleTraitCor, nrow(exprMatrix))

corDF <- melt(moduleTraitCor)
pvalDF <- melt(moduleTraitPvalues)
heatmapDF <- merge(corDF, pvalDF, by = c("Var1", "Var2"))
colnames(heatmapDF) <- c("Module", "Trait", "Correlation", "Pvalue")

pdf("Module_Trait_Heatmap.pdf", width = 7, height = 6)
ggplot(heatmapDF, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name = "Correlation") +
  geom_text(aes(label = sprintf("%.2f\n(%.1e)", Correlation, Pvalue)), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  ggtitle("Module‑Trait Correlation Heatmap")
dev.off()

write.table(heatmapDF, file = "Module_Trait_Correlation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Module membership (MM) & Gene significance (GS) -----------------------
geneModuleMembership <- as.data.frame(cor(exprMatrix, moduleEigengenes, use = "p"))
geneMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exprMatrix)))

geneTraitSignificance <- as.data.frame(cor(exprMatrix, clinicalData, use = "p"))
geneGSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exprMatrix)))

geneNames <- colnames(exprMatrix)
geneInfo <- data.frame(Gene = geneNames, Module = moduleColors)

for (mod in colnames(moduleEigengenes)) {
  geneInfo[, paste("MM_", mod, sep = "")] <- geneModuleMembership[, mod]
  geneInfo[, paste("p.MM_", mod, sep = "")] <- geneMMPvalue[, mod]
}

for (trait in colnames(clinicalData)) {
  geneInfo[, paste("GS_", trait, sep = "")] <- geneTraitSignificance[, trait]
  geneInfo[, paste("p.GS_", trait, sep = "")] <- geneGSPvalue[, trait]
}

geneInfo <- geneInfo[order(geneInfo$Module), ]
write.table(geneInfo, file = "GeneInfo_Modules.txt", sep = "\t", row.names = FALSE, quote = FALSE)
step <- step + 1
setTxtProgressBar(pb, step)

# ------------------------ Export gene list per module ------------------------
uniqueModules <- unique(moduleColors)
for (mod in uniqueModules) {
  moduleGenes <- geneNames[moduleColors == mod]
  outFileName <- paste0("ModuleGenes_", mod, ".txt")
  write.table(moduleGenes, file = outFileName, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
step <- step + 1
setTxtProgressBar(pb, step)

rm(rawData, exprData, exprDataCopy, numericData, normalizedData, filteredData)
close(pb)

moduleSizes <- table(moduleColors)
moduleSizeDF <- as.data.frame(moduleSizes)
colnames(moduleSizeDF) <- c("Module", "GeneCount")

pdf("Module_Gene_Counts.pdf", width = 8, height = 6)
ggplot(moduleSizeDF, aes(x = Module, y = GeneCount, fill = Module)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Gene Counts per Module") +
  xlab("Module") +
  ylab("Gene Count") +
  theme(text = element_text(size = 12))
dev.off()

targetTrait <- "Disease"
allModules <- unique(moduleColors)
for (mod in allModules) {
  pdf(file = paste0("MM_vs_GS_", mod, ".pdf"), width = 6, height = 6)
  moduleGenes <- (moduleColors == mod)
  mmColName <- paste0("ME", mod)
  gsColName <- targetTrait

  MM <- as.numeric(geneModuleMembership[moduleGenes, mmColName])
  GS <- as.numeric(geneTraitSignificance[moduleGenes, gsColName])

  corTest <- cor.test(MM, GS)
  corVal <- corTest$estimate
  pVal <- corTest$p.value

  plot(MM, GS,
       xlab = paste("Module Membership in", mod, "module"),
       ylab = paste("Gene Significance for", targetTrait),
       main = paste0("Module: ", mod,
                     "\ncor=", signif(corVal, 3),
                     ", p=", format(pVal, scientific = TRUE, digits = 2)),
       pch = 21, bg = adjustcolor(mod, alpha.f = 0.6),
       col = "black", cex = 1.5)
  abline(lm(GS ~ MM), col = "blue", lwd = 2, lty = 2)
  abline(v = 0.8, h = 0.2, col = "orange", lty = 3, lwd = 1.5)
  dev.off()
}

df <- data.frame(
  Module = moduleColors,
  MM = sapply(1:length(moduleColors), function(i) geneModuleMembership[i, paste0("ME", moduleColors[i])]),
  GS = geneTraitSignificance[, "Disease"]
)

pdf("Boxplot_MM.pdf", width = 7, height = 6)
p_box <- ggplot(df, aes(x = Module, y = MM, fill = Module)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Boxplot of Module Membership across Modules", y = "Module Membership") +
  theme(legend.position = "none")
print(p_box)
dev.off()

pdf("ViolinPlot_GS.pdf", width = 7, height = 6)
p_violin <- ggplot(df, aes(x = Module, y = GS, fill = Module)) +
  geom_violin(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Violin Plot of Gene Significance across Modules", y = "Gene Significance") +
  theme(legend.position = "none")
print(p_violin)
dev.off()

# ------------------------ Module‑module correlation heatmap ---------------------------
moduleCor <- cor(moduleEigengenes, use = "p")
moduleCor_melt <- melt(moduleCor)

pdf("Module_Module_Correlation.pdf", width = 7, height = 6)
ggplot(moduleCor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  theme_minimal() +
  ggtitle("Module‑to‑Module Correlation Heatmap") +
  xlab("Module Eigengenes") +
  ylab("Module Eigengenes") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

# TOM plot for black module
targetModule <- "black"
inModule <- (moduleColors == targetModule)
modGenes <- geneNames[inModule]
modTOM <- TOMMatrix[inModule, inModule]
geneDendroModule <- hclust(as.dist(1 - modTOM), method = "average")

pdf(paste0("TOMplot_", targetModule, ".pdf"), width = 10, height = 10)
TOMplot(modTOM, geneDendroModule, main = paste("TOM Plot for", targetModule, "Module"))
dev.off()

mod_colors_palette <- c(
  "brown"="#A0522D","blue"="#1E90FF","yellow"="#FFD700","turquoise"="#20B2AA",
  "grey"="#BEBEBE","black"="#000000","magenta"="#FF00FF",
  "green"="#228B22","red"="#E64B35","pink"="#FFB6C1",
  "purple"="#6A5ACD","salmon"="#FA8072","orange"="#FFA500",
  "cyan"="#00CED1","lightyellow"="#FFFACD",
  "darkorange"="#FF8C00","royalblue"="#4169E1","darkgrey"="#A9A9A9"
)

all_modules <- unique(moduleColors)

for(targetModule in all_modules){
  cat("Processing module: ", targetModule, "\n")
  mod_genes <- colnames(exprMatrix)[moduleColors == targetModule]
  if(length(mod_genes)<2) {
    cat("Module",targetModule,"has fewer than 2 genes, skip.\n")
    next
  }
  mod_color <- ifelse(targetModule %in% names(mod_colors_palette), mod_colors_palette[targetModule], "#888888")
  sample_list <- rownames(exprMatrix)
  me_col <- paste0("ME", targetModule)
  if(!(me_col %in% colnames(MEs))) {
    cat("Eigengene for module", targetModule,"not found, skip.\n")
    next
  }
  mod_expr <- exprMatrix[sample_list, mod_genes, drop=FALSE]
  eigengene <- MEs[sample_list, me_col, drop=TRUE]
  bar_df <- data.frame(
    Sample = factor(sample_list, levels=sample_list),
    Eigengene = eigengene
  )

  heatmap_grob <- grid.grabExpr({
    pheatmap(
      t(mod_expr),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_rownames = FALSE,
      show_colnames = FALSE,
      scale = "row",
      color = colorRampPalette(c("skyblue", "white", "orange"))(100),
      border_color = NA,
      legend = FALSE,
      fontsize = 10
    )
  })

  bar_p <- ggplot(bar_df, aes(x = Sample, y = Eigengene)) +
    geom_bar(stat = "identity", width = 1, fill = mod_color) +
    labs(y = "Eigengene", x = "Sample") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 6),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      plot.margin = margin(0, 5, 0, 5),
      axis.line = element_line(linewidth = 0.6)
    )

  pdf_file <- paste0("Module_", targetModule, "_heatmap_eigengene.pdf")
  pdf(pdf_file, width = 9, height = 7)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1,
                                             heights = unit(c(0.58, 0.38), "npc"),
                                             widths = unit(1, "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_grob)
  grid.text(targetModule, x = 0.5, y = unit(0.98, "npc"),
            gp = gpar(fontface = "bold", fontsize = 20, col = mod_color))
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(bar_p, newpage=FALSE)
  popViewport()
  dev.off()
  cat("Saved PDF:", pdf_file, "\n")
}

cat("★ All‑module output finished!\n")
getwd()
