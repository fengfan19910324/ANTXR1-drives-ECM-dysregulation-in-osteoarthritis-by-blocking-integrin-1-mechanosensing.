# ===================================
# Step 0: Parameter setup & package loading
# ===================================
# Analysis parameters
filterP <- 0.05
filterQ <- 0.05
colorBy <- "p.adjust"
GO_Dir <- "D:\\"
geneFile <- "IntersectionGenes.csv"

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(pbapply)
suppressMessages(library(stringr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

# ==============================
# Step 1: Directory check & input reading
# ==============================
if (!dir.exists(GO_Dir)) {
  stop("Directory does not exist. Please verify GO_Dir.")
}
setwd(GO_Dir)

if (!file.exists(geneFile)) {
  stop("Input gene file not found. Check path or filename.")
}

file_content <- readLines(geneFile, n = 1)
sep_used <- ifelse(str_detect(file_content, "\t"), "\t", ",")
geneTbl <- read.table(geneFile, header = FALSE, sep = sep_used,
                      stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(geneTbl) < 2) {
  stop("Too few genes loaded. Please check input file content.")
}
cat("Step 1/6: Input data loaded\n")

# ==============================
# Step 2: Gene‑ID conversion & filtering
# ==============================
allGenes <- unique(as.character(geneTbl[, 1]))
gene2entrez <- mget(allGenes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezVec <- unlist(gene2entrez)
entrezVec <- as.character(entrezVec)
entrezVec <- entrezVec[!is.na(entrezVec) & entrezVec != "NA"]

if (length(entrezVec) == 0) {
  stop("No valid Entrez IDs retrieved. Check input gene symbols.")
}
cat("Step 2/6: Gene‑ID conversion done, valid genes:", length(entrezVec), "\n")

# ==============================
# Step 3: GO enrichment analysis
# ==============================
pboptions(type = "timer")
cat("Step 3/6: Running GO enrichment, please wait...\n")

goResList <- pblapply(1, function(x) {
  enrichGO(gene = entrezVec,
           OrgDb = org.Hs.eg.db,
           pvalueCutoff = 1,
           qvalueCutoff = 1,
           ont = "all",
           readable = TRUE)
})
GOres <- goResList[[1]]
GOresDF <- as.data.frame(GOres)

if (nrow(GOresDF) == 0) {
  stop("No GO enrichment terms returned.")
}

GOresDF_filt <- GOresDF %>% filter(pvalue < filterP & p.adjust < filterQ)
if (nrow(GOresDF_filt) == 0) {
  warning("No significant GO terms after filtering; full result will be exported.")
}

write.table(GOresDF_filt, file = "GO_Analysis_Results_Filtered.csv",
            sep = ",", quote = FALSE, row.names = FALSE)
write.table(GOresDF, file = "GO_Analysis_Results_All.csv",
            sep = ",", quote = FALSE, row.names = FALSE)
cat("Step 3/6: GO enrichment finished and saved\n")

# ==============================
# Step 4: Barplot & dotplot
# ==============================
cat("Step 4/6: Generating GO barplot and bubble plot...\n")

pdf("GO_Barplot_v2.pdf", width = 8, height = 7)
barp <- tryCatch({
  barplot(GOres, drop = TRUE, showCategory = 10, label_format = 80,
          split = "ONTOLOGY", color = colorBy) +
    facet_grid(ONTOLOGY ~ ., scales = "free")
}, error = function(e) {
  plot.new()
  title("Barplot rendering failed")
})
print(barp)
dev.off()

pdf("GO_BubblePlot_v2.pdf", width = 8, height = 7)
bubblep <- tryCatch({
  dotplot(GOres, showCategory = 10, orderBy = "GeneRatio", label_format = 80,
          split = "ONTOLOGY", color = colorBy) +
    facet_grid(ONTOLOGY ~ ., scales = "free")
}, error = function(e) {
  plot.new()
  title("Bubble‑plot rendering failed")
})
print(bubblep)
dev.off()
cat("Step 4/6: Barplot / bubble plot exported\n")

# ==============================
# Step 5: Custom barplot
# ==============================
cat("Step 5/6: Generating custom GO barplot...\n")

dat_show <- GOresDF_filt %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
if (nrow(dat_show) == 0) {
  dat_show <- GOresDF %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
}

pdf("GO_CustomBarplot.pdf", width = 11, height = 8)
p1 <- tryCatch({
  ggbarplot(dat_show, x = "Description", y = "Count", fill = "ONTOLOGY",
            color = "white", palette = "jco", legend = "right",
            sort.val = "desc", sort.by.groups = TRUE) +
    rotate_x_text(65) +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(size = 10, color = "black")) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0))
}, error = function(e) {
  plot.new()
  title("Custom barplot rendering failed")
})
print(p1)
dev.off()
cat("Step 5/6: Custom barplot saved\n")

# ==============================
# Step 6: GO circular plot
# ==============================
cat("Step 6/6: Generating GO circular plot...\n")

ontColors <- c("#0093D1", "#ECB64B", "#95B1AE")
GO_sorted <- GOresDF_filt[order(GOresDF_filt$p.adjust), ]
if (nrow(GO_sorted) == 0) {
  GO_sorted <- GOresDF[order(GOresDF$p.adjust), ]
}

BP_data <- head(GO_sorted[GO_sorted$ONTOLOGY == "BP", , drop = FALSE], 6)
CC_data <- head(GO_sorted[GO_sorted$ONTOLOGY == "CC", , drop = FALSE], 6)
MF_data <- head(GO_sorted[GO_sorted$ONTOLOGY == "MF", , drop = FALSE], 6)
GO_circ <- rbind(BP_data, CC_data, MF_data)

if (nrow(GO_circ) == 0) {
  warning("Insufficient terms for circular plot, skipped.")
} else {
  bgGenes <- as.numeric(sapply(strsplit(GO_circ$BgRatio, "/"), `[`, 1))
  selGenes <- as.numeric(sapply(strsplit(GO_circ$GeneRatio, "/"), `[`, 1))
  richFactor <- selGenes / bgGenes
  logP <- -log10(GO_circ$pvalue)
  reds <- brewer.pal(8, "Reds")
  logP_colFun <- colorRamp2(c(0, 2, 4, 6, 8, 10, 15, 20), reds)
  barCol <- logP_colFun(logP)

  circ_df <- data.frame(GO = GO_circ$ID, start = 1, end = max(bgGenes))
  rownames(circ_df) <- circ_df$GO

  pdf("GO_CircularPlot_v2.pdf", width = 10, height = 10)
  par(omi = c(0.1, 0.1, 0.1, 1.5))
  circos.clear()
  circos.par(track.margin = c(0.01, 0.01))
  circos.genomicInitialize(circ_df, plotType = "none")

  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    secInd <- get.cell.meta.data("sector.index")
    circos.text(mean(get.cell.meta.data("xlim")),
                mean(get.cell.meta.data("ylim")),
                secInd, cex = 0.75, facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.09, bg.border = NA,
  bg.col = ontColors[as.numeric(as.factor(GO_circ$ONTOLOGY))])

  for (si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
                major.at = seq(0, max(bgGenes), by = 100), labels.facing = "clockwise")
  }
  circos.clear()

  midLegend <- Legend(
    labels = c('Total number of genes', 'Number of selected genes', 'Rich Factor (0‑1)'),
    type = "points", pch = c(15, 15, 17),
    legend_gp = gpar(col = c('#ECB64B', '#FF3B6D', ontColors[1])),
    title = "", nrow = 3, size = unit(3, "mm"))
  sz <- unit(1, "snpc")
  draw(midLegend, x = sz * 0.42)

  mainLeg <- Legend(labels = c("BP", "CC", "MF"), type = "points", pch = 15,
                    legend_gp = gpar(col = ontColors), title_position = "topcenter",
                    title = "GO Category", nrow = 3, size = unit(3, "mm"),
                    grid_height = unit(5, "mm"), grid_width = unit(5, "mm"))

  logPLeg <- Legend(
    labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
    type = "points", pch = 16, legend_gp = gpar(col = reds),
    title = "-log10(P)", title_position = "topcenter",
    grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), size = unit(3, "mm"))

  leg_combine <- packLegend(mainLeg, logPLeg)
  draw(leg_combine, x = sz * 0.85, y = sz * 0.55, just = "left")
  dev.off()
  cat("Circular plot finished\n")
}

# =========================
# Step 7: Enhanced GO enrichment circle plot
# =========================
cat("Step 7/7: Generating enhanced GO enrichment circle plot...\n")

go_type_colors <- c("#5EC2AC", "#B6B129", "#F98E71")
go_res_ranked <- GOresDF[order(GOresDF$p.adjust), ]
go_res_sig <- go_res_ranked[go_res_ranked$pvalue < 0.05, , drop = FALSE]

go_bp <- head(go_res_sig[go_res_sig$ONTOLOGY == "BP", , drop = FALSE], 6)
go_cc <- head(go_res_sig[go_res_sig$ONTOLOGY == "CC", , drop = FALSE], 6)
go_mf <- head(go_res_sig[go_res_sig$ONTOLOGY == "MF", , drop = FALSE], 6)
go_ring_data <- rbind(go_bp, go_cc, go_mf)

if (nrow(go_ring_data) == 0) {
  warning("No significant GO terms for enhanced circle plot, skipped.")
} else {
  bg_gene_count <- as.numeric(sapply(strsplit(go_ring_data$BgRatio, "/"), `[`, 1))
  selected_gene_count <- as.numeric(sapply(strsplit(go_ring_data$GeneRatio, "/"), `[`, 1))
  rich_factor <- selected_gene_count / bg_gene_count
  log_pval <- -log10(go_ring_data$pvalue)
  pval_palette <- brewer.pal(n = 8, name = "Reds")
  logpval_map <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10, 15, 20), colors = pval_palette)
  bg_color_map <- logpval_map(log_pval)

  circ_df <- data.frame(GOID = go_ring_data$ID, start = 1, end = max(bg_gene_count))
  rownames(circ_df) <- circ_df$GOID

  track_bed_bg <- data.frame(GOID = go_ring_data$ID, start = 1, end = bg_gene_count,
                             n_bg = bg_gene_count, col_bg = bg_color_map)
  track_bed_sel <- data.frame(GOID = go_ring_data$ID, start = 1, end = selected_gene_count,
                              n_sel = selected_gene_count)
  track_bed_ratio <- data.frame(GOID = go_ring_data$ID, start = 1, end = max(bg_gene_count),
                                rf = rich_factor,
                                col = go_type_colors[as.numeric(as.factor(go_ring_data$ONTOLOGY))])
  track_bed_ratio$rf <- track_bed_ratio$rf / max(track_bed_ratio$rf) * 9.5

  pdf("GO_Rich_Circlize_Enhanced.pdf", width = 10, height = 10)
  par(omi = c(0.1, 0.1, 0.1, 1.5))
  circos.clear()
  circos.par(track.margin = c(0.01, 0.01))
  circos.genomicInitialize(circ_df, plotType = "none")

  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.idx <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.idx, cex = 0.8,
                facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.09, bg.border = NA,
  bg.col = go_type_colors[as.numeric(as.factor(go_ring_data$ONTOLOGY))])

  for (si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.65, sector.index = si, track.index = 1,
                major.at = seq(0, max(bg_gene_count), by = 100), labels.facing = "clockwise")
  }

  circos.genomicTrack(track_bed_bg, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                           col = value[, "col_bg"], border = NA, ...)
                        circos.genomicText(region, value, y = 0.5, labels = value[, "n_bg"],
                                           cex = 0.7, adj = 0, ...)
                      })

  circos.genomicTrack(track_bed_sel, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                           col = '#6D72B9', border = NA, ...)
                        circos.genomicText(region, value, y = 0.6, labels = value[, "n_sel"],
                                           cex = 0.75, adj = 0, ...)
                      })

  circos.genomicTrack(track_bed_ratio, ylim = c(0, 10), track.height = 0.35,
                      bg.border = "white", bg.col = "#9DADBC",
                      panel.fun = function(region, value, ...) {
                        cell.xlim <- get.cell.meta.data("cell.xlim")
                        cell.ylim <- get.cell.meta.data("cell.ylim")
                        for (j in 1:9) {
                          y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                          circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                        }
                        circos.genomicRect(region, value, ytop = 0, ybottom = value[, "rf"],
                                           col = value[, "col"], border = NA, ...)
                      })
  circos.clear()

  legend_middle <- Legend(
    labels = c('Total number of genes', 'Number of selected genes', 'Rich Factor (0‑1)'),
    type = "points", pch = c(15, 15, 17),
    legend_gp = gpar(col = c('pink', '#BA55D3', go_type_colors[1])),
    title = "", nrow = 3, size = unit(3, "mm"))
  circle_sz <- unit(1, "snpc")
  draw(legend_middle, x = circle_sz * 0.42)

  legend_main <- Legend(
    labels = c("Biological Process", "Cellular Component", "Molecular Function"),
    type = "points", pch = 15,
    legend_gp = gpar(col = go_type_colors), title_position = "topcenter",
    title = "GO Category", nrow = 3, size = unit(3, "mm"),
    grid_height = unit(5, "mm"), grid_width = unit(5, "mm"))

  legend_pval <- Legend(
    labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
    type = "points", pch = 16, legend_gp = gpar(col = pval_palette),
    title = "-log10(Pvalue)", title_position = "topcenter",
    grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), size = unit(3, "mm"))

  legend_all <- packLegend(legend_main, legend_pval)
  draw(legend_all, x = circle_sz * 0.85, y = circle_sz * 0.55, just = "left")
  dev.off()
  cat("Enhanced GO circle plot saved to GO_Rich_Circlize_Enhanced.pdf\n")
}

cat("GO enrichment pipeline completed\n")
