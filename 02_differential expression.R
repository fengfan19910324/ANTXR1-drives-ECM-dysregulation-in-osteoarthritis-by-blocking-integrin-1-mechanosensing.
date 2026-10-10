# ========================== 0. Global parameters ==========================
# Please set your working directory before running
# setwd("your_local_work_directory")

threshold_logFC <- 0.585
threshold_adjP  <- 0.05
max_display_genes <- 50
result_dir <- "result_csv"
if(!dir.exists(result_dir)) dir.create(result_dir)

step1_start <- Sys.time()
cat("[Step 1] Start running...\n")

# Load required packages
library(limma)

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

step1_end <- Sys.time()
cat("[Step 1] Finished, elapsed:", round(difftime(step1_end, step1_start, units="secs"), 2), "sec\n\n")

# ========================== 2. File path configuration ==========================
step2_start <- Sys.time()
cat("[Step 2] Parameter and file‑path configuration...\n")

file_expr  <- "Sample Type Matrix.csv"

step2_end <- Sys.time()
cat("[Step 2] Finished, elapsed:", round(difftime(step2_end, step2_start, units="secs"), 2), "sec\n\n")

# ========================== Read expression matrix & ensure numeric ==========================
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")

expr_raw <- read.table(file_expr, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
rownames(expr_raw) <- expr_raw[,1]
expr_mat <- expr_raw[,-1, drop=FALSE]

expr_mat <- as.matrix(expr_mat)
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- rownames(expr_raw)
colnames(expr_mat) <- colnames(expr_raw)[-1]
combined_expr <- expr_mat

if(!is.numeric(expr_mat)) stop("Non‑numeric columns remain in expression matrix!")
cat("Expression matrix dimension:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# ========================== Auto‑detect sample groups ==========================
sample_names <- colnames(expr_mat)
group_info <- ifelse(grepl("_Con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_OA$", sample_names, ignore.case = TRUE), "Osteoarthritis", "Unknown"))

if(any(group_info == "Unknown")) stop("Unknown‑group samples detected, check sample name suffix!")

num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Osteoarthritis")
cat("Group info: Control =", num_ctrl, ", Osteoarthritis =", num_treat, "\n")

# ========================== limma differential expression analysis ==========================
group_labels <- factor(group_info, levels = c("Control", "Osteoarthritis"))
design_mat <- model.matrix(~0 + group_labels)
colnames(design_mat) <- c("Control", "Osteoarthritis")

fit <- lmFit(expr_mat, design_mat)
contrast_mat <- makeContrasts(Osteoarthritis - Control, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

all_diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
cat("DE analysis complete. Genes with FDR<0.05:", sum(all_diff_results$adj.P.Val < 0.05), "\n")

write.table(cbind(Gene=rownames(all_diff_results), all_diff_results),
            file="DE_results.csv", sep=",", quote=FALSE, row.names=FALSE)

# ------------------------- Step 6: Filter significant DEGs & export table -------------------------
step6_start <- Sys.time()
cat("[Step 6] Filter significant DEGs and export tables...\n")

significant_DEGs <- all_diff_results[with(all_diff_results,
                                          (abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP)), ]

output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)

SE <- ifelse(as.numeric(output_DEGs[, "t"]) != 0,
             abs(as.numeric(output_DEGs[, "logFC"]) / as.numeric(output_DEGs[, "t"])),
             NA)
output_DEGs <- cbind(output_DEGs, SE = SE)

output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors = FALSE)
desired_order <- c("Gene", "logFC", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]

write.table(output_DEGs, file = "DE_significant_genes.xls", sep = "\t", quote = FALSE, row.names = FALSE)

step6_end <- Sys.time()
cat("[Step 6] Finished, elapsed:", round(difftime(step6_end, step6_start, units = "secs"), 2), "sec\n\n")

# ========================== 7. DEG heatmap ==========================
step7_start <- Sys.time()
cat("[Step 7] Plot DEG heatmap\n")

ordered_DEGs <- significant_DEGs[order(as.numeric(as.vector(significant_DEGs$logFC))), ]
ordered_gene_names <- rownames(ordered_DEGs)
total_DEG_count <- length(ordered_gene_names)

if (total_DEG_count > (max_display_genes * 2)) {
  selected_gene_set <- ordered_gene_names[c(1:max_display_genes,
                                            (total_DEG_count - max_display_genes + 1):total_DEG_count)]
} else {
  selected_gene_set <- ordered_gene_names
}

heatmap_expr <- combined_expr[selected_gene_set, ]

sample_annotation <- data.frame(Group = factor(c(rep("Control", num_ctrl),
                                                 rep("Osteoarthritis", num_treat))))
rownames(sample_annotation) <- colnames(combined_expr)

annotation_colors <- list(
  Group = c("Control" = "#66C2A5",
            "Osteoarthritis"   = "#FC8D62")
)

color_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(255)

heatmap_result <- pheatmap(
  mat               = heatmap_expr,
  annotation_col    = sample_annotation,
  annotation_colors = annotation_colors,
  color             = color_palette,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  scale             = "row",
  fontsize          = 12,
  fontsize_row      = 7,
  fontsize_col      = 10,
  border_color      = NA,
  main              = "Differential Expression Heatmap",
  silent            = TRUE
)

control_count <- num_ctrl
treatment_count <- num_treat
sample_text <- paste("Control:", control_count, "| Osteoarthritis:", treatment_count)

library(grid)
pdf(file = "DE_heatmap.pdf", height = 8, width = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2,
                                           heights = unit(c(0.93, 0.07), "npc"))))
pushViewport(viewport(layout.pos.row = 1))
grid.draw(heatmap_result$gtable)
popViewport()
pushViewport(viewport(layout.pos.row = 2))
grid.text(label = sample_text,
          x = 0.5,
          y = 0.5,
          gp = gpar(fontsize = 13, fontface = "bold", col = "#333333"))
popViewport()
popViewport()
dev.off()

step7_end <- Sys.time()
cat("[Step 7] Finished, elapsed:", round(difftime(step7_end, step7_start, units="secs"), 2), "sec\n\n")

# ========================== 8. Volcano plot ==========================
step8_start <- Sys.time()
cat("[Step 8] Plot volcano plot...\n")

significance_status <- ifelse(
  (all_diff_results$adj.P.Val < threshold_adjP & abs(all_diff_results$logFC) > threshold_logFC),
  ifelse(all_diff_results$logFC > threshold_logFC, "Up regulated", "Down regulated"),
  "Not Significant"
)
all_diff_results$Significance <- significance_status

volcano_plot <- ggplot(all_diff_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Down regulated"  = "#1E90FF",
                                "Not Significant" = "#808080",
                                "Up regulated"    = "#FF4500")) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(threshold_adjP),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P‑value") +
  theme_minimal(base_size = 14) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, color = "#2F4F4F"),
        axis.title    = element_text(face = "bold", color = "#2F4F4F"),
        axis.text     = element_text(color = "#2F4F4F"),
        panel.grid.major = element_line(color = "#D3D3D3", linetype = "dotted"),
        panel.grid.minor = element_blank())

pdf(file = "DE_volcano.pdf", width = 6, height = 6)
print(volcano_plot)
dev.off()

step8_end <- Sys.time()
cat("[Step 8] Finished, elapsed:", round(difftime(step8_end, step8_start, units="secs"), 2), "sec\n\n")

# ========================== 9. PCA analysis ==========================
step9_start <- Sys.time()
cat("[Step 9] Running PCA analysis...\n")

pca_result <- prcomp(t(combined_expr), scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(combined_expr),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = factor(c(rep("Control", num_ctrl), rep("Osteoarthritis", num_treat)))
)

pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 1) +
  geom_point(size = 4, alpha = 0.8, shape = 16) +
  geom_text_repel(aes(label = Sample), size = 4, fontface = "bold", max.overlaps = 15) +
  scale_color_manual(values = c("Control" = "#0072B2",
                                "Osteoarthritis" = "#E69F00")) +
  labs(title = "PCA Analysis of Gene Expression",
       x = paste("PC1 (", pca_var_perc[1], "%)", sep = ""),
       y = paste("PC2 (", pca_var_perc[2], "%)", sep = ""),
       color = "") +
  theme_classic(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 18, color = "#333333"),
    axis.title    = element_text(face = "bold", size = 16, color = "#333333"),
    axis.text     = element_text(size = 14, color = "#555555"),
    legend.text   = element_text(size = 14, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_line(color = "#DDDDDD", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )

ggsave("DE_PCA.pdf", pca_plot, width = 7, height = 7, dpi = 300)

step9_end <- Sys.time()
cat("[Step 9] PCA finished, elapsed:", round(difftime(step9_end, step9_start, units="secs"), 2), "sec\n\n")

# ========================== Volcano plot labeling top genes ==========================
de_table <- all_diff_results
sig_degs <- de_table[with(de_table, abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP), ]

up_top <- rownames(sig_degs)[order(-sig_degs$logFC)][1:25]
down_top <- rownames(sig_degs)[order(sig_degs$logFC)][1:25]
label_genes <- unique(c(up_top, down_top))

de_table$Gene <- rownames(de_table)
de_table$label <- ifelse(de_table$Gene %in% label_genes, de_table$Gene, "")

de_table$Group <- "Not Significant"
de_table$Group[de_table$logFC > threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Up regulated"
de_table$Group[de_table$logFC < -threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Down regulated"

group_counts <- table(de_table$Group)
up_count   <- ifelse("Up regulated"   %in% names(group_counts), as.integer(group_counts["Up regulated"]), 0)
down_count <- ifelse("Down regulated" %in% names(group_counts), as.integer(group_counts["Down regulated"]), 0)
not_count  <- ifelse("Not Significant"%in% names(group_counts), as.integer(group_counts["Not Significant"]), 0)

legend_labels <- c(
  paste0("Up regulated (", up_count, ")"),
  paste0("Down regulated (", down_count, ")"),
  paste0("Not Significant (", not_count, ")")
)

volcano_labeled_plot <- ggplot(de_table, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(aes(label=label),
                  size=3,
                  max.overlaps=50,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color="grey50",
                  show.legend = FALSE) +
  scale_color_manual(
    values = c("Up regulated" = "#FF4500",
               "Down regulated" = "#1E90FF",
               "Not Significant" = "#808080"),
    breaks = c("Up regulated", "Down regulated", "Not Significant"),
    labels = legend_labels
  ) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(threshold_adjP), linetype="dashed", color="black") +
  labs(title="Volcano Plot with Top Genes Labeled",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P‑value") +
  theme_minimal(base_size = 14) +
  theme(plot.title=element_text(face="bold", hjust=0.5, color="#2F4F4F"))

ggsave("DE_volcano_with_labels.pdf", volcano_labeled_plot, width=8, height=8, dpi=300)

# Export DEG gene list
diff_gene_list <- data.frame(gene = rownames(significant_DEGs))
write.table(diff_gene_list,
            file = "DEG_geneList.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

cat("All steps completed successfully!\n")
