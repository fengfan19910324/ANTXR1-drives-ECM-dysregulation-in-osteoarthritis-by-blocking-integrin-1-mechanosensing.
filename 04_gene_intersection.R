
# ---- Packages ----
library(ggvenn)
library(RColorBrewer)
library(viridis)
library(ComplexUpset)
library(ggplot2)

# ---- 1. Working directory (edit to your local path) ----
setwd("D:\\")

# ---- 2. Output folder ----
output_folder <- "output_folder"
if (!dir.exists(output_folder)) dir.create(output_folder)

# ---- 3. Discover all .txt files ----
txt_files <- list.files(pattern = "\\.txt$", full.names = TRUE)
if (length(txt_files) == 0) {
  stop("No .txt files found in the working directory. Check path / file suffix.")
}

# ---- 4. Read each file into a named list ----
gene_list <- list()
for (file in txt_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  rt <- read.table(file, header = FALSE, sep = "\t", check.names = FALSE)
  genes_vec <- unique(as.vector(rt[, 1]))     # de-duplicate within file
  gene_list[[file_name]] <- genes_vec
  cat(file_name, "genes:", length(gene_list[[file_name]]), "\n")
}

# ---- 5. Colors ----
nColors <- length(gene_list)
myColors <- viridis(nColors)

# ---- 6. Venn diagram ----
# Note: ggvenn is designed for 2-4 sets. Beyond 4 sets the diagram
# becomes unreadable; use the UpSet plot instead.
pdf(file = file.path(output_folder, "venn.pdf"), width = 10, height = 10)
ggvenn(gene_list,
       show_percentage = TRUE,
       stroke_color = "white",
       stroke_size = 1.5,
       fill_color = myColors,
       set_name_color = "black",
       set_name_size = 8,
       text_size = 6,
       text_color = "black")
dev.off()
cat("Venn diagram saved to:", file.path(output_folder, "venn.pdf"), "\n")

# ---- 7. Global intersection across ALL lists ----
intersect_genes <- Reduce(intersect, gene_list)
cat("Genes shared across ALL lists:", length(intersect_genes), "\n")

global_intersect_df <- data.frame(Gene = intersect_genes,
                                  stringsAsFactors = FALSE)
write.csv(global_intersect_df,
          file = file.path(output_folder, "IntersectionGenes.csv"),
          row.names = FALSE, quote = FALSE)
cat("Global intersection saved to:",
    file.path(output_folder, "IntersectionGenes.csv"), "\n")

# ---- 8. Pairwise intersection table ----
pairwise_intersections <- data.frame(
  File1 = character(),
  File2 = character(),
  Intersection_Count = integer(),
  Intersection_Genes = character(),
  stringsAsFactors = FALSE
)
files <- names(gene_list)
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    file1 <- files[i]
    file2 <- files[j]
    common_genes <- intersect(gene_list[[file1]], gene_list[[file2]])
    pairwise_intersections <- rbind(
      pairwise_intersections,
      data.frame(
        File1 = file1,
        File2 = file2,
        Intersection_Count = length(common_genes),
        Intersection_Genes = paste(common_genes, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
  }
}
write.csv(pairwise_intersections,
          file = file.path(output_folder, "PairwiseIntersectionGenes.csv"),
          row.names = FALSE, quote = FALSE)
cat("Pairwise intersection table saved to:",
    file.path(output_folder, "PairwiseIntersectionGenes.csv"), "\n")

# ---- 9. Build membership table for UpSet ----
all_genes <- unique(unlist(gene_list))
membership_df <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
for (set_name in names(gene_list)) {
  membership_df[[set_name]] <- as.integer(membership_df$Gene %in% gene_list[[set_name]])
}
head(membership_df)

# ---- 10. UpSet plot ----
upset_plot <- upset(
  membership_df,
  intersect = names(gene_list),
  sort_intersections_by = "degree",
  base_annotations = list(
    "Intersection Size" = intersection_size(text = list(size = 12))
  )
)
print(upset_plot)

ggsave(
  filename = file.path(output_folder, "ComplexUpSetPlot.pdf"),
  plot = upset_plot,
  width = 10,
  height = 6,
  device = "pdf"
)
cat("UpSet plot saved.\n")
