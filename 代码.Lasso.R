# ---- Load packages ----
library(glmnet)
library(ggplot2)
library(caret)
library(dplyr)
library(pROC)

# ---- Optional: Excel export (auto-install on first run) ----
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)

# ---- Working directory (edit to your local path) ----
setwd("D:\\")

# Fix random seed for reproducibility
set.seed(123)

# ------------------------ Data preparation ---------------------------
cat("Reading expression matrix...\n")
expression_matrix <- read.csv("Sample Type Matrix.csv", row.names = 1, check.names = FALSE)
cat("Expression matrix dimension:", dim(expression_matrix), "\n")

cat("Reading sample lists...\n")
control_samples <- readLines("Control.txt")
oa_samples <- readLines("OA.txt")
cat("Number of Control samples:", length(control_samples), "\n")
cat("Number of OA samples:", length(oa_samples), "\n")

# Keep only samples present in the expression matrix
control_samples <- control_samples[control_samples %in% colnames(expression_matrix)]
oa_samples <- oa_samples[oa_samples %in% colnames(expression_matrix)]
cat("Control samples found in matrix:", length(control_samples), "\n")
cat("OA samples found in matrix:", length(oa_samples), "\n")

if (length(control_samples) == 0 | length(oa_samples) == 0) {
  stop("No matching samples found in the expression matrix.")
}

# Response variable: 0 = Control, 1 = OA
all_samples <- c(control_samples, oa_samples)
y <- c(rep(0, length(control_samples)), rep(1, length(oa_samples)))

# Extract expression data for all samples
X_all <- t(expression_matrix[, all_samples])

# ------------------------ Load & subset overlap genes ---------------------------
cat("Reading overlap genes...\n")
intersection_genes <- read.csv("IntersectionGenes.csv", header = TRUE)

# Locate the gene column
if ("Gene" %in% colnames(intersection_genes)) {
  overlap_genes <- intersection_genes$Gene
} else if ("GeneSymbol" %in% colnames(intersection_genes)) {
  overlap_genes <- intersection_genes$GeneSymbol
} else {
  overlap_genes <- intersection_genes[, 1]
}
cat("Number of overlap genes:", length(overlap_genes), "\n")

# Keep only genes present in the matrix
overlap_genes <- overlap_genes[overlap_genes %in% rownames(expression_matrix)]
cat("Overlap genes found in matrix:", length(overlap_genes), "\n")
if (length(overlap_genes) == 0) {
  stop("No overlap genes found in the expression matrix.")
}

X <- X_all[, overlap_genes]
cat("Final data dimension - samples:", nrow(X), "genes:", ncol(X), "\n")
cat("Class distribution - Control:", sum(y == 0), "OA:", sum(y == 1), "\n")

# ------------------------ Data preprocessing ---------------------------
# Remove zero-variance genes
zero_var_genes <- apply(X, 2, var) == 0
if (any(zero_var_genes)) {
  cat("Removing zero-variance genes:", sum(zero_var_genes), "\n")
  X <- X[, !zero_var_genes]
}

# Standardize expression per gene
X_scaled <- scale(X)

# ------------------------ LASSO regression ---------------------------
cat("Running LASSO regression...\n")

# Cross-validated lambda selection (10-fold, AUC metric)
cv_lasso <- cv.glmnet(X_scaled, y,
                      family = "binomial",
                      alpha = 1,                       # LASSO
                      nfolds = min(10, nrow(X)),       # 10-fold CV
                      type.measure = "auc",            # evaluate by AUC
                      standardize = FALSE,             # already standardized
                      parallel = FALSE)                # set TRUE to parallelize (needs doParallel)

lambda_min <- cv_lasso$lambda.min
lambda_1se <- cv_lasso$lambda.1se
cat("Best lambda (min):", lambda_min, "\n")
cat("Best lambda (1se):", lambda_1se, "\n")
cat("Maximum AUC:", max(cv_lasso$cvm), "\n")

# ------------------------ Gene selection (objective) ---------------------------
cat("\nSelecting genes under both lambda criteria:\n")

# 1. lambda.1se (more conservative, fewer genes)
lasso_model_1se <- glmnet(X_scaled, y,
                          family = "binomial", alpha = 1,
                          lambda = lambda_1se, standardize = FALSE)
coefficients_1se <- as.matrix(coef(lasso_model_1se))
selected_genes_1se <- coefficients_1se[coefficients_1se != 0, ]
selected_genes_1se <- selected_genes_1se[names(selected_genes_1se) != "(Intercept)"]
cat("Genes selected by lambda.1se:", length(selected_genes_1se), "\n")

# 2. lambda.min (better prediction, more genes)
lasso_model_min <- glmnet(X_scaled, y,
                          family = "binomial", alpha = 1,
                          lambda = lambda_min, standardize = FALSE)
coefficients_min <- as.matrix(coef(lasso_model_min))
selected_genes_min <- coefficients_min[coefficients_min != 0, ]
selected_genes_min <- selected_genes_min[names(selected_genes_min) != "(Intercept)"]
cat("Genes selected by lambda.min:", length(selected_genes_min), "\n")

# Default to lambda.1se result as the primary outcome (more robust)
final_genes <- selected_genes_1se
cat("Final selected gene count:", length(final_genes), "\n")

# ------------------------ Visualization ---------------------------
pdf("LASSO_Cross_Validation.pdf", width = 8, height = 6)
plot(cv_lasso)
title("LASSO cross-validation", line = 2.5)
dev.off()

# Export CV curve data (matches LASSO_Cross_Validation plot)
cv_df <- data.frame(
  Lambda = cv_lasso$lambda,
  Mean_Metric = cv_lasso$cvm
)
cv_df$Lambda_Type <- "Other"
cv_df$Lambda_Type[cv_df$Lambda == lambda_min] <- "lambda.min"
cv_df$Lambda_Type[cv_df$Lambda == lambda_1se] <- "lambda.1se"
write.csv(cv_df, "LASSO_CrossValidation_Data.csv", row.names = FALSE)
write_xlsx(cv_df, "LASSO_Cross_Validation_01.xlsx")

# Coefficient path plot
pdf("LASSO_Coefficient_Path.pdf", width = 10, height = 6)
plot(cv_lasso$glmnet.fit, xvar = "lambda", label = TRUE)
abline(v = log(lambda_1se), col = "red", lty = 2)
abline(v = log(lambda_min), col = "blue", lty = 2)
legend("topright", legend = c("lambda.1se", "lambda.min"),
       col = c("red", "blue"), lty = 2)
title("LASSO coefficient path")
dev.off()

# ------------------------ Model evaluation ---------------------------
if (length(selected_genes_1se) > 0) {
  X_1se <- X_scaled[, names(selected_genes_1se), drop = FALSE]
  model_1se <- glm(y ~ ., data = data.frame(X_1se), family = binomial)
  predictions_1se <- predict(model_1se, type = "response")
  roc_1se <- roc(y, predictions_1se)
  auc_1se <- auc(roc_1se)
} else {
  auc_1se <- NA
}

if (length(selected_genes_min) > 0) {
  X_min <- X_scaled[, names(selected_genes_min), drop = FALSE]
  model_min <- glm(y ~ ., data = data.frame(X_min), family = binomial)
  predictions_min <- predict(model_min, type = "response")
  roc_min <- roc(y, predictions_min)
  auc_min <- auc(roc_min)
} else {
  auc_min <- NA
}

# ROC comparison
pdf("ROC_Curve_Comparison.pdf", width = 8, height = 6)
if (!is.na(auc_1se) & !is.na(auc_min)) {
  plot(roc_1se, col = "blue", main = "ROC comparison across lambda criteria")
  lines(roc_min, col = "red")
  legend("bottomright",
         legend = c(paste("lambda.1se (AUC =", round(auc_1se, 3), ")", length(selected_genes_1se), "genes"),
                    paste("lambda.min (AUC =", round(auc_min, 3), ")", length(selected_genes_min), "genes")),
         col = c("blue", "red"), lwd = 2)
} else if (!is.na(auc_1se)) {
  plot(roc_1se, col = "blue", main = "ROC curve (lambda.1se)")
  legend("bottomright", legend = paste("AUC =", round(auc_1se, 3)), col = "blue", lwd = 2)
} else if (!is.na(auc_min)) {
  plot(roc_min, col = "red", main = "ROC curve (lambda.min)")
  legend("bottomright", legend = paste("AUC =", round(auc_min, 3)), col = "red", lwd = 2)
}
dev.off()

# ------------------------ Result export ---------------------------
results_list <- list()

if (length(selected_genes_1se) > 0) {
  result_1se <- data.frame(
    Gene = names(selected_genes_1se),
    Coefficient = selected_genes_1se,
    Abs_Coefficient = abs(selected_genes_1se),
    Lambda_Type = "lambda.1se",
    AUC = round(auc_1se, 4)
  )
  result_1se <- result_1se[order(-result_1se$Abs_Coefficient), ]
  results_list[["lambda_1se"]] <- result_1se
}

if (length(selected_genes_min) > 0) {
  result_min <- data.frame(
    Gene = names(selected_genes_min),
    Coefficient = selected_genes_min,
    Abs_Coefficient = abs(selected_genes_min),
    Lambda_Type = "lambda.min",
    AUC = round(auc_min, 4)
  )
  result_min <- result_min[order(-result_min$Abs_Coefficient), ]
  results_list[["lambda_min"]] <- result_min
}

if (length(results_list) > 0) {
  all_results <- do.call(rbind, results_list)
  write.csv(all_results, "LASSO_All_Results.csv", row.names = FALSE)
  write_xlsx(all_results, "LASSO_All_Results.xlsx")

  for (lambda_type in names(results_list)) {
    filename_csv <- paste0("LASSO_Results_", lambda_type, ".csv")
    filename_xlsx <- paste0("LASSO_Results_", lambda_type, ".xlsx")
    write.csv(results_list[[lambda_type]], filename_csv, row.names = FALSE)
    write_xlsx(results_list[[lambda_type]], filename_xlsx)
  }

  if (length(selected_genes_1se) > 0) {
    writeLines(names(selected_genes_1se), "Recommended_Genes_lambda_1se.txt")
  }
  if (length(selected_genes_min) > 0) {
    writeLines(names(selected_genes_min), "Alternative_Genes_lambda_min.txt")
  }
} else {
  empty_result <- data.frame(
    Gene = character(0),
    Coefficient = numeric(0),
    Abs_Coefficient = numeric(0),
    Lambda_Type = character(0),
    AUC = numeric(0)
  )
  write.csv(empty_result, "LASSO_All_Results.csv", row.names = FALSE)
  write_xlsx(empty_result, "LASSO_All_Results.xlsx")
}

# ------------------------ Stability analysis ---------------------------
cat("Running stability analysis...\n")
stability_results <- list()
seeds <- c(123, 456, 789, 321, 654)

for (seed in seeds) {
  set.seed(seed)
  cv_temp <- cv.glmnet(X_scaled, y,
                       family = "binomial", alpha = 1,
                       nfolds = min(10, nrow(X)),
                       type.measure = "auc", standardize = FALSE)
  lasso_temp <- glmnet(X_scaled, y,
                       family = "binomial", alpha = 1,
                       lambda = cv_temp$lambda.1se, standardize = FALSE)
  coef_temp <- as.matrix(coef(lasso_temp))
  genes_temp <- coef_temp[coef_temp != 0, ]
  genes_temp <- genes_temp[names(genes_temp) != "(Intercept)"]
  stability_results[[as.character(seed)]] <- names(genes_temp)
}

# Gene selection frequency across seeds
all_genes <- unique(unlist(stability_results))
gene_frequency <- sapply(all_genes, function(gene) {
  sum(sapply(stability_results, function(x) gene %in% x))
})

stability_df <- data.frame(
  Gene = all_genes,
  Frequency = gene_frequency,
  Percentage = round(gene_frequency / length(seeds) * 100, 1)
)
stability_df <- stability_df[order(-stability_df$Frequency), ]
write.csv(stability_df, "Gene_Selection_Stability.csv", row.names = FALSE)
write_xlsx(stability_df, "Gene_Selection_Stability.xlsx")
cat("Stability analysis completed.\n")

# ------------------------ Console summary ---------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("            LASSO logistic regression - result summary\n")
cat(strrep("=", 60), "\n", sep = "")
cat("Input gene count:", length(overlap_genes), "\n")
cat("Genes selected by lambda.1se:", length(selected_genes_1se), "\n")
cat("Genes selected by lambda.min:", length(selected_genes_min), "\n")
if (length(selected_genes_1se) > 0) cat("AUC (lambda.1se):", round(auc_1se, 3), "\n")
if (length(selected_genes_min) > 0) cat("AUC (lambda.min):", round(auc_min, 3), "\n")

cat("\nRecommended: use lambda.1se result (more robust, fewer genes)\n")
if (length(selected_genes_1se) > 0) {
  cat("Genes selected by lambda.1se:\n")
  for (i in 1:min(10, length(selected_genes_1se))) {
    cat(i, ". ", names(selected_genes_1se)[i], " (coef: ",
        round(selected_genes_1se[i], 4), ")\n", sep = "")
  }
}

cat("\nStability analysis:\n")
if (nrow(stability_df) > 0) {
  stable_genes <- stability_df[stability_df$Frequency >= 3, ]
  cat("Genes stable across", length(seeds), "repeats (>=60%):", nrow(stable_genes), "\n")
  if (nrow(stable_genes) > 0) {
    for (i in 1:min(5, nrow(stable_genes))) {
      cat("  - ", stable_genes$Gene[i], " (", stable_genes$Percentage[i], "%)\n", sep = "")
    }
  }
}

cat("\nOutput files:\n")
cat("  - LASSO_Cross_Validation.pdf: cross-validation plot\n")
cat("  - LASSO_Coefficient_Path.pdf: coefficient path plot\n")
cat("  - ROC_Curve_Comparison.pdf: ROC comparison\n")
cat("  - LASSO_All_Results.csv/xlsx: combined results\n")
cat("  - LASSO_Results_lambda_1se.csv/xlsx: lambda.1se details\n")
cat("  - LASSO_Results_lambda_min.csv/xlsx: lambda.min details\n")
cat("  - Recommended_Genes_lambda_1se.txt: recommended genes\n")
cat("  - Alternative_Genes_lambda_min.txt: alternative genes\n")
cat("  - Gene_Selection_Stability.csv/xlsx: stability analysis\n")
cat(strrep("=", 60), "\n", sep = "")

# ------------------------ Feature importance plots ---------------------------
if (length(final_genes) > 0) {
  if (length(selected_genes_1se) > 0) {
    importance_plot_1se <- ggplot(result_1se,
                                  aes(x = reorder(Gene, Abs_Coefficient), y = Abs_Coefficient)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(title = "Gene importance (lambda.1se)",
           x = "Gene", y = "Absolute coefficient") +
      theme_minimal()
    ggsave("Gene_Importance_lambda_1se.pdf", importance_plot_1se, width = 8, height = 6)
    write_xlsx(result_1se, "Gene_Importance_lambda_1se.xlsx")
  }

  if (length(selected_genes_min) > 0) {
    importance_plot_min <- ggplot(result_min,
                                  aes(x = reorder(Gene, Abs_Coefficient), y = Abs_Coefficient)) +
      geom_col(fill = "darkred") +
      coord_flip() +
      labs(title = "Gene importance (lambda.min)",
           x = "Gene", y = "Absolute coefficient") +
      theme_minimal()
    ggsave("Gene_Importance_lambda_min.pdf", importance_plot_min, width = 8, height = 6)
    write_xlsx(result_min, "Gene_Importance_lambda_min_01.xlsx")
  }
}

# ------------------------ Save processed data ---------------------------
save(X, y, overlap_genes, selected_genes_1se, selected_genes_min, stability_df,
     file = "LASSO_Analysis_Data_Objective.RData")

cat("\nAnalysis completed. Choose the gene set based on your goal:\n")
cat("  - Robust, interpretable model: lambda.1se result\n")
cat("  - Highest prediction accuracy: lambda.min result\n")
cat("  - Most stable genes: stability analysis\n")
cat("  - All key results exported to Excel as well\n")
