# -------------------- Fixed paths (do NOT modify upstream project paths) --------------------
ROOT_DIR <- "D:\\"
MODEL_DIR <- file.path(ROOT_DIR, "14.Machine‑Learning‑Diagnostic‑Model")
WORK_DIR <- file.path(ROOT_DIR, "23.SHAP‑Model‑Interpretation")

TRAIN_FILE <- file.path(MODEL_DIR, "train.csv")
TEST_FILE <- file.path(MODEL_DIR, "test.csv")

PROJECT21_FILE <- file.path(
  ROOT_DIR, "21.ROC‑External‑Validation‑Dataset2",
  "Sample Type Matrix.csv"
)
PROJECT22_FILE <- file.path(
  ROOT_DIR, "22.ROC‑External‑Validation‑Dataset3",
  "Sample Type Matrix.csv"
)

RESULT_DIR <- file.path(WORK_DIR, "SHAP_Result")
PLOT_DIR   <- file.path(RESULT_DIR, "Figures")
TABLE_DIR  <- file.path(RESULT_DIR, "Tables")
RDS_DIR    <- file.path(RESULT_DIR, "RDS")

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# Frozen six‑gene signature
FINAL_GENES <- c(
  "TBC1D8B", "FAM65B", "B4GALNT1",
  "ITGBL1", "ANTXR1", "CNIH3"
)

# Keep NA normally: read final algorithm determined by nested‑CV in Project14.
# Override ONLY when Project14 lacks model.selectionRule.txt / Selected_Final_Model.
FINAL_METHOD_OVERRIDE <- NA_character_

# TRUE: generate extra Project22 sensitivity plots. Project22 is NOT treated as independent validation.
RUN_PROJECT22_SENSITIVITY <- TRUE

# Exact SHAP for six‑gene model. Max background samples from training set.
BACKGROUND_MAX_N <- 50L
STRICT_PREDICTION_AUDIT <- TRUE

# -------------------- Dependencies --------------------
required_packages <- c(
  "kernelshap", "shapviz", "ggplot2", "patchwork",
  "data.table", "pheatmap", "RColorBrewer", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "), "\n",
    "Please run dependency installation script: ", file.path(WORK_DIR, "00.Install_SHAP_Dependencies.R")
  )
}

suppressPackageStartupMessages({
  library(kernelshap)
  library(shapviz)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# Global SHAP‑oriented theme
theme_shap <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E7E7E7", linewidth = 0.35),
      axis.text = element_text(color = "#222222"),
      axis.title = element_text(color = "#222222", face = "bold"),
      plot.title = element_text(face = "bold", color = "#16324F"),
      plot.subtitle = element_text(color = "#4D4D4D"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#EAF2F8", color = NA),
      strip.text = element_text(face = "bold", color = "#16324F")
    )
}
theme_set(theme_shap())

# -------------------- Utility functions --------------------
assert_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
}

canonical_gsm <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("GSM[0-9]+", x, ignore.case = TRUE),
                sub(".*(GSM[0-9]+).*", "\\1", x, ignore.case = TRUE), x)
  toupper(out)
}

parse_type <- function(x) {
  z <- tolower(as.character(x))
  ans <- rep(NA_integer_, length(z))
  ans[grepl("(^|[_.-])(con|control|normal)([_.-]|$)", z)] <- 0L
  ans[grepl("(^|[_.-])(oa|tra|treat|case|disease)([_.-]|$)", z)] <- 1L
  ans
}

as_binary_type <- function(x, sample_names = NULL) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    y <- rep(NA_integer_, length(x))
    lx <- tolower(trimws(x))
    y[lx %in% c("0", "con", "control", "normal")] <- 0L
    y[lx %in% c("1", "oa", "case", "disease", "treat", "tra")] <- 1L
  } else {
    y <- as.integer(x)
  }
  if (anyNA(y) && !is.null(sample_names)) {
    from_name <- parse_type(sample_names)
    y[is.na(y)] <- from_name[is.na(y)]
  }
  if (anyNA(y) || !all(y %in% c(0L, 1L))) {
    stop("Unrecognized Control/OA labels detected.")
  }
  y
}

read_sample_by_gene <- function(path, default_cohort) {
  assert_file(path, "sample‑by‑gene matrix")
  dat <- read.csv(path, check.names = FALSE, row.names = 1)
  missing_genes <- setdiff(FINAL_GENES, colnames(dat))
  if (length(missing_genes)) {
    stop(basename(path), " missing genes: ", paste(missing_genes, collapse = ", "))
  }
  if (!"Type" %in% colnames(dat)) stop(basename(path), " lacks 'Type' column.")
  ids <- rownames(dat)
  cohort <- rep(default_cohort, length(ids))
  cohort[grepl("GSE114007", ids, ignore.case = TRUE)] <- "GSE114007_cartilage"
  cohort[grepl("GSE63359",  ids, ignore.case = TRUE)] <- "GSE63359_blood"

  x <- as.matrix(dat[, FINAL_GENES, drop = FALSE])
  storage.mode(x) <- "double"
  meta <- data.frame(
    Sample = ids,
    Cohort = cohort,
    Type = as_binary_type(dat$Type, ids),
    Group = factor(as_binary_type(dat$Type, ids), levels = c(0, 1),
                   labels = c("Control", "OA")),
    stringsAsFactors = FALSE,
    row.names = ids
  )
  list(x = x, meta = meta)
}

read_gene_by_sample <- function(path, cohort_name) {
  assert_file(path, "gene‑by‑sample matrix")
  dat <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  gene_col <- colnames(dat)[1]
  genes <- as.character(dat[[gene_col]])
  idx <- match(FINAL_GENES, genes)
  if (anyNA(idx)) {
    stop(basename(path), " missing genes: ",
         paste(FINAL_GENES[is.na(idx)], collapse = ", "))
  }
  mat <- as.matrix(dat[idx, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- FINAL_GENES
  x <- t(mat)
  ids <- rownames(x)
  y <- parse_type(ids)
  if (anyNA(y)) stop(cohort_name, ": cannot infer group from sample names.")
  meta <- data.frame(
    Sample = ids,
    Cohort = cohort_name,
    Type = y,
    Group = factor(y, levels = c(0, 1), labels = c("Control", "OA")),
    stringsAsFactors = FALSE,
    row.names = ids
  )
  list(x = x[, FINAL_GENES, drop = FALSE], meta = meta)
}

fit_preprocessor <- function(x) {
  med <- apply(x, 2, median, na.rm = TRUE)
  xi <- x
  for (j in seq_len(ncol(xi))) xi[!is.finite(xi[, j]), j] <- med[j]
  mu    <- colMeans(xi)
  sigma <- apply(xi, 2, sd)
  sigma[!is.finite(sigma) | sigma == 0] <- 1
  list(median = med, mean = mu, sd = sigma, genes = colnames(x))
}

apply_preprocessor <- function(x, pp) {
  needed <- as.character(pp$genes)
  missing_genes <- setdiff(needed, colnames(x))
  if (length(missing_genes)) {
    stop("Input data lacks pre‑processing genes: ",
         paste(missing_genes, collapse = ", "))
  }
  z <- x[, needed, drop = FALSE]
  for (j in seq_len(ncol(z))) {
    bad <- !is.finite(z[, j])
    if (any(bad)) z[bad, j] <- pp$median[j]
  }
  z <- sweep(z, 2, pp$mean, "-")
  z <- sweep(z, 2, pp$sd, "/")
  z
}

read_selected_method <- function(model_dir, available_methods = NULL) {
  if (!is.na(FINAL_METHOD_OVERRIDE) && nzchar(FINAL_METHOD_OVERRIDE)) {
    return(FINAL_METHOD_OVERRIDE)
  }
  rule_file <- file.path(model_dir, "model.selectionRule.txt")
  if (file.exists(rule_file)) {
    txt <- readLines(rule_file, warn = FALSE)
    hit <- grep("^Final selected model:", txt, value = TRUE)
    if (length(hit)) return(trimws(sub("^Final selected model:", "", hit[1])))
  }
  report_file <- file.path(model_dir, "model.performance_detailed_report.txt")
  if (file.exists(report_file)) {
    report <- read.delim(report_file, check.names = FALSE)
    if ("Selected_Final_Model" %in% colnames(report)) {
      flag <- tolower(as.character(report$Selected_Final_Model)) %in%
        c("true", "t", "1", "yes")
      if (sum(flag, na.rm = TRUE) == 1L) return(as.character(report$Method[flag]))
    }
  }
  inner_file <- file.path(model_dir, "model.innerCV_method_comparison.txt")
  if (file.exists(inner_file)) {
    inner <- read.delim(inner_file, check.names = FALSE)
    if ("Method" %in% colnames(inner) && nrow(inner)) {
      return(as.character(inner$Method[1]))
    }
  }
  if (!is.null(available_methods)) {
    stop(
      "Cannot auto‑detect final model from Project14. Set FINAL_METHOD_OVERRIDE.\n",
      "Available candidates: ", paste(available_methods, collapse = ", ")
    )
  }
  NA_character_
}

load_frozen_model <- function(model_dir) {
  selected_file <- file.path(model_dir, "model.selectedModel.rds")
  all_file     <- file.path(model_dir, "model.MLmodel.rds")
  if (file.exists(selected_file)) {
    fit <- readRDS(selected_file)
    method <- read_selected_method(model_dir)
    if (is.na(method)) method <- class(fit)[1]
    return(list(fit = fit, method = method, source = selected_file))
  }
  assert_file(all_file, "Project14 model file")
  models <- readRDS(all_file)
  if (!is.list(models) || is.null(names(models))) {
    stop("model.MLmodel.rds is not a named list of models.")
  }
  method <- read_selected_method(model_dir, names(models))
  if (!method %in% names(models)) {
    stop("Detected final algorithm not present in model.MLmodel.rds: ", method)
  }
  list(fit = models[[method]], method = method, source = all_file)
}

load_required_model_namespace <- function(fit) {
  class_to_package <- c(
    lognet = "glmnet", "svm.formula" = "e1071", train = "caret",
    glmboost = "mboost", plsRglmmodel = "plsRglm",
    rfsrc = "randomForestSRC", gbm = "gbm",
    "xgb.Booster" = "xgboost", naiveBayes = "e1071"
  )
  cls <- class(fit)[1]
  pkg <- unname(class_to_package[cls])
  if (!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
    stop("Model class is ", cls, ", required prediction package missing: ", pkg)
  }
  invisible(pkg)
}

predict_oa_probability <- function(object, newdata) {
  z <- as.data.frame(newdata, check.names = FALSE)
  vars <- object$subFeature
  if (is.null(vars) || !length(vars)) {
    stop("Frozen model lacks subFeature; cannot guarantee predictor order.")
  }
  vars <- as.character(vars)
  missing_vars <- setdiff(vars, colnames(z))
  if (length(missing_vars)) {
    stop("Prediction data missing model predictors: ", paste(missing_vars, collapse = ", "))
  }
  z <- z[, vars, drop = FALSE]
  cls <- class(object)[1]
  out <- switch(
    cls,
    lognet = predict(object, type = "response", as.matrix(z)),
    glm    = predict(object, type = "response", newdata = z),
    "svm.formula" = {
      pred <- predict(object, z, probability = TRUE)
      prob <- attr(pred, "probabilities")
      if (is.null(prob)) stop("SVM did not return probabilities.")
      positive <- if ("1" %in% colnames(prob)) "1" else colnames(prob)[ncol(prob)]
      prob[, positive]
    },
    train = {
      prob <- predict(object, z, type = "prob")
      positive <- if ("1" %in% colnames(prob)) "1" else colnames(prob)[ncol(prob)]
      prob[[positive]]
    },
    glmboost     = predict(object, newdata = z, type = "response"),
    plsRglmmodel = predict(object, newdata = z, type = "response"),
    rfsrc = {
      pred <- predict(object, newdata = z)$predicted
      positive <- if ("1" %in% colnames(pred)) "1" else colnames(pred)[ncol(pred)]
      pred[, positive]
    },
    gbm          = predict(object, newdata = z, type = "response"),
    "xgb.Booster"= predict(object, as.matrix(z)),
    naiveBayes = {
      prob <- predict(object, newdata = z, type = "raw")
      positive <- if ("1" %in% colnames(prob)) "1" else colnames(prob)[ncol(prob)]
      prob[, positive]
    },
    stop("Unsupported model class: ", cls)
  )
  out <- as.numeric(out)
  if (length(out) != nrow(z) || any(!is.finite(out))) {
    stop("Frozen model returns invalid prediction probabilities.")
  }
  if (any(out < -1e‑8 | out > 1 + 1e‑8)) {
    stop("Predictions out of 0‑1 range; output may not be probability.")
  }
  pmin(pmax(out, 0), 1)
}

prediction_audit <- function(pred, ids, method, model_dir) {
  risk_file <- file.path(model_dir, "model.riskMatrix.txt")
  if (!file.exists(risk_file)) {
    warning("Project14 model.riskMatrix.txt missing; skip prediction‑reproduction audit.")
    return(invisible(NULL))
  }
  risk <- read.delim(risk_file, check.names = FALSE)
  if (!all(c("id", method) %in% colnames(risk))) {
    warning("riskMatrix does not contain final‑model column; skip audit: ", method)
    return(invisible(NULL))
  }
  mi <- match(canonical_gsm(ids), canonical_gsm(risk$id))
  ok <- !is.na(mi)
  if (sum(ok) < 5L) {
    warning("Too few matched training samples; skip prediction‑reproduction audit.")
    return(invisible(NULL))
  }
  old <- as.numeric(risk[mi[ok], method])
  new <- pred[ok]
  audit <- data.frame(
    N_matched = sum(ok),
    Correlation = suppressWarnings(cor(old, new, use = "complete.obs")),
    Mean_absolute_difference = mean(abs(old - new), na.rm = TRUE),
    Max_absolute_difference = max(abs(old - new), na.rm = TRUE)
  )
  write.csv(audit, file.path(TABLE_DIR, "Prediction_reproduction_audit.csv"),
            row.names = FALSE)
  if (STRICT_PREDICTION_AUDIT &&
      (is.na(audit$Correlation) || audit$Correlation < 0.999 ||
       audit$Max_absolute_difference > 1e‑4)) {
    stop(
      "Recomputed probabilities differ from frozen Project14 results. SHAP halted to avoid mis‑interpretation.\n",
      "Ensure model.selectedModel.rds, model.preprocessing.rds and riskMatrix come from identical nested‑CV run.\n",
      "Audit table: ", file.path(TABLE_DIR, "Prediction_reproduction_audit.csv")
    )
  }
  invisible(audit)
}

save_plot_both <- function(plot_object, stem, width, height) {
  ggsave(file.path(PLOT_DIR, paste0(stem, ".pdf")), plot_object,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(PLOT_DIR, paste0(stem, ".tiff")), plot_object,
         width = width, height = height, dpi = 600,
         compression = "lzw", bg = "white")
}

# -------------------- Load frozen model & training preprocessor --------------------
assert_file(TRAIN_FILE, "Project14 train.csv")
train_obj <- read_sample_by_gene(TRAIN_FILE, "GSE51588_subchondral_bone")
frozen <- load_frozen_model(MODEL_DIR)
fit <- frozen$fit
selected_method <- frozen$method
load_required_model_namespace(fit)

if (is.null(fit$subFeature)) {
  stop("Fitted model lacks subFeature field; cannot safely verify predictors.")
}
if (!all(as.character(fit$subFeature) %in% FINAL_GENES)) {
  stop("Frozen model uses variables outside six‑gene signature: ",
       paste(setdiff(as.character(fit$subFeature), FINAL_GENES), collapse = ", "))
}

pp_file <- file.path(MODEL_DIR, "model.preprocessing.rds")
if (file.exists(pp_file)) {
  pp <- readRDS(pp_file)
  pp_source <- pp_file
} else {
  warning(
    "Project14 model.preprocessing.rds missing; reconstruct median/mean/sd from train.csv.\n",
    "Prediction‑reproduction audit will block SHAP if mismatch occurs."
  )
  pp <- fit_preprocessor(train_obj$x)
  pp_source <- "Reconstructed_from_train.csv"
}

train_scaled <- apply_preprocessor(train_obj$x, pp)
train_prob <- predict_oa_probability(fit, train_scaled)
prediction_audit(train_prob, rownames(train_scaled), selected_method, MODEL_DIR)

# SHAP background is strictly sampled from training‑set scaled data
if (nrow(train_scaled) > BACKGROUND_MAX_N) {
  set.seed(20260929L)
  bg_index <- sample(seq_len(nrow(train_scaled)), BACKGROUND_MAX_N)
  background <- train_scaled[bg_index, , drop = FALSE]
} else {
  background <- train_scaled
}

# -------------------- Assemble cohorts for interpretation --------------------
cohorts <- list(
  GSE51588_subchondral_bone = list(
    x = train_scaled,
    meta = train_obj$meta,
    role = "Training‑set explanation"
  )
)

if (file.exists(TEST_FILE)) {
  test_obj <- read_sample_by_gene(TEST_FILE, "External_test")
  test_scaled <- apply_preprocessor(test_obj$x, pp)
  for (cohort_name in unique(test_obj$meta$Cohort)) {
    idx <- which(test_obj$meta$Cohort == cohort_name)
    cohorts[[cohort_name]] <- list(
      x = test_scaled[idx, , drop = FALSE],
      meta = test_obj$meta[idx, , drop = FALSE],
      role = "Independent external validation"
    )
  }
}

if (file.exists(PROJECT21_FILE)) {
  p21 <- read_gene_by_sample(PROJECT21_FILE, "GSE55457_synovium_OA_subset")
  cohorts[["GSE55457_synovium_OA_subset"]] <- list(
    x = apply_preprocessor(p21$x, pp),
    meta = p21$meta,
    role = "Independent synovium validation"
  )
}

if (RUN_PROJECT22_SENSITIVITY && file.exists(PROJECT22_FILE)) {
  p22 <- read_gene_by_sample(PROJECT22_FILE, "GSE12021_synovium_reanalysis")
  cohorts[["GSE12021_synovium_reanalysis"]] <- list(
    x = apply_preprocessor(p22$x, pp),
    meta = p22$meta,
    role = "Sensitivity‑analysis; overlaps with GSE55457"
  )
}
# Note: Project20 excluded: contains mixed RA samples mis‑annotated as OA.

# -------------------- Compute exact SHAP per cohort --------------------
shap_objects    <- list()
shapviz_objects <- list()
importance_list <- list()
direction_list  <- list()
prediction_list <- list()

for (cohort_name in names(cohorts)) {
  cat("Calculating SHAP for cohort: ", cohort_name, "\n")
  obj  <- cohorts[[cohort_name]]
  x    <- obj$x
  meta <- obj$meta
  prob <- predict_oa_probability(fit, x)

  ks <- kernelshap::kernelshap(
    object = fit,
    X = as.data.frame(x, check.names = FALSE),
    bg_X = as.data.frame(background, check.names = FALSE),
    pred_fun = predict_oa_probability,
    exact = TRUE,
    verbose = FALSE
  )
  sv <- shapviz::shapviz(ks)

  shap_objects[[cohort_name]]    <- ks
  shapviz_objects[[cohort_name]] <- sv
  saveRDS(ks, file.path(RDS_DIR, paste0(cohort_name, "_kernelshap.rds")))

  shap_matrix <- as.matrix(ks$S)
  shap_matrix <- shap_matrix[, FINAL_GENES, drop = FALSE]
  rownames(shap_matrix) <- rownames(x)

  shap_wide <- data.frame(
    Sample = rownames(shap_matrix),
    Cohort = cohort_name,
    Group  = as.character(meta[rownames(shap_matrix), "Group"]),
    Predicted_probability = prob,
    shap_matrix,
    check.names = FALSE
  )
  write.csv(
    shap_wide,
    file.path(TABLE_DIR, paste0(cohort_name, "_SHAP_values.csv")),
    row.names = FALSE
  )

  imp <- colMeans(abs(shap_matrix))
  importance_list[[cohort_name]] <- data.frame(
    Cohort = cohort_name,
    Role   = obj$role,
    Gene   = names(imp),
    Mean_abs_SHAP = as.numeric(imp),
    stringsAsFactors = FALSE
  )

  direction_list[[cohort_name]] <- do.call(rbind, lapply(FINAL_GENES, function(g) {
    data.frame(
      Cohort = cohort_name,
      Gene   = g,
      Mean_SHAP_Control = mean(shap_matrix[meta$Type == 0L, g], na.rm = TRUE),
      Mean_SHAP_OA      = mean(shap_matrix[meta$Type == 1L, g], na.rm = TRUE),
      Difference_OA_minus_Control =
        mean(shap_matrix[meta$Type == 1L, g], na.rm = TRUE) -
        mean(shap_matrix[meta$Type == 0L, g], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  prediction_list[[cohort_name]] <- data.frame(
    Sample = rownames(x), Cohort = cohort_name,
    Group = as.character(meta$Group), Probability = prob,
    stringsAsFactors = FALSE
  )

  # Global SHAP bar + bee‑swarm
  p_bar <- shapviz::sv_importance(
    sv, kind = "bar", max_display = length(FINAL_GENES), show_numbers = TRUE
  ) +
    scale_fill_gradient(low = "#8ECAE6", high = "#023047") +
    labs(
      title = paste0(cohort_name, ": global SHAP importance"),
      subtitle = "Mean absolute contribution to frozen‑model OA probability",
      x = "Mean |SHAP value|", y = NULL
    ) + theme_shap(12)

  p_bee <- shapviz::sv_importance(
    sv, kind = "bee", max_display = length(FINAL_GENES)
  ) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B") +
    labs(
      title = paste0(cohort_name, ": SHAP summary"),
      subtitle = "Red: higher standardized expression; blue: lower expression",
      x = "SHAP value for OA probability", y = NULL
    ) + theme_shap(12)

  combined <- p_bar + p_bee + patchwork::plot_annotation(tag_levels = "A")
  save_plot_both(combined, paste0(cohort_name, "_01_global_SHAP"), 13, 6.5)

  # Dependence plots for each gene
  dep_plots <- lapply(FINAL_GENES, function(g) {
    shapviz::sv_dependence(
      sv, v = g, color_var = "auto", alpha = 0.8, size = 2.2
    ) +
      labs(title = g, x = "Training‑standardized expression",
           y = "SHAP value for OA probability") +
      theme_shap(10)
  })
  dep_panel <- patchwork::wrap_plots(dep_plots, ncol = 3) +
    patchwork::plot_annotation(
      title = paste0(cohort_name, ": six‑gene SHAP dependence"),
      subtitle = "Association with model contribution; not a causal effect"
    )
  save_plot_both(dep_panel, paste0(cohort_name, "_02_dependence"), 13, 8)

  # Representative waterfall plots
  oa_index  <- which(meta$Type == 1L)
  con_index <- which(meta$Type == 0L)
  high_oa <- if (length(oa_index)) oa_index[which.max(prob[oa_index])] else which.max(prob)
  low_con <- if (length(con_index)) con_index[which.min(prob[con_index])] else which.min(prob)

  p_water_oa <- shapviz::sv_waterfall(
    sv, row_id = high_oa, max_display = length(FINAL_GENES)
  ) + labs(title = paste0("Representative high‑risk OA: ", rownames(x)[high_oa])) +
    theme_shap(11)

  p_water_con <- shapviz::sv_waterfall(
    sv, row_id = low_con, max_display = length(FINAL_GENES)
  ) + labs(title = paste0("Representative low‑risk control: ", rownames(x)[low_con])) +
    theme_shap(11)

  water_panel <- p_water_oa / p_water_con +
    patchwork::plot_annotation(
      title = paste0(cohort_name, ": individual prediction explanations"),
      subtitle = "Contributions relative to shared training‑set background"
    )
  save_plot_both(water_panel, paste0(cohort_name, "_03_waterfall"), 11, 10)
}

# -------------------- Cross‑cohort comparison --------------------
importance_df <- do.call(rbind, importance_list)
direction_df  <- do.call(rbind, direction_list)
prediction_df <- do.call(rbind, prediction_list)

write.csv(importance_df, file.path(TABLE_DIR, "All_cohorts_SHAP_importance.csv"), row.names = FALSE)
write.csv(direction_df,  file.path(TABLE_DIR, "All_cohorts_SHAP_direction.csv"),  row.names = FALSE)
write.csv(prediction_df, file.path(TABLE_DIR, "All_cohorts_predictions.csv"),    row.names = FALSE)

cohort_order <- names(cohorts)
importance_df$Cohort <- factor(importance_df$Cohort, levels = cohort_order)
importance_df$Gene   <- factor(importance_df$Gene,   levels = rev(FINAL_GENES))

p_heat <- ggplot(importance_df, aes(Cohort, Gene, fill = Mean_abs_SHAP)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Mean_abs_SHAP)), size = 3.2) +
  scale_fill_gradientn(
    colors = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"),
    name = "Mean |SHAP|"
  ) +
  labs(
    title = "Cross‑cohort stability of six‑gene model contributions",
    subtitle = "Project22 is repeated‑cohort sensitivity analysis, not independent validation",
    x = NULL, y = NULL
  ) +
  theme_shap(11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot_both(p_heat, "04_cross_cohort_SHAP_heatmap", 11.5, 6.8)

direction_long <- rbind(
  data.frame(
    Cohort = direction_df$Cohort, Gene = direction_df$Gene,
    Group = "Control", Mean_SHAP = direction_df$Mean_SHAP_Control
  ),
  data.frame(
    Cohort = direction_df$Cohort, Gene = direction_df$Gene,
    Group = "OA", Mean_SHAP = direction_df$Mean_SHAP_OA
  )
)
direction_long$Gene <- factor(direction_long$Gene, levels = FINAL_GENES)

p_direction <- ggplot(
  direction_long,
  aes(Gene, Mean_SHAP, color = Group, group = Group)
) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4) +
  facet_wrap(~Cohort, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(Control = "#2166AC", OA = "#B2182B")) +
  labs(
    title = "Mean signed SHAP contribution by phenotype and cohort",
    subtitle = "Positive values push frozen model toward OA prediction",
    x = NULL, y = "Mean SHAP value", color = NULL
  ) +
  theme_shap(11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top")
save_plot_both(p_direction, "05_cross_cohort_signed_SHAP", 12, 9)

# Main manuscript figure: training beeswarm + cross‑cohort heatmap
training_sv <- shapviz_objects[["GSE51588_subchondral_bone"]]
p_main_bee <- shapviz::sv_importance(
  training_sv, kind = "bee", max_display = length(FINAL_GENES)
) +
  scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B") +
  labs(
    title = "Training‑set SHAP summary",
    x = "SHAP value for OA probability", y = NULL
  ) + theme_shap(12)

main_figure <- p_main_bee / p_heat +
  patchwork::plot_annotation(
    title = paste0("Frozen six‑gene diagnostic model interpretation: ", selected_method),
    subtitle = "Exact model‑agnostic SHAP values using a common training background",
    tag_levels = "A"
  )
save_plot_both(main_figure, "Figure_SHAP_main", 13, 12)

# -------------------- Reproducibility records --------------------
cohort_summary <- do.call(rbind, lapply(names(cohorts), function(nm) {
  m <- cohorts[[nm]]$meta
  data.frame(
    Cohort = nm,
    Role = cohorts[[nm]]$role,
    N_expression_profiles = nrow(m),
    N_Control = sum(m$Type == 0L),
    N_OA = sum(m$Type == 1L),
    stringsAsFactors = FALSE
  )
}))
write.csv(cohort_summary, file.path(TABLE_DIR, "Cohort_summary.csv"), row.names = FALSE)

model_info <- data.frame(
  Item = c(
    "Selected_method", "Model_class", "Model_source",
    "Preprocessing_source", "Model_features", "SHAP_background",
    "Interpretation_boundary"
  ),
  Value = c(
    selected_method,
    paste(class(fit), collapse = ";"),
    frozen$source,
    pp_source,
    paste(as.character(fit$subFeature), collapse = ";"),
    paste0(nrow(background), " training expression profiles"),
    "Predictive contribution only; not causal or mechanistic effect"
  )
)
write.csv(model_info, file.path(TABLE_DIR, "SHAP_model_information.csv"), row.names = FALSE)

# English readme
writeLines(
  c(
    "Project23 SHAP interpretation finished.",
    paste0("Frozen final algorithm: ", selected_method),
    paste0("Model predictors: ", paste(fit$subFeature, collapse = ", ")),
    "SHAP background: sampled exclusively from Project14 training‑set scaled data.",
    "Project20: excluded due to mixed RA/OA annotation.",
    "Project22: patient overlap with Project21; used only for sensitivity analysis.",
    "SHAP values describe model prediction behaviour; DO NOT interpret as gene‑level causal effects for OA."
  ),
  file.path(RESULT_DIR, "README_results_interpretation.txt")
)

capture.output(sessionInfo(), file = file.path(RESULT_DIR, "sessionInfo.txt"))

cat("\n============================================================\n")
cat("Project23 SHAP analysis completed.\n")
cat("Final algorithm: ", selected_method, "\n")
cat("Result root directory: ", RESULT_DIR, "\n")
cat("Manuscript main figure: ", file.path(PLOT_DIR, "Figure_SHAP_main.pdf"), "\n")
cat("============================================================\n")
