# ---- Parameters ----
work_dir <- "D:\\"
gene_file_path <- "IntersectionGenes.txt"
output_train <- "train.csv"
output_test  <- "test.csv"

setwd(work_dir)

library(limma)
set.seed(20260929)

# ---- Step 1: Declare train / validation files ----
train_file        <- "Sample Type Matrix.csv"
validation_files  <- c("GSE114007.csv", "GSE63359.csv")

required_files <- c(train_file, validation_files, gene_file_path)
missing_files  <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

cat("Step 1/6: train / validation files identified.\n")
cat("Training file:        ", train_file, "\n")
cat("Validation files:    ", paste(validation_files, collapse = ", "), "\n")

# ---- Step 2: I/O helpers ----
read_expression_matrix <- function(in_file) {
  df <- read.csv(in_file, header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE)
  if (nrow(df) < 2 || ncol(df) < 2) {
    stop(in_file, ": expression matrix dimension is abnormal.")
  }
  gene_id <- trimws(as.character(df[[1]]))
  if (anyNA(gene_id) || any(!nzchar(gene_id))) {
    warning(in_file, ": empty gene names in first column, those rows are dropped.")
  }
  mat_df <- df[, -1, drop = FALSE]
  mat <- suppressWarnings(
    as.matrix(data.frame(lapply(mat_df, as.numeric), check.names = FALSE))
  )
  rownames(mat) <- gene_id
  colnames(mat) <- colnames(mat_df)
  keep <- !is.na(rownames(mat)) & nzchar(rownames(mat)) &
          rowSums(!is.na(mat)) > 0
  mat <- mat[keep, , drop = FALSE]
  if (anyDuplicated(rownames(mat))) {
    mat <- limma::avereps(mat)
  }
  if (any(!is.finite(mat), na.rm = TRUE)) {
    stop(in_file, ": non-finite values (Inf/-Inf) detected in expression values.")
  }
  mat
}

add_cohort_prefix <- function(mat, in_file) {
  cohort_tag <- tools::file_path_sans_ext(basename(in_file))
  colnames(mat) <- paste0(cohort_tag, "_", colnames(mat))
  mat
}

# ---- Step 3: Read train and validation separately ----
cat("Step 2/6: Reading training and external validation matrices.\n")

train_expr <- read_expression_matrix(train_file)
train_expr <- add_cohort_prefix(train_expr, train_file)

validation_expr_list <- lapply(validation_files, function(in_file) {
  x <- read_expression_matrix(in_file)
  add_cohort_prefix(x, in_file)
})
names(validation_expr_list) <-
  tools::file_path_sans_ext(basename(validation_files))

# No merging here; no shared preprocessing parameters are estimated.

# ---- Step 4: Align genes across cohorts ----
cat("Step 3/6: Aligning genes shared across cohorts.\n")

common_genes <- Reduce(
  intersect,
  c(list(rownames(train_expr)),
    lapply(validation_expr_list, rownames))
)
if (length(common_genes) < 5) {
  stop("Fewer than 5 genes shared across cohorts; check platform annotation.")
}

target_genes <- read.table(
  gene_file_path, header = FALSE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)[, 1]
target_genes <- unique(trimws(as.character(target_genes)))
target_genes <- target_genes[!is.na(target_genes) & nzchar(target_genes)]

select_genes <- target_genes[target_genes %in% common_genes]
if (length(select_genes) < 2) {
  stop("Fewer than 2 overlapping genes between IntersectionGenes.txt and matrices.")
}

cat(sprintf("Shared genes: %d; final candidate genes: %d.\n",
            length(common_genes), length(select_genes)))

# ---- Step 5: Build train / validation matrices ----
cat("Step 4/6: Building train and external validation matrices.\n")
cat("No joint ComBat between train and validation.\n")

train_mat <- t(train_expr[select_genes, , drop = FALSE])

validation_mats <- lapply(validation_expr_list, function(x) {
  t(x[select_genes, , drop = FALSE])
})

# Row-bind validation cohorts only after independent read + gene alignment.
# No shared mean / SD / ComBat parameters are computed.
test_mat <- do.call(rbind, validation_mats)

if (!identical(colnames(train_mat), colnames(test_mat))) {
  stop("Gene order differs between train and validation matrices.")
}
if (anyDuplicated(rownames(train_mat))) stop("Duplicated sample names in training set.")
if (anyDuplicated(rownames(test_mat))) stop("Duplicated sample names in validation set.")
if (length(intersect(rownames(train_mat), rownames(test_mat))) > 0) {
  stop("Overlapping sample names between training and validation sets.")
}

# ---- Step 6: Parse Con/OA labels from sample names ----
cat("Step 5/6: Parsing Con/OA labels from sample names.\n")

parse_type <- function(ids) {
  suffix <- tolower(sub("^.*[_.-]", "", ids))
  control_labels <- c("con", "control", "normal", "healthy")
  case_labels    <- c("oa", "case", "disease", "tra")
  out <- rep(NA_integer_, length(ids))
  out[suffix %in% control_labels] <- 0L
  out[suffix %in% case_labels]    <- 1L
  if (anyNA(out)) {
    stop(
      "Cannot infer Con/OA label for: ",
      paste(ids[is.na(out)], collapse = ", "),
      ". Unknown samples are NOT auto-coded as OA."
    )
  }
  out
}

train_type <- parse_type(rownames(train_mat))
test_type  <- parse_type(rownames(test_mat))

if (!identical(sort(unique(train_type)), c(0L, 1L))) {
  stop("Training set must contain both Con and OA samples.")
}

train_out <- cbind(train_mat, Type = train_type)
test_out  <- cbind(test_mat,  Type = test_type)

rownames(train_out) <- gsub("merge_", "Train.", rownames(train_out))
rownames(test_out)  <- rownames(test_out)

# Split validation labels by cohort (reuse already-parsed labels)
val_split_idx  <- cumsum(vapply(validation_mats, nrow, integer(1)))
val_type_list  <- split(
  test_type,
  rep(seq_along(validation_mats),
      times = diff(c(0, val_split_idx)))
)

# ---- Step 7: Backup existing outputs and write new ones ----
cat("Step 6/6: Writing train.csv and test.csv.\n")

backup_existing_file <- function(path) {
  if (file.exists(path)) {
    backup_path <- paste0(
      tools::file_path_sans_ext(path),
      "_backup_before_fix_20260929.",
      tools::file_ext(path)
    )
    if (!file.exists(backup_path)) {
      file.copy(path, backup_path, overwrite = FALSE)
      cat("Backed up: ", backup_path, "\n")
    }
  }
}
backup_existing_file(output_train)
backup_existing_file(output_test)

# Keep Type numeric; first column = sample_id
write_sample_csv <- function(mat, outfile) {
  df <- as.data.frame(mat, check.names = FALSE)
  df <- cbind(sample_id = rownames(df), df)
  write.table(df, file = outfile, sep = ",", quote = FALSE, row.names = FALSE)
}
write_sample_csv(train_out, output_train)
write_sample_csv(test_out,  output_test)

# ---- Cohort manifest ----
cohort_manifest <- data.frame(
  File = c(train_file, validation_files),
  Cohort = c(
    tools::file_path_sans_ext(basename(train_file)),
    tools::file_path_sans_ext(basename(validation_files))
  ),
  Role = c(
    "Training",
    rep("External validation", length(validation_files))
  ),
  Samples = c(
    nrow(train_mat),
    vapply(validation_mats, nrow, integer(1))
  ),
  Controls = c(
    sum(train_type == 0),
    vapply(val_type_list, function(x) sum(x == 0L), integer(1))
  ),
  OA = c(
    sum(train_type == 1),
    vapply(val_type_list, function(x) sum(x == 1L), integer(1))
  ),
  stringsAsFactors = FALSE
)
write.table(cohort_manifest, file = "cohort_manifest.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
cat(sprintf("Training set: %d samples (Con=%d, OA=%d).\n",
            nrow(train_mat), sum(train_type == 0), sum(train_type == 1)))
cat(sprintf("Validation set: %d samples (Con=%d, OA=%d).\n",
            nrow(test_mat), sum(test_type == 0), sum(test_type == 1)))
cat("Training and validation are processed independently; no joint ComBat.\n")
