suppressPackageStartupMessages({
  library(randomForestSRC)
  library(glmnet)
  library(plsRglm)
  library(gbm)
  library(caret)
  library(mboost)
  library(e1071)
  library(MASS)
  library(xgboost)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(pROC)
  library(circlize)
})

input_files <- c("train.csv", "test.csv", "refer.txt")
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing input files in working directory: ", paste(missing_files, collapse = ", "))
}

RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  method = gsub(" ", "", method)
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)

  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )

  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")

  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "mode" = mode,
              "classVar" = classVar)
  args = c(args, method_param)

  obj <- do.call(what = paste0("Run", method_name),
                 args = args)

  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{message("\n")}
  return(obj)
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  inner_nfolds <- max(3L, min(10L, as.integer(min(table(Train_label[[classVar]])))))
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha,
                     nfolds = inner_nfolds, type.measure = "auc")
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  fit <- step(glm(formula = Train_label[[classVar]] ~ .,
                  family = "binomial",
                  data = as.data.frame(Train_set)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLDA <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  lda_folds <- max(3L, min(5L, as.integer(min(table(data[[classVar]])))))
  fit = train(eval(parse(text = paste(classVar, "~."))),
              data = data,
              method="lda",
              trControl = trainControl(method = "cv", number = lda_folds))
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial())
  boost_folds <- max(3L, min(5L, as.integer(min(table(data[[classVar]])))))
  cvm <- cvrisk(fit, papply = lapply,
                folds = cv(model.weights(fit), type = "kfold", B = boost_folds))
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = max(mstop(cvm), 40)))
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  fit <- plsRglm(Train_label[[classVar]],
                 as.data.frame(Train_set),
                 modele = "pls-glm-logistic",
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  rf_nodesize = 5
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
               data = cbind(Train_set, Train_label[classVar]),
               ntree = 1000, nodesize = rf_nodesize,
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunGBM <- function(Train_set, Train_label, mode, classVar){
  gbm_folds <- max(3L, min(5L, as.integer(min(table(Train_label[[classVar]])))))
  gbm_cores <- max(1L, min(2L, parallel::detectCores(logical = FALSE)))
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = gbm_folds, n.cores = gbm_cores)
  best <- max(1L, which.min(fit$cv.error))
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = gbm_cores)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  xgb_folds <- max(3L, min(5L, as.integer(min(table(Train_label[[classVar]])))))
  indexes = createFolds(factor(Train_label[[classVar]]), k = xgb_folds, list=T)
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgb.DMatrix(data = Train_set[-pt, ],
                         label = Train_label[[classVar]][-pt])
    dtest = xgb.DMatrix(data = Train_set[pt, ],
                        label = Train_label[[classVar]][pt])
    watchlist <- list(train=dtrain, test=dtest)

    bst <- xgb.train(data=dtrain,
                     max.depth=2, eta=1, nthread = 2, nrounds=10,
                     watchlist=watchlist,
                     objective = "binary:logistic", verbose = F)
    which.min(bst$evaluation_log$test_logloss)
  }))

  nround <- as.numeric(names(which.max(table(CV))))
  fit <- xgboost(data = Train_set,
                 label = Train_label[[classVar]],
                 max.depth = 2, eta = 1, nthread = 2, nrounds = nround,
                 objective = "binary:logistic", verbose = F)
  fit$subFeature = colnames(Train_set)

  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))),
                    data = data)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

standarize.fun <- function(indata, centerFlag, scaleFlag) {
  scale(indata, center=centerFlag, scale=scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }

  if (is.null(centerFlags)){
    centerFlags = F; message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags)==1){
    centerFlags = rep(centerFlags, length(data)); message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data); message("match centerFlags with cohort by order\n")
  }

  if (is.null(scaleFlags)){
    scaleFlags = F; message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags)==1){
    scaleFlags = rep(scaleFlags, length(data)); message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data); message("match scaleFlags with cohort by order\n")
  }

  centerFlags <- centerFlags[names(data)]; scaleFlags <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data, centerFlag = centerFlags, scaleFlag = scaleFlags, SIMPLIFY = F)
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}

ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)],
    "glm" = names(coef(fit)),
    "svm.formula" = fit$subFeature,
    "train" = fit$coefnames,
    "glmboost" = names(coef(fit)[abs(coef(fit))>0]),
    "plsRglmmodel" = {
      cc <- as.matrix(fit$Coeffs)
      rownames(cc)[rowSums(abs(cc), na.rm = TRUE) > 0]
    },
    "rfsrc" = {
      imp <- fit$importance
      if (is.matrix(imp) || is.data.frame(imp)) imp <- rowMeans(as.matrix(imp), na.rm = TRUE)
      names(imp)[is.finite(imp) & imp > 0.01]
    },
    "gbm" = {
      ss <- summary.gbm(fit, plotit = FALSE)
      rownames(ss)[ss$rel.inf > 0]
    },
    "xgb.Booster" = fit$subFeature,
    "naiveBayes" = fit$subFeature
  ))

  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}

ReportFeatures <- function(fit) {
  if (!is.null(fit$selectedFeature)) return(fit$selectedFeature)
  fit$subFeature
}

PredictSVMProbability <- function(fit, new_data) {
  pred <- predict(fit, as.data.frame(new_data), probability = TRUE)
  prob <- attr(pred, "probabilities")
  if (is.null(prob)) stop("SVM did not return probability, check probability=TRUE.")
  positive_col <- if ("1" %in% colnames(prob)) "1" else colnames(prob)[ncol(prob)]
  as.numeric(prob[, positive_col])
}

CalPredictScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, fit$subFeature]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = predict(fit, type = 'response', as.matrix(new_data)),
    "glm"         = predict(fit, type = 'response', as.data.frame(new_data)),
    "svm.formula" = PredictSVMProbability(fit, new_data),
    "train"       = predict(fit, new_data, type = "prob")[[2]],
    "glmboost"    = predict(fit, type = "response", as.data.frame(new_data)),
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
    "gbm"          = predict(fit, type = 'response', as.data.frame(new_data)),
    "xgb.Booster" = predict(fit, as.matrix(new_data)),
    "naiveBayes" = predict(object = fit, type = "raw", newdata = new_data)[, "1"]
  ))
  RS = as.numeric(as.vector(RS))
  names(RS) = rownames(new_data)
  return(RS)
}

PredictClass <- function(fit, new_data){
  new_data <- new_data[, fit$subFeature]
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = predict(fit, type = 'class', as.matrix(new_data)),
    "glm"         = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                           yes = "1", no = "0"),
    "svm.formula" = predict(fit, as.data.frame(new_data), decision.values = T),
    "train"       = predict(fit, new_data, type = "raw"),
    "glmboost"    = predict(fit, type = "class", as.data.frame(new_data)),
    "plsRglmmodel" = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
    "gbm"          = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"),
    "xgb.Booster" = ifelse(test = predict(fit, as.matrix(new_data))>0.5,
                           yes = "1", no = "0"),
    "naiveBayes" = predict(object = fit, type = "class", newdata = new_data)
  ))
  label = as.character(as.vector(label))
  names(label) = rownames(new_data)
  return(label)
}

RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){

  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }

  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(Train_set[, fit$subFeature],
                                 Test_set[, fit$subFeature])

    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    Test_label <- rbind.data.frame(Train_label[,c(cohortVar, classVar)],
                                   Test_label[,c(cohortVar, classVar)])
    Test_label[,1] <- factor(Test_label[,1],
                             levels = c(unique(Train_label[,cohortVar]), setdiff(unique(Test_label[,cohortVar]),unique(Train_label[,cohortVar]))))
  } else {
    new_data <- Test_set[, fit$subFeature]
  }

  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))

  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  unlist(lapply(Predict.out, function(data){
    as.numeric(auc(suppressMessages(roc(data[[classVar]], data$RS))))
  }))
}

SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)

  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = T, numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)

  Heatmap(as.matrix(Cindex_mat), name = "AUC",
          right_annotation = row_ha,
          top_annotation = col_ha,
          col = c("#4195C1", "#FFFFFF", "#FFBC90"),
          rect_gp = gpar(col = "black", lwd = 1),
          cluster_columns = cluster_columns, cluster_rows = cluster_rows,
          show_column_names = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) {
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}

Train_data <- read.table("train.csv", header = TRUE, sep = ",", check.names = FALSE,
                         row.names = 1, stringsAsFactors = FALSE)
Test_data <- read.table("test.csv", header = TRUE, sep = ",", check.names = FALSE,
                        row.names = 1, stringsAsFactors = FALSE)
Train_expr <- as.matrix(Train_data[, 1:(ncol(Train_data)-1), drop = FALSE])
Test_expr <- as.matrix(Test_data[, 1:(ncol(Test_data)-1), drop = FALSE])
storage.mode(Train_expr) <- "double"
storage.mode(Test_expr) <- "double"
Train_class <- data.frame(Type = as.integer(as.character(Train_data[, ncol(Train_data)])),
                          row.names = rownames(Train_data))
Test_class <- data.frame(Type = as.integer(as.character(Test_data[, ncol(Test_data)])),
                         row.names = rownames(Test_data))
Test_class$Cohort <- sub("_.*$", "", rownames(Test_class))
Test_class <- Test_class[, c("Cohort", "Type"), drop = FALSE]
if (!all(Train_class$Type %in% c(0L, 1L)) || !all(Test_class$Type %in% c(0L, 1L)))
  stop("Type must be strictly encoded as 0/1.")
if (length(unique(Train_class$Type)) != 2) stop("Training set must contain both 0 and 1.")
if (anyDuplicated(rownames(Train_expr)) || anyDuplicated(rownames(Test_expr)))
  stop("Duplicated sample names exist.")
if (length(intersect(rownames(Train_expr), rownames(Test_expr))))
  stop("Duplicated samples between training and test set.")

model_genes <- colnames(Train_expr)
missing_in_test <- setdiff(model_genes, colnames(Test_expr))
if (length(missing_in_test))
  stop("Test set missing features from training set: ", paste(missing_in_test, collapse = ", "))
Test_expr <- Test_expr[, model_genes, drop = FALSE]

fit_preprocessor <- function(x) {
  med <- apply(x, 2, median, na.rm = TRUE)
  xi <- x
  for (j in seq_len(ncol(xi))) xi[!is.finite(xi[, j]), j] <- med[j]
  mu <- colMeans(xi)
  sigma <- apply(xi, 2, sd)
  sigma[!is.finite(sigma) | sigma == 0] <- 1
  list(median = med, mean = mu, sd = sigma, genes = colnames(x))
}

apply_preprocessor <- function(x, pp) {
  x <- x[, pp$genes, drop = FALSE]
  for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- pp$median[j]
  x <- sweep(x, 2, pp$mean, "-")
  sweep(x, 2, pp$sd, "/")
}

fit_one_method <- function(method, x, label, classVar = "Type") {
  parts <- strsplit(method, "\\+")[[1]]
  if (length(parts) == 1) parts <- c("simple", parts)
  selector <- parts[1]
  learner <- parts[2]
  vars <- if (selector == "simple") colnames(x) else
    RunML(method = selector, Train_set = x, Train_label = label,
          mode = "Variable", classVar = classVar)
  vars <- intersect(unique(vars), colnames(x))
  if (length(vars) <= min.selected.var) return(NULL)
  fit <- RunML(method = learner, Train_set = x[, vars, drop = FALSE],
               Train_label = label, mode = "Model", classVar = classVar)
  selected_vars <- tryCatch(
    intersect(unique(ExtractVar(fit)), vars),
    error = function(e) character(0)
  )
  if (length(selected_vars) <= min.selected.var) return(NULL)
  fit$inputFeature <- vars
  fit$selectedFeature <- selected_vars
  fit$subFeature <- vars
  fit
}

fixed_auc <- function(y, score) {
  ok <- !is.na(y) & is.finite(score)
  if (sum(ok) < 5 || length(unique(y[ok])) < 2)
    return(c(AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_))
  rr <- pROC::roc(y[ok], score[ok], levels = c(0, 1),
                  direction = "<", quiet = TRUE)
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(rr)))
  c(AUC = as.numeric(pROC::auc(rr)), CI_low = ci[1], CI_high = ci[3])
}

training_threshold <- function(y, score) {
  ok <- !is.na(y) & is.finite(score)
  if (sum(ok) < 5L || length(unique(y[ok])) < 2L) return(NA_real_)
  rr <- pROC::roc(y[ok], score[ok], levels = c(0, 1),
                  direction = "<", quiet = TRUE)
  coord <- pROC::coords(
    rr, "best", best.method = "youden", ret = "threshold",
    transpose = FALSE
  )
  if (is.data.frame(coord) || is.matrix(coord)) {
    if ("threshold" %in% colnames(coord)) return(as.numeric(coord[1, "threshold"]))
    return(as.numeric(coord[1, 1]))
  }
  if (is.list(coord)) {
    if (!is.null(coord$threshold)) return(as.numeric(unlist(coord$threshold))[1])
    return(as.numeric(unlist(coord, use.names = FALSE))[1])
  }
  as.numeric(coord)[1]
}

safe_rbind <- function(x) {
  if (!length(x)) return(NULL)
  do.call(rbind, x)
}

evaluate_method_cv <- function(method, raw_x, label, train_indices,
                               stage, seed_offset = 0L) {
  n <- nrow(raw_x)
  pred_sum <- rep(0, n)
  pred_count <- integer(n)
  selection_rows <- list()
  error_rows <- list()
  for (fold_name in names(train_indices)) {
    tr <- train_indices[[fold_name]]
    va <- setdiff(seq_len(n), tr)
    pp <- fit_preprocessor(raw_x[tr, , drop = FALSE])
    x_tr <- apply_preprocessor(raw_x[tr, , drop = FALSE], pp)
    x_va <- apply_preprocessor(raw_x[va, , drop = FALSE], pp)
    y_tr <- label[tr, , drop = FALSE]
    set.seed(GLOBAL_SEED + seed_offset + match(fold_name, names(train_indices)))
    fit <- tryCatch(fit_one_method(method, x_tr, y_tr, "Type"),
                    error = function(e) e)
    if (inherits(fit, "error") || is.null(fit)) {
      error_rows[[length(error_rows)+1L]] <- data.frame(
        Stage = stage, Fold = fold_name, Method = method,
        Error = if (inherits(fit, "error")) conditionMessage(fit)
        else "Too few selected variables",
        stringsAsFactors = FALSE)
      next
    }
    pr <- tryCatch(CalPredictScore(fit, x_va), error = function(e) e)
    if (inherits(pr, "error")) {
      error_rows[[length(error_rows)+1L]] <- data.frame(
        Stage = stage, Fold = fold_name, Method = method,
        Error = conditionMessage(pr), stringsAsFactors = FALSE)
      next
    }
    pred_sum[va] <- pred_sum[va] + as.numeric(pr)
    pred_count[va] <- pred_count[va] + 1L
    selection_rows[[length(selection_rows)+1L]] <- data.frame(
      Stage = stage, Fold = fold_name, Method = method,
      Gene = ReportFeatures(fit), stringsAsFactors = FALSE)
  }
  pred <- pred_sum / ifelse(pred_count == 0, NA, pred_count)
  met <- fixed_auc(label$Type, pred)
  list(Method = method, Predictions = pred, Counts = pred_count,
       AUC = unname(met["AUC"]), CI_low = unname(met["CI_low"]),
       CI_high = unname(met["CI_high"]),
       Selections = safe_rbind(selection_rows), Errors = safe_rbind(error_rows))
}

compare_methods_inner_cv <- function(raw_x, label, methods, folds, stage) {
  results <- vector("list", length(methods)); names(results) <- methods
  metrics <- list(); selections <- list(); errors <- list()
  for (i in seq_along(methods)) {
    method <- methods[i]
    cat(stage, ": ", i, "/", length(methods), " ", method, "\n", sep = "")
    flush.console()
    ev <- evaluate_method_cv(method, raw_x, label, folds, stage,
                             seed_offset = i * 1000L)
    results[[method]] <- ev
    metrics[[method]] <- data.frame(Method = method, AUC = ev$AUC,
                                    CI_low = ev$CI_low, CI_high = ev$CI_high,
                                    OOF_N = sum(is.finite(ev$Predictions)), stringsAsFactors = FALSE)
    if (!is.null(ev$Selections)) selections[[method]] <- ev$Selections
    if (!is.null(ev$Errors)) errors[[method]] <- ev$Errors
  }
  metric_df <- safe_rbind(metrics)
  metric_df <- metric_df[is.finite(metric_df$AUC), , drop = FALSE]
  metric_df <- metric_df[order(-metric_df$AUC, metric_df$Method), , drop = FALSE]
  list(Results = results, Metrics = metric_df,
       Selections = safe_rbind(selections), Errors = safe_rbind(errors))
}

methodRT <- read.table("refer.txt", header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
methods <- unique(gsub("-| ", "", methodRT$Model))
classVar <- "Type"
cat("Detected ", length(methods), " candidate algorithm combinations. Strict Nested‑CV involves heavy repeated fitting. Long pause on one algorithm is not equivalent to program hang.\n", sep="")
flush.console()

set.seed(GLOBAL_SEED)
outer_indices <- caret::createMultiFolds(
  factor(Train_class$Type), k = OUTER_FOLDS, times = OUTER_REPEATS)
n <- nrow(Train_expr)
nested_sum <- rep(0, n)
nested_count <- integer(n)
outer_selection <- list()
outer_error_log <- list()
outer_gene_log <- list()
outer_fold_prediction <- list()

for (outer_i in seq_along(outer_indices)) {
  outer_name <- names(outer_indices)[outer_i]
  cat("\n===== Outer fold: ", outer_name, " =====\n", sep = "")
  outer_tr <- outer_indices[[outer_i]]
  outer_va <- setdiff(seq_len(n), outer_tr)
  x_outer_raw <- Train_expr[outer_tr, , drop = FALSE]
  y_outer <- Train_class[outer_tr, , drop = FALSE]
  inner_k <- min(INNER_FOLDS, as.integer(min(table(y_outer$Type))))
  if (inner_k < 3L) stop("Minor‑class sample size in outer training fold insufficient for inner CV.")
  set.seed(GLOBAL_SEED + outer_i * 10000L)
  inner_indices <- caret::createMultiFolds(
    factor(y_outer$Type), k = inner_k, times = INNER_REPEATS)
  inner <- compare_methods_inner_cv(
    x_outer_raw, y_outer, methods, inner_indices,
    stage = paste0("Outer_", outer_name, "_InnerCV"))
  if (!nrow(inner$Metrics)) stop("All inner models failed for outer fold ", outer_name, ".")
  if (!is.null(inner$Errors)) outer_error_log[[outer_name]] <- inner$Errors
  if (!is.null(inner$Selections)) outer_gene_log[[outer_name]] <- inner$Selections

  ranked_methods <- inner$Metrics$Method
  selected_fit <- NULL; selected_method <- NA_character_; selected_prediction <- NULL
  pp_outer <- fit_preprocessor(x_outer_raw)
  x_outer_train <- apply_preprocessor(x_outer_raw, pp_outer)
  x_outer_valid <- apply_preprocessor(Train_expr[outer_va, , drop = FALSE], pp_outer)
  for (candidate in ranked_methods) {
    set.seed(GLOBAL_SEED + outer_i * 100000L + match(candidate, methods))
    candidate_fit <- tryCatch(
      fit_one_method(candidate, x_outer_train, y_outer, classVar),
      error = function(e) e)
    if (inherits(candidate_fit, "error") || is.null(candidate_fit)) next
    candidate_pred <- tryCatch(CalPredictScore(candidate_fit, x_outer_valid),
                               error = function(e) e)
    if (inherits(candidate_pred, "error")) next
    selected_fit <- candidate_fit
    selected_method <- candidate
    selected_prediction <- as.numeric(candidate_pred)
    break
  }
  if (is.null(selected_fit)) stop("Cannot fit any candidate model for outer fold ", outer_name, ".")
  nested_sum[outer_va] <- nested_sum[outer_va] + selected_prediction
  nested_count[outer_va] <- nested_count[outer_va] + 1L
  outer_selection[[outer_name]] <- data.frame(
    OuterFold = outer_name, SelectedMethod = selected_method,
    InnerAUC = inner$Metrics$AUC[match(selected_method, inner$Metrics$Method)],
    SelectedGenes = paste(ReportFeatures(selected_fit), collapse = ";"),
    stringsAsFactors = FALSE)
  outer_fold_prediction[[outer_name]] <- data.frame(
    OuterFold = outer_name, Sample = rownames(Train_expr)[outer_va],
    Truth = Train_class$Type[outer_va], Probability = selected_prediction,
    SelectedMethod = selected_method, stringsAsFactors = FALSE)

  saveRDS(list(
    LastCompletedOuterFold = outer_name,
    NestedSum = nested_sum, NestedCount = nested_count,
    OuterSelection = outer_selection, OuterErrors = outer_error_log,
    OuterGenes = outer_gene_log, OuterPredictions = outer_fold_prediction
  ), "model.nestedCV_checkpoint.rds")
}

nested_probability <- nested_sum / ifelse(nested_count == 0, NA, nested_count)
nested_metric <- fixed_auc(Train_class$Type, nested_probability)
nested_performance <- data.frame(
  Evaluation = "Strict nested CV outer OOF",
  AUC = nested_metric["AUC"], CI_low = nested_metric["CI_low"],
  CI_high = nested_metric["CI_high"], N = sum(is.finite(nested_probability)),
  OuterFolds = OUTER_FOLDS, OuterRepeats = OUTER_REPEATS,
  InnerFolds = INNER_FOLDS, InnerRepeats = INNER_REPEATS)
write.table(data.frame(id = rownames(Train_expr), Type = Train_class$Type,
                       Probability = nested_probability,
                       PredictionCount = nested_count),
            "model.nestedCV_OOF_predictions.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(nested_performance, "model.nestedCV_performance.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(safe_rbind(outer_selection), "model.nestedCV_selectedAlgorithms.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(safe_rbind(outer_fold_prediction), "model.nestedCV_foldPredictions.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

set.seed(GLOBAL_SEED + 900000L)
final_inner_indices <- caret::createMultiFolds(
  factor(Train_class$Type), k = FINAL_INNER_FOLDS,
  times = FINAL_INNER_REPEATS)
final_inner <- compare_methods_inner_cv(
  Train_expr, Train_class, methods, final_inner_indices,
  stage = "FullTraining_InnerCV")
if (!nrow(final_inner$Metrics)) stop("All algorithms failed in full‑training‑set inner CV.")
best_method <- final_inner$Metrics$Method[1]
cat("Final algorithm determined only by full‑training‑set inner CV: ", best_method, "\n", sep="")
write.table(final_inner$Metrics, "model.innerCV_method_comparison.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

full_pp <- fit_preprocessor(Train_expr)
Train_set <- apply_preprocessor(Train_expr, full_pp)
Test_set <- apply_preprocessor(Test_expr, full_pp)
model <- list(); final_fit_errors <- list()
for (method in final_inner$Metrics$Method) {
  set.seed(GLOBAL_SEED + match(method, methods))
  fit <- tryCatch(fit_one_method(method, Train_set, Train_class, classVar),
                  error = function(e) e)
  if (inherits(fit, "error") || is.null(fit)) {
    final_fit_errors[[method]] <- data.frame(
      Stage = "FinalFullTrainingFit", Fold = "FullTraining", Method = method,
      Error = if (inherits(fit, "error")) conditionMessage(fit)
      else "Too few selected variables", stringsAsFactors = FALSE)
  } else model[[method]] <- fit
}
if (!length(model)) stop("No model fitted successfully on full training set.")
if (!best_method %in% names(model)) {
  best_method <- final_inner$Metrics$Method[
    final_inner$Metrics$Method %in% names(model)][1]
  warning("Top inner‑CV model failed on full training set, fallback to: ", best_method)
}
selected_model <- model[[best_method]]
saveRDS(model, "model.MLmodel.rds")
saveRDS(selected_model, "model.selectedModel.rds")
saveRDS(full_pp, "model.preprocessing.rds")

logisticmodel <- lapply(model, function(fit) {
  vars <- ReportFeatures(fit)
  tryCatch({
    tmp <- glm(Train_class[[classVar]] ~ ., family = "binomial",
               data = as.data.frame(Train_set[, vars, drop = FALSE]))
    tmp$subFeature <- vars; tmp
  }, error = function(e) NULL)
})
logisticmodel <- logisticmodel[!vapply(logisticmodel, is.null, logical(1))]
saveRDS(logisticmodel, "model.logisticmodel.rds")

valid_methods <- intersect(final_inner$Metrics$Method, names(model))
all_ids <- c(rownames(Train_set), rownames(Test_set))
risk_mat <- matrix(NA_real_, nrow = length(all_ids), ncol = length(valid_methods),
                   dimnames = list(all_ids, valid_methods))
class_mat <- matrix(NA_character_, nrow = length(all_ids), ncol = length(valid_methods),
                    dimnames = list(all_ids, valid_methods))
thresholds <- setNames(rep(NA_real_, length(valid_methods)), valid_methods)
external_rows <- list()
inner_oof_matrix <- matrix(NA_real_, nrow = nrow(Train_expr),
                           ncol = length(valid_methods),
                           dimnames = list(rownames(Train_expr), valid_methods))

for (method in valid_methods) {
  inner_oof_matrix[, method] <- final_inner$Results[[method]]$Predictions
  thresholds[method] <- training_threshold(Train_class$Type,
                                            inner_oof_matrix[, method])
  fit <- model[[method]]
  train_score <- CalPredictScore(fit, Train_set)
  test_score <- CalPredictScore(fit, Test_set)
  risk_mat[rownames(Train_set), method] <- train_score
  risk_mat[rownames(Test_set), method] <- test_score
  class_mat[, method] <- ifelse(risk_mat[, method] >= thresholds[method], "1", "0")
  for (cohort in unique(Test_class$Cohort)) {
    ids <- rownames(Test_class)[Test_class$Cohort == cohort]
    met <- fixed_auc(Test_class[ids, "Type"], test_score[ids])
    pred <- ifelse(test_score[ids] >= thresholds[method], 1L, 0L)
    truth <- Test_class[ids, "Type"]
    external_rows[[length(external_rows)+1L]] <- data.frame(
      Method = method, Selected_Final_Model = method == best_method,
      Cohort = cohort, AUC = met["AUC"], CI_low = met["CI_low"],
      CI_high = met["CI_high"], Threshold = thresholds[method],
      Sensitivity = sum(pred == 1 & truth == 1) / sum(truth == 1),
      Specificity = sum(pred == 0 & truth == 0) / sum(truth == 0),
      N = length(ids), stringsAsFactors = FALSE)
  }
}

external_metrics <- safe_rbind(external_rows)
riskTab <- data.frame(id = rownames(risk_mat), risk_mat, check.names = FALSE)
classTab <- data.frame(id = rownames(class_mat), class_mat, check.names = FALSE)
write.table(riskTab, "model.riskMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(classTab, "model.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(id = rownames(inner_oof_matrix), Type = Train_class$Type,
                       inner_oof_matrix, check.names = FALSE),
            "model.OOFpredictions.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(external_metrics, "model.externalValidationMetrics.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(Method = names(thresholds), Threshold = thresholds),
            "model.training_thresholds.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

cohorts <- unique(Test_class$Cohort)
AUC_mat <- matrix(NA_real_, nrow = length(valid_methods),
                  ncol = 1 + length(cohorts),
                  dimnames = list(valid_methods, c("Train_InnerCV", cohorts)))
AUC_mat[, "Train_InnerCV"] <- final_inner$Metrics$AUC[
  match(valid_methods, final_inner$Metrics$Method)]
for (cohort in cohorts) {
  d <- external_metrics[external_metrics$Cohort == cohort, ]
  AUC_mat[d$Method, cohort] <- d$AUC
}
AUC_mat <- AUC_mat[order(-AUC_mat[, "Train_InnerCV"]), , drop = FALSE]
write.table(cbind(Method = rownames(AUC_mat), AUC_mat),
            "model.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

fea_list <- lapply(model[valid_methods], ReportFeatures)
fea_df <- safe_rbind(lapply(names(fea_list), function(method) {
  data.frame(features = fea_list[[method]], algorithm = method,
             Selected_Final_Model = method == best_method,
             stringsAsFactors = FALSE)
}))
write.table(fea_df, "model.genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

all_selection <- safe_rbind(c(outer_gene_log,
                              list(FullTraining = final_inner$Selections)))
if (!is.null(all_selection)) {
  selection_frequency <- aggregate(Fold ~ Stage + Method + Gene,
                                   all_selection,
                                   function(x) length(unique(x)))
  names(selection_frequency)[4] <- "Selected_Folds"
  write.table(selection_frequency, "model.gene_selection_frequency.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
}

performance_report <- final_inner$Metrics[
  match(rownames(AUC_mat), final_inner$Metrics$Method), ]
performance_report$External_Mean_AUC_Descriptive <-
  rowMeans(AUC_mat[, cohorts, drop = FALSE], na.rm = TRUE)
performance_report$Rank_by_FullTraining_InnerCV <- seq_len(nrow(performance_report))
performance_report$Selected_Final_Model <- performance_report$Method == best_method
performance_report$Nested_Pipeline_AUC <- unname(nested_metric["AUC"])
write.table(performance_report, "model.performance_detailed_report.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

outer_method_frequency <- as.data.frame(table(
  safe_rbind(outer_selection)$SelectedMethod), stringsAsFactors = FALSE)
colnames(outer_method_frequency) <- c("Method", "Selected_Outer_Folds")
outer_method_frequency$Selection_Frequency <-
  outer_method_frequency$Selected_Outer_Folds / length(outer_indices)
write.table(outer_method_frequency, "model.nestedCV_algorithm_frequency.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(c(
  paste0("Final selected model: ", best_method),
  "Unbiased internal performance: model.nestedCV_performance.txt.",
  "Outer folds were used only to evaluate the complete model‑selection pipeline.",
  "Within each outer training set, inner CV repeated preprocessing, feature selection, tuning and algorithm comparison.",
  "The final algorithm was selected by CV using the complete training set only.",
  "External validation cohorts were used only after preprocessing, algorithm, features and thresholds were frozen.",
  "External AUCs were not used for model selection or ranking.",
  "No random noise was added. Validation data used training median/mean/SD only."
), "model.selectionRule.txt")

all_errors <- safe_rbind(c(outer_error_log,
                           list(FullTrainingInner = final_inner$Errors), final_fit_errors))
if (!is.null(all_errors))
  write.table(all_errors, "model.error_log.txt", sep = "\t",
              row.names = FALSE, quote = FALSE)

if (nrow(AUC_mat) > 0 && ncol(AUC_mat) > 0) {
  CohortCol <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(ncol(AUC_mat)),
                        colnames(AUC_mat))
  hm <- SimpleHeatmap(Cindex_mat = AUC_mat,
                      avg_Cindex = AUC_mat[, "Train_InnerCV"],
                      CohortCol = CohortCol, barCol = "#4C78A8",
                      cellwidth = 1, cellheight = 0.5,
                      cluster_columns = FALSE, cluster_rows = FALSE)
  pdf("model.AUCheatmap.pdf", width = ncol(AUC_mat) + 6,
      height = max(5, 0.5 * nrow(AUC_mat) * 0.45 + 2))
  draw(hm); dev.off()
}

cat("Strict Nested‑CV pipeline finished.\n")
cat("Unbiased internal AUC: ", unname(nested_metric["AUC"]), "\n", sep="")
cat("Final algorithm selected by full‑training‑set inner CV: ", best_method, "\n", sep="")
cat("External validation was never involved in preprocessing, feature‑selection, tuning or algorithm selection.\n")

final_model <- readRDS("model.selectedModel.rds")
