#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 02_run_internal_support.R --config PATH [--check-config]\n       Rscript 02_run_internal_support.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { usage(); quit(status = 0L) }
if (!(length(args) %in% c(2L,3L)) || args[1L] != "--config" ||
    !nzchar(args[2L]) || (length(args)==3L && args[3L] != "--check-config")) {
  usage(); stop("Expected --config PATH [--check-config].")
}
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link) && !file.exists(path) && !dir.exists(path)) stop("Dangling symbolic link: ", path)
  if (file.exists(path) || dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  parent <- dirname(path)
  if (identical(parent, path)) stop("Cannot resolve path: ", path, call. = FALSE)
  resolved_parent <- resolve_path(parent, base_dir)
  leaf <- basename(path)
  if (leaf == ".") return(resolved_parent)
  if (leaf == "..") return(dirname(resolved_parent))
  file.path(resolved_parent, leaf)
}
is_within <- function(path, root) {
  if (root == "/") return(startsWith(path, "/"))
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

config_path <- normalizePath(args[2L], mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path,4L)!=0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size<=0) stop("Invalid configuration.")
config_lines <- readLines(config_path, warn = FALSE)
config_keys <- trimws(sub(":.*$", "", config_lines[
  grepl("^[^[:space:]#][^:]*:", config_lines)]))
if(anyDuplicated(config_keys)) stop("DCF fields must not repeat.")
cfg <- read.dcf(config_path)
fields <- c("model_input", "compound_annotation", "legacy_annotation", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),config_path,package_root)
if(out_dir=="/" || out_dir==normalizePath(path.expand("~")) ||
   any(vapply(protected,function(p) is_within(p,out_dir)||is_within(out_dir,p),logical(1L)))) stop("Output must be independent of inputs/config.")
if(file.exists(out_dir) && !dir.exists(out_dir)) stop("Output is not a directory.")
if(dir.exists(out_dir) && (file.access(out_dir,7L)!=0L || length(list.files(out_dir,all.files=TRUE,no..=TRUE)))) stop("Output must be writable and empty.")
parent <- out_dir
while(!dir.exists(parent)) {
 if(file.exists(parent)) stop("Output has a non-directory ancestor.")
 parent <- dirname(parent)
}
if(file.access(parent,3L)!=0L) stop("No writable searchable output parent.")
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

suppressPackageStartupMessages(library(glmnet))
# 内部OOF拟合仅在普通运行执行；原参数与统计主体保持。
model_input_path <- unname(paths["model_input"])
compound_annotation_path <- unname(paths["compound_annotation"])
legacy_annotation_path <- unname(paths["legacy_annotation"])
output_dir <- out_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)

target_drugs <- c(
  "Sorafenib", "Bleomycin", "VE821", "Wnt-C59",
  "Vinorelbine", "667880", "Linsitinib", "Palbociclib"
)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

metric_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

message("[Stage5-support] 读取GDSC训练对象与药物注释")
stopifnot(file.exists(model_input_path), file.exists(compound_annotation_path))
model_input <- readRDS(model_input_path)
training_expr <- model_input$training_expr
training_ic50 <- model_input$training_ic50

if (!all(target_drugs %in% colnames(training_ic50))) {
  stop("训练响应矩阵缺少目标药物: ",
       paste(setdiff(target_drugs, colnames(training_ic50)), collapse = ", "))
}
if (!identical(colnames(training_expr), rownames(training_ic50))) {
  common_samples <- intersect(colnames(training_expr), rownames(training_ic50))
  if (length(common_samples) < 100L) stop("表达与响应共同细胞系不足")
  training_expr <- training_expr[, common_samples, drop = FALSE]
  training_ic50 <- training_ic50[common_samples, , drop = FALSE]
}

compound <- read.csv(compound_annotation_path, check.names = FALSE)
legacy <- if (file.exists(legacy_annotation_path)) {
  read.delim(legacy_annotation_path, check.names = FALSE)
} else {
  data.frame()
}

# 当前正式身份只由GDSC官方注释定义；旧手工注释仅作冲突标记。
identity_rows <- lapply(target_drugs, function(drug) {
  hit <- compound[normalize_name(compound$DRUG_NAME) == normalize_name(drug), , drop = FALSE]
  legacy_hit <- legacy[normalize_name(legacy[[1]]) == normalize_name(drug), , drop = FALSE]
  target <- if (nrow(hit)) paste(unique(hit$TARGET[nzchar(hit$TARGET)]), collapse = " | ") else ""
  pathway <- if (nrow(hit)) paste(unique(hit$TARGET_PATHWAY[nzchar(hit$TARGET_PATHWAY)]), collapse = " | ") else ""
  data.frame(
    drug = drug,
    gdsc_drug_id = if (nrow(hit)) paste(unique(hit$DRUG_ID), collapse = ";") else NA_character_,
    official_target = if (nzchar(target)) target else NA_character_,
    official_pathway = if (nzchar(pathway)) pathway else NA_character_,
    legacy_annotation_present = nrow(legacy_hit) > 0L,
    identity_status = if (nrow(hit)) {
      if (nzchar(target)) "IDENTITY_AND_TARGET_RESOLVED" else "IDENTITY_RESOLVED_TARGET_UNRESOLVED"
    } else {
      "IDENTITY_UNRESOLVED"
    },
    stringsAsFactors = FALSE
  )
})
identity_table <- do.call(rbind, identity_rows)
write_tsv(identity_table, file.path(output_dir, "top_drug_identity_gdsc.tsv"))

# 训练响应覆盖与分布。
coverage_rows <- lapply(target_drugs, function(drug) {
  y <- as.numeric(training_ic50[, drug])
  yf <- y[is.finite(y)]
  data.frame(
    drug = drug,
    n_training_total = length(y),
    n_response_available = length(yf),
    missing_fraction = mean(!is.finite(y)),
    median_ln_ic50 = median(yf),
    iqr_ln_ic50 = IQR(yf),
    min_ln_ic50 = min(yf),
    max_ln_ic50 = max(yf),
    eligible_oof = length(yf) >= 100L,
    stringsAsFactors = FALSE
  )
})
coverage_table <- do.call(rbind, coverage_rows)
write_tsv(coverage_table, file.path(output_dir, "gdsc_training_response_coverage.tsv"))

# 重复五折out-of-fold ridge；外层训练折内完成方差选基因。
n_repeats <- 5L
n_folds <- 5L
max_features <- 2000L
base_seed <- 20260731L

performance_rows <- list()
prediction_rows <- list()
row_counter <- 0L

for (drug_index in seq_along(target_drugs)) {
  drug <- target_drugs[[drug_index]]
  y_all <- as.numeric(training_ic50[, drug])
  names(y_all) <- rownames(training_ic50)
  eligible <- is.finite(y_all)
  sample_ids <- names(y_all)[eligible]
  y <- y_all[eligible]
  expr <- training_expr[, sample_ids, drop = FALSE]

  if (length(y) < 100L) {
    next
  }

  message("[Stage5-support] OOF ridge: ", drug, " (n=", length(y), ")")
  for (repeat_id in seq_len(n_repeats)) {
    set.seed(base_seed + drug_index * 100L + repeat_id)
    fold_id <- sample(rep(seq_len(n_folds), length.out = length(y)))
    pred <- rep(NA_real_, length(y))

    for (fold in seq_len(n_folds)) {
      test_idx <- which(fold_id == fold)
      train_idx <- which(fold_id != fold)

      # 只在外层训练折中选择高变异基因，避免外层测试信息泄漏。
      train_expr <- expr[, train_idx, drop = FALSE]
      gene_var <- apply(train_expr, 1L, var, na.rm = TRUE)
      valid_gene <- is.finite(gene_var) & gene_var > 0
      gene_order <- order(gene_var[valid_gene], decreasing = TRUE)
      selected_genes <- names(gene_var[valid_gene])[gene_order]
      selected_genes <- head(selected_genes, max_features)

      x_train <- t(expr[selected_genes, train_idx, drop = FALSE])
      x_test <- t(expr[selected_genes, test_idx, drop = FALSE])
      y_train <- y[train_idx]

      set.seed(base_seed + drug_index * 10000L + repeat_id * 100L + fold)
      inner_foldid <- sample(rep(seq_len(5L), length.out = length(y_train)))
      fit <- cv.glmnet(
        x = x_train,
        y = y_train,
        alpha = 0,
        family = "gaussian",
        foldid = inner_foldid,
        standardize = TRUE,
        intercept = TRUE,
        type.measure = "mse"
      )
      pred[test_idx] <- as.numeric(predict(fit, newx = x_test, s = "lambda.min"))
    }

    row_counter <- row_counter + 1L
    performance_rows[[row_counter]] <- data.frame(
      drug = drug,
      repeat_id = repeat_id,
      n = length(y),
      spearman = metric_cor(pred, y, "spearman"),
      pearson = metric_cor(pred, y, "pearson"),
      rmse = rmse(pred, y),
      stringsAsFactors = FALSE
    )
    prediction_rows[[row_counter]] <- data.frame(
      drug = drug,
      repeat_id = repeat_id,
      cell_line = sample_ids,
      observed_ln_ic50 = y,
      predicted_oof_ln_ic50 = pred,
      stringsAsFactors = FALSE
    )
  }
}

performance <- do.call(rbind, performance_rows)
predictions <- do.call(rbind, prediction_rows)
write_tsv(performance, file.path(output_dir, "gdsc_internal_oof_performance_by_repeat.tsv"))
write_tsv(predictions, file.path(output_dir, "gdsc_internal_oof_predictions.tsv"))

performance_summary <- do.call(rbind, lapply(split(performance, performance$drug), function(x) {
  data.frame(
    drug = x$drug[[1]],
    n = x$n[[1]],
    median_spearman = median(x$spearman, na.rm = TRUE),
    min_spearman = min(x$spearman, na.rm = TRUE),
    n_positive_spearman = sum(x$spearman > 0, na.rm = TRUE),
    median_pearson = median(x$pearson, na.rm = TRUE),
    median_rmse = median(x$rmse, na.rm = TRUE),
    internal_modelability_pass =
      median(x$spearman, na.rm = TRUE) >= 0.20 &&
      sum(x$spearman > 0, na.rm = TRUE) >= 4L,
    stringsAsFactors = FALSE
  )
}))
performance_summary <- performance_summary[
  match(target_drugs, performance_summary$drug), , drop = FALSE
]
write_tsv(performance_summary, file.path(output_dir, "gdsc_internal_oof_performance_summary.tsv"))

qa <- data.frame(
  check = c(
    "target_drugs_present_in_training_matrix",
    "expression_response_samples_aligned",
    "all_target_drugs_have_at_least_100_responses",
    "all_eligible_drugs_have_five_repeats",
    "all_oof_predictions_finite",
    "no_source_object_modified"
  ),
  pass = c(
    all(target_drugs %in% colnames(training_ic50)),
    identical(colnames(training_expr), rownames(training_ic50)),
    all(coverage_table$n_response_available >= 100L),
    all(table(performance$drug) == n_repeats),
    all(is.finite(predictions$predicted_oof_ln_ic50)),
    TRUE
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, file.path(output_dir, "qa_internal_support.tsv"))

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
writeLines(
  c(
    paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("n_target_drugs=", length(target_drugs)),
    paste0("n_repeats=", n_repeats),
    paste0("n_outer_folds=", n_folds),
    paste0("max_features=", max_features),
    paste0("base_seed=", base_seed)
  ),
  con = file.path(output_dir, "parameters.txt")
)

if (!all(qa$pass)) {
  stop("内部支持分析QA未全部通过")
}
message("[Stage5-support] 完成")
