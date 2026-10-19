#!/usr/bin/env Rscript
# Purpose: Adjust the three reported drug predictions for the other drugs shared PC1.
# Task: drug_shared_axis; step: s04
# Inputs (DCF fields): drug_complete,shared_patient,prediction_audit,prediction_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task drug_shared_axis --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 04_run_shared_axis_adjustment.R --config PATH [--check-config]\n       Rscript 04_run_shared_axis_adjustment.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("drug_complete", "shared_patient", "prediction_audit", "prediction_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("drug_complete", "shared_patient", "prediction_audit")]
required_files <- unname(input_paths)

# 仅验证预测目录元信息，不读取audit内容或枚举预测文件。
prediction_dir <- unname(paths["prediction_dir"])
if (!dir.exists(prediction_dir) || file.access(prediction_dir, 5L) != 0L) stop("Prediction directory must be readable and searchable.")

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),prediction_dir,config_path,package_root)
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

options(stringsAsFactors = FALSE)
output_dir <- out_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  read.delim(
    con, sep = "\t", header = TRUE, quote = "\"", comment.char = "",
    check.names = FALSE, na.strings = c("NA", "")
  )
}

write_tsv <- function(x, name) {
  write.table(
    x, file.path(output_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

cliffs_delta <- function(x, y) {
  n1 <- length(x)
  n0 <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u1 <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n0) - 1
}

# 三药来自当前稿件固定比较集合，不再依赖已退出投稿的内部证据评级。
d <- lapply(input_paths, read_tsv)
strict_drugs <- c("Vinorelbine", "Linsitinib", "Palbociclib")
raw_drug <- d$drug_complete[
  d$drug_complete$drug %in% strict_drugs & d$drug_complete$status == "ok", ,
  drop = FALSE
]
formal_cancers <- sort(unique(raw_drug$Cancer[raw_drug$significant %in% TRUE]))
raw_drug <- raw_drug[raw_drug$Cancer %in% formal_cancers, , drop = FALSE]
if (length(formal_cancers) != 18L || nrow(raw_drug) != 54L) {
  stop("Expected 18 formal cancers and 54 Cancer x Drug units.")
}

# 每个药物使用排除自身后的285药物PC1，避免当前药物直接进入共享轴。
residual_rows <- list()
rr <- 0L
for (cancer in formal_cancers) {
  audit_row <- d$prediction_audit[d$prediction_audit$Cancer == cancer, , drop = FALSE]
  if (nrow(audit_row) != 1L) stop("Prediction source unresolved for ", cancer)
  # 内嵌历史路径只提供文件名；在配置预测目录内解析并保护输出边界。
  source_name <- basename(as.character(audit_row$source_file))
  if (length(source_name) != 1L || is.na(source_name) || !nzchar(source_name) ||
      source_name %in% c(".", "..") || grepl("[/\\\\]", source_name)) {
    stop("Prediction basename unresolved for ", cancer)
  }
  prediction_path <- resolve_path(file.path(prediction_dir, source_name), prediction_dir)
  if (!is_within(prediction_path, prediction_dir) || is_within(prediction_path, out_dir) ||
      !file_test("-f", prediction_path) || file.access(prediction_path, 4L) != 0L ||
      is.na(file.info(prediction_path)$size) || file.info(prediction_path)$size <= 0L) {
    stop("Prediction source must be a protected readable nonempty file for ", cancer)
  }
  pred <- read.csv(
    prediction_path, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_col <- names(pred)[1L]
  pred$Sample15 <- substr(as.character(pred[[sample_col]]), 1L, 15L)
  shared <- d$shared_patient[d$shared_patient$Cancer == cancer, , drop = FALSE]
  idx <- match(shared$Sample15, pred$Sample15)
  if (anyNA(idx)) stop("Shared-score patients missing from prediction matrix: ", cancer)
  pred <- pred[idx, , drop = FALSE]
  if (!identical(pred$Sample15, shared$Sample15)) stop("Sample order mismatch: ", cancer)

  drug_cols <- setdiff(names(pred), c(sample_col, "Sample15"))
  if (!all(strict_drugs %in% drug_cols)) stop("Strict drug column missing: ", cancer)

  cancer_rows <- list()
  for (drug in strict_drugs) {
    other_cols <- setdiff(drug_cols, drug)
    other_mat <- as.matrix(data.frame(lapply(pred[other_cols], as.numeric), check.names = FALSE))
    keep <- apply(other_mat, 2L, function(x) all(is.finite(x)) && stats::sd(x) > 0)
    other_mat <- other_mat[, keep, drop = FALSE]
    if (ncol(other_mat) < 280L) stop("Too few nonconstant leave-one-drug-out columns: ", cancer)
    pca <- prcomp(other_mat, center = TRUE, scale. = TRUE, rank. = 1L)
    pc1 <- pca$x[, 1L]
    pc1_var <- pca$sdev[1L]^2 / sum(pca$sdev^2)

    y <- as.numeric(pred[[drug]])
    fit <- lm(y ~ pc1)
    resid_y <- residuals(fit)
    is_c1 <- shared$C_subtype == "C1"
    x <- resid_y[is_c1]
    y4 <- resid_y[!is_c1]
    wt <- suppressWarnings(wilcox.test(x, y4, exact = FALSE))
    delta <- cliffs_delta(x, y4)
    median_diff <- median(x) - median(y4)
    raw_row <- raw_drug[raw_drug$Cancer == cancer & raw_drug$drug == drug, , drop = FALSE]

    rr <- rr + 1L
    cancer_rows[[drug]] <- data.frame(
      drug = drug,
      Cancer = cancer,
      n_C1 = sum(is_c1),
      n_C4 = sum(!is_c1),
      n_other_drugs_for_PC1 = ncol(other_mat),
      loo_PC1_variance_fraction = pc1_var,
      abs_spearman_drug_vs_loo_PC1 = abs(suppressWarnings(cor(
        as.numeric(pred[[drug]]), pc1, method = "spearman"
      ))),
      raw_cliffs_delta_C1_minus_C4 = raw_row$cliffs_delta_C1_minus_C4,
      raw_more_sensitive_group = raw_row$more_sensitive_group,
      raw_FDR = raw_row$FDR,
      residual_median_diff_C1_minus_C4 = median_diff,
      residual_cliffs_delta_C1_minus_C4 = delta,
      residual_p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
  cancer_df <- do.call(rbind, cancer_rows)
  cancer_df$residual_FDR_within_cancer_3_drugs <- p.adjust(
    cancer_df$residual_p_value, method = "BH"
  )
  cancer_df$residual_more_sensitive_group <- ifelse(
    cancer_df$residual_median_diff_C1_minus_C4 < 0, "C1",
    ifelse(cancer_df$residual_median_diff_C1_minus_C4 > 0, "C4", "tie")
  )
  cancer_df$residual_significant <- (
    cancer_df$residual_FDR_within_cancer_3_drugs < 0.05 &
      abs(cancer_df$residual_cliffs_delta_C1_minus_C4) >= 0.20
  )
  cancer_df$direction_preserved_after_residual <- (
    cancer_df$residual_more_sensitive_group == cancer_df$raw_more_sensitive_group
  )
  residual_rows[[cancer]] <- cancer_df
}
residual <- do.call(rbind, residual_rows)
rownames(residual) <- NULL
write_tsv(residual, "strict_drug_shared_axis_residual_by_cancer.tsv")

if (any(residual$n_other_drugs_for_PC1 < 285L)) stop("Expected all 285 other drug columns.")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
