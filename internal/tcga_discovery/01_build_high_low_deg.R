#!/usr/bin/env Rscript
# Purpose: Identify acidification-associated immune pathway enrichment in TCGA.
# Task: tcga_immune_discovery; step: s01
# Inputs (DCF fields): processed_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task tcga_immune_discovery --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 01_build_high_low_deg.R --config PATH [--check-config]\n",
      "       Rscript 01_build_high_low_deg.R --help\n\n",
      "DCF fields: processed_data, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation calculates per-cancer High versus Low differential expression.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (!(length(args) %in% c(2L, 3L)) || args[1L] != "--config" || !nzchar(args[2L]) ||
    (length(args) == 3L && args[3L] != "--check-config")) {
  usage()
  stop("Expected --config PATH, optionally followed by --check-config.", call. = FALSE)
}
flags <- if (length(args) == 3L) args[3L] else character()

# Resolve DCF-relative paths through existing ancestors and symlinks.
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
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
config_path <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path, mode = 4L) != 0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size <= 0) {
  stop("Configuration must be a readable nonempty regular file.", call. = FALSE)
}
cfg <- read.dcf(config_path)
required_fields <- c("processed_data", "output_dir")
if (nrow(cfg) != 1L || !setequal(colnames(cfg), required_fields) ||
    ncol(cfg) != length(required_fields)) {
  stop("Configuration must contain one DCF record with exactly: ",
       paste(required_fields, collapse = ", "), call. = FALSE)
}
cfg <- setNames(as.character(cfg[1L, required_fields]), required_fields)
if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME", cfg, fixed = TRUE))) {
  stop("Replace all example placeholders with explicit paths.", call. = FALSE)
}
paths <- vapply(cfg, resolve_path, character(1L), base_dir = dirname(config_path))
processed_data <- unname(paths["processed_data"])
out_dir <- unname(paths["output_dir"])
required_files <- processed_data

# Protect inputs and existing outputs; checks create no directories.
protected_paths <- c(dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
if ("--check-config" %in% flags) {
  cat("Configuration paths passed metadata checks. No inputs loaded; no analysis or outputs created.\n")
  quit(status = 0L)
}

# Retrieve only the named input object; preserve configured variables.
load_named_object <- function(path, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(path, envir = input_env)
  if (!exists(object_name, envir = input_env, inherits = FALSE)) {
    stop("Missing input object: ", object_name, call. = FALSE)
  }
  get(object_name, envir = input_env, inherits = FALSE)
}
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})
pancancer_data_list <- load_named_object(processed_data, "pancancer_data_list")
main_output_dir <- out_dir
full_results_dir <- file.path(main_output_dir, "full_degs_for_gsea")
sig_results_dir <- file.path(main_output_dir, "significant_degs")
for (path in c(out_dir, full_results_dir, sig_results_dir)) {
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(path)) {
    stop("Cannot create output directory.", call. = FALSE)
  }
}

cancer_types <- names(pancancer_data_list)
for (cancer in cancer_types) {
  message(paste("正在处理:", cancer))

  # 提取当前癌种的数据
  current_data <- pancancer_data_list[[cancer]]
  expr_matrix_full <- current_data$expression_data
  clinical_info_full <- current_data$clinical_data

  if (is.null(expr_matrix_full) || is.null(clinical_info_full) || nrow(clinical_info_full) == 0 || ncol(expr_matrix_full) == 0) {
    warning(paste("跳过", cancer, "，因为数据不完整。"))
    next
  }

  # 筛选出仅包含 "High" 和 "Low" 组的样本
  clinical_info <- clinical_info_full %>%
    filter(Group %in% c("High", "Low"))
  if(nrow(clinical_info) < 2 || length(unique(clinical_info$Group)) < 2) {
    warning(paste("跳过", cancer, "，因为过滤后High/Low组样本不足。"))
    next
  }
  expr_matrix <- expr_matrix_full[, clinical_info$sample]

  # 创建分组信息
  group_list <- factor(clinical_info$Group, levels = c("Low", "High"))
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- c("Low", "High")

  message("...正在进行差异分析")
  v <- voom(expr_matrix + 1, design, plot = FALSE)
  fit <- lmFit(v, design)
  contrast_matrix <- makeContrasts(High_vs_Low = High - Low, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  full_results <- topTable(fit2, number = Inf)
  full_results$gene <- rownames(full_results)

  output_path_full <- file.path(full_results_dir,
                                paste0(cancer, "_DEG_full_results.csv"))
  write.csv(full_results, file = output_path_full, row.names = FALSE)

  # 保存显著差异基因
  adj_p_value_threshold <- 0.05
  logfc_threshold <- 1.0

  significant_degs <- full_results %>%
    filter(adj.P.Val < adj_p_value_threshold, abs(logFC) > logfc_threshold)

  if (nrow(significant_degs) > 0) {
    output_path_sig <- file.path(sig_results_dir,
                                 paste0(cancer, "_DEG_significant.csv"))
    write.csv(significant_degs, file = output_path_sig, row.names = FALSE)
    message(paste("...处理完成！找到", nrow(significant_degs), "个显著差异基因。"))
  } else {
    message("...处理完成！未找到符合条件的显著差异基因。")
  }
}
