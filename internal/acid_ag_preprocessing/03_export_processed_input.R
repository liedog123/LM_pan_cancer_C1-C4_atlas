#!/usr/bin/env Rscript
# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Task: acid_ag_groups; step: s03
# Inputs (DCF fields): acid_scores,clinical_groups,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_ag_groups --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 03_export_processed_input.R --config PATH [--check-config]\n",
      "       Rscript 03_export_processed_input.R --help\n\n",
      "DCF fields: acid_scores, clinical_groups, tpm_root, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary invocation applies source cohort eligibility and exports an independent processed-input object.\n",
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

# Resolve paths relative to the DCF directory, including existing ancestors and symlinks.
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
required_fields <- c("acid_scores", "clinical_groups", "tpm_root", "output_dir")
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
acid_scores_path <- unname(paths["acid_scores"])
clinical_groups_path <- unname(paths["clinical_groups"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- c(acid_scores_path, clinical_groups_path)
# Inspect TPM path metadata only; no serialized objects are loaded.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
tpm_metadata_files <- list.files(tpm_root, pattern = "_mrna_expr_tpm\\.rdata$", recursive = TRUE, full.names = TRUE)
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM input files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
required_files <- c(required_files, tpm_real_paths)


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
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

# Isolate each serialized input and retrieve only the source's named object.
load_named_object <- function(path, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(path, envir = input_env)
  if (!exists(object_name, envir = input_env, inherits = FALSE)) {
    stop("Required object not found: ", object_name, " in ", path, call. = FALSE)
  }
  get(object_name, envir = input_env, inherits = FALSE)
}
# Adapt only data.table clinical items for the source's data.frame column selection.
as_merge_clinical_tables <- function(clinical_data_cluster_list) {
  lapply(clinical_data_cluster_list, function(x) {
    if (inherits(x, "data.table")) as.data.frame(x) else x
  })
}
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
acid_scores_list <- load_named_object(acid_scores_path, "acid_scores_list")
clinical_data_cluster_list <- load_named_object(clinical_groups_path, "clinical_data_cluster_list")
clinical_data_cluster_list <- as_merge_clinical_tables(clinical_data_cluster_list)

# 获取并预处理 TCGA 表达数据的函数
get_tcga_data_tpm <- function(project) {
  project_dir <- file.path(tpm_root, project)
  if (!dir.exists(project_dir)) {
    stop("项目文件夹不存在：", project_dir)
  }
  tpm_file <- file.path(project_dir, paste0(project, "_mrna_expr_tpm.rdata"))
  if (!file.exists(tpm_file)) {
    stop("TPM 表达数据文件不存在：", tpm_file)
  }
  mrna_expr_tpm <- load_named_object(tpm_file, "mrna_expr_tpm")
  # 标准化样本ID
  colnames(mrna_expr_tpm) <- gsub("\\.", "-", substr(colnames(mrna_expr_tpm), 1, 15))
  # 筛选肿瘤样本：TCGA 样本ID的第14-15位 < 10 代表肿瘤样本
  tumor_samples_indices <- which(as.numeric(substr(colnames(mrna_expr_tpm), 14, 15)) < 10)
  tpm_tumor <- mrna_expr_tpm[, tumor_samples_indices, drop = FALSE]
  return(tpm_tumor)
}

# 筛选样本数 > 30 的癌症类型
eligible_projects <- c()
for (project in names(clinical_data_cluster_list)) {
  clinical_data <- clinical_data_cluster_list[[project]]
  if (!"Group" %in% colnames(clinical_data)) {
    warning(paste("项目", project, "的临床数据中缺少 'Group' 列，已跳过。"))
    next
  }
  group_counts <- table(clinical_data$Group)
  if (all(c("High", "Low") %in% names(group_counts)) &&
      group_counts["High"] > 30 && group_counts["Low"] > 30) {
    eligible_projects <- c(eligible_projects, project)
  }
}

# 需要保留的临床数据列
columns_to_keep <- c(
  "sample", "_PATIENT", "cancer type abbreviation",
  "age_at_initial_pathologic_diagnosis", "gender", "race",
  "ajcc_pathologic_tumor_stage", "clinical_stage",
  "histological_grade", "OS", "OS.time", "DSS", "DSS.time",
  "DFI", "DFI.time", "PFI", "PFI.time", "Group"
)

# 整合数据
pancancer_data_list <- list()
for (project in eligible_projects) {
  message(paste("正在处理项目:", project))

  # 获取当前项目的所有原始数据
  current_clinical <- clinical_data_cluster_list[[project]]
  current_scores <- acid_scores_list[[project]]
  current_expression <- get_tcga_data_tpm(project)

  colnames(current_scores) <- gsub("\\.", "-", substr(colnames(current_scores), 1, 15))

  # 寻找三大数据集共有的样本
  common_samples <- intersect(current_clinical$sample, colnames(current_expression))
  common_samples <- intersect(common_samples, colnames(current_scores))
  common_samples <- sort(common_samples)

  if (length(common_samples) == 0) {
    warning(paste("项目", project, "没有找到任何共同样本，已跳过。"))
    next
  }
  message(paste("...找到", length(common_samples), "个共同样本。"))

  # 根据共同样本过滤并对齐所有数据
  # 过滤临床数据
  clinical_final <- current_clinical[current_clinical$sample %in% common_samples, ]
  clinical_final <- clinical_final[match(common_samples, clinical_final$sample), ]

  # 只保留指定的列
  existing_columns <- intersect(columns_to_keep, colnames(clinical_final))
  clinical_final <- clinical_final[, existing_columns, drop = FALSE]

  # 过滤表达矩阵和得分矩阵
  expression_final <- current_expression[, common_samples, drop = FALSE]
  scores_final <- current_scores[, common_samples, drop = FALSE]

  # 存入最终列表
  project_name_cleaned <- sub("TCGA-", "", project)

  pancancer_data_list[[project_name_cleaned]] <- list(
    clinical_data = clinical_final,
    expression_data = as.matrix(expression_final),
    acid_scores = scores_final
  )
}

save(pancancer_data_list, file = file.path(out_dir, "pancancer_full_data_list.rdata"))
