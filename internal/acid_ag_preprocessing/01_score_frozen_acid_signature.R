#!/usr/bin/env Rscript
# Continue acid scoring, AG grouping or processed-input export from fixed study inputs.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 01_score_frozen_acid_signature.R --config PATH [--check-config]\n",
      "       Rscript 01_score_frozen_acid_signature.R --help\n\n",
      "DCF fields: final_genes, clinical_table, tpm_root, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary invocation computes tumor acid scores from the frozen signature and writes clinical/score objects.\n",
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
required_fields <- c("final_genes", "clinical_table", "tpm_root", "output_dir")
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
final_genes_path <- unname(paths["final_genes"])
clinical_table_path <- unname(paths["clinical_table"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- c(final_genes_path, clinical_table_path)
# Inspect TPM path metadata only; no serialized objects are loaded.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
expected_projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP",  "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)
tpm_metadata_files <- file.path(tpm_root, expected_projects, paste0(expected_projects, "_mrna_expr_tpm.rdata"))
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
suppressPackageStartupMessages({
  library(GSVA)
  library(data.table)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
final_selected_genes <- load_named_object(final_genes_path, "final_selected_genes")

acid_gene_set <- list( "acid_genes" = as.matrix(final_selected_genes$Gene))   # 定义酸性基因集
#获取临床数据
all_clinical_data <- fread(clinical_table_path, sep = "\t", header = TRUE)
all_clinical_data <- all_clinical_data[!duplicated(all_clinical_data$sample), ]

#选择癌症类型
projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP",  "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)

# 获取 TCGA 数据的函数
get_tcga_data_tpm <- function(project) {
  project_dir <- file.path(tpm_root, project)
  if (!dir.exists(project_dir)) {
    stop("项目文件夹不存在：", project_dir)
  }
  message("开始处理项目 ", project)
  tpm_file <- file.path(project_dir, paste0(project, "_mrna_expr_tpm.rdata"))
  if (file.exists(tpm_file)) {
    mrna_expr_tpm <- load_named_object(tpm_file, "mrna_expr_tpm")
  } else {
    stop("TPM 表达数据文件不存在：", tpm_file)
  }
  return(mrna_expr_tpm)
}

#整理临床数据和表达数据
common_sample_select <- function(all_clinical_data){
  tpm_matrix <- get_tcga_data_tpm(project)
  colnames(tpm_matrix) <- gsub("\\.", "-", substr(colnames(tpm_matrix), 1, 15))

  tumor_samples <- which(as.numeric(substr(colnames(tpm_matrix), 14, 15)) < 10)
  tpm_matrix <- tpm_matrix[, tumor_samples]

  common_samples <- intersect(all_clinical_data$sample, colnames(tpm_matrix))
  clinical_data <- all_clinical_data[all_clinical_data$sample %in% common_samples, ]

  return(clinical_data)
}

clinical_data_list <- list()
for (project in projects) {
  clinical_data <- common_sample_select(all_clinical_data)
  clinical_data_list[[project]] <- clinical_data
}
#保留含有表达矩阵的临床样本
save(clinical_data_list, file = file.path(out_dir, "clinical_data_list.rdata"))

#遍历癌症类型计算酸性得分
acid_scores_list <- list()
for (project in projects) {
  # 获取数据
  tpm_matrix <- get_tcga_data_tpm(project)
  colnames(tpm_matrix) <- gsub("\\.", "-", substr(colnames(tpm_matrix), 1, 15))
  all_clinical_data <- clinical_data_list[[project]]

  # 筛选样本
  tumor_samples <- which(as.numeric(substr(colnames(tpm_matrix), 14, 15)) < 10)
  tpm_matrix <- tpm_matrix[, tumor_samples]
  common_samples <- intersect(all_clinical_data$sample, colnames(tpm_matrix))
  tpm_matrix <- tpm_matrix[, common_samples]
  tpm_matrix <- log2(tpm_matrix+1)

  # ssgsea分析：酸性基因集得分
  ssgsea_param <- ssgseaParam(exprData = as.matrix(tpm_matrix),
                              geneSets = acid_gene_set,
                              normalize = TRUE)
  gsva_results <- gsva(ssgsea_param, verbose = FALSE)
  acid_scores_list[[project]] <- gsva_results
}

save(acid_scores_list, file = file.path(out_dir, "per_acid_scores_list.rdata"))
