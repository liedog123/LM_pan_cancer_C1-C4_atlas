#!/usr/bin/env Rscript
# 固定正式12成员，接既有九基因评分计算逐癌Spearman相关。

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 04_correlate_formal_ph_panel.R --config PATH [--check-config]\n",
      "       Rscript 04_correlate_formal_ph_panel.R --help\n\n",
      "DCF fields: acid_scores, tpm_root, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary invocation correlates existing acid scores with the formal 12-gene pH panel; rho output only.\n",
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
required_fields <- c("acid_scores", "tpm_root", "output_dir")
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
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- acid_scores_path
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
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
acid_scores_list <- load_named_object(acid_scores_path, "acid_scores_list")
projects <- expected_projects

# 正式成员由作者显式固定；源文件随后覆盖为7成员的语句不纳入公开副本。
phGene <- c("PDK1", "ACSS2", "DHCR7", "ELOVL5", "FDFT1", "HMGCR",
            "HMGCS1", "IDI1", "LDLR", "MSMO1", "MTHFR", "SQLE")

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

#计算酸性得分和ph基因的相关性
# 初始化一个空矩阵来存储所有癌症类型的相关性结果
correlation_matrix <- matrix(NA, nrow = length(phGene), ncol = length(projects))
rownames(correlation_matrix) <- phGene
colnames(correlation_matrix) <- projects

# 遍历每个癌症类型
for (project in projects) {
  # 提取数据
  expr <- get_tcga_data_tpm(project)
  colnames(expr) <- gsub("\\.", "-", substr(colnames(expr), 1, 15))
  acid_scores <- acid_scores_list[[project]]

  # 确保表达矩阵和酸性评分有共同的样本
  common_samples <- intersect(colnames(expr), colnames(acid_scores))
  expr <- expr[, common_samples, drop = FALSE]
  acid_scores <- acid_scores[, common_samples]

  # 提取 pH 相关基因的表达数据
  ph_matrix <- expr[rownames(expr) %in% phGene, , drop = FALSE]

  # 遍历每个 pH 基因，计算与酸性评分的相关性
  for (i in 1:nrow(ph_matrix)) {
    ph_gene_expression <- as.numeric(ph_matrix[i, ])
    cor_test <- cor.test(acid_scores, ph_gene_expression, method = "spearman")  # 使用 Spearman 相关性
    correlation_matrix[rownames(ph_matrix)[i], project] <- cor_test$estimate    # 存储相关性系数
  }
}
colnames(correlation_matrix) <- sub("^TCGA-", "", colnames(correlation_matrix))


# 仅导出rho矩阵；不保存P值、不新增校正或绘图。
save(correlation_matrix, file = file.path(out_dir, "ph_acid_correlation_matrix.rdata"))
write.table(correlation_matrix, file = file.path(out_dir, "ph_acid_correlation_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
