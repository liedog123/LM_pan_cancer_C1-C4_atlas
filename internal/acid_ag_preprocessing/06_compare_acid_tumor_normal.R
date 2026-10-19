#!/usr/bin/env Rscript
# Purpose: Summarise Acid_Score across cancers and tumour-normal comparisons.
# Task: acid_distributions; step: s06
# Inputs (DCF fields): validation_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_distributions --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 06_compare_acid_tumor_normal.R --config PATH [--check-config]\n",
      "       Rscript 06_compare_acid_tumor_normal.R --help\n\n",
      "DCF fields: validation_data, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "Compares Acid_Score in the matching prepared tumor-normal object; raw P values and separate median summaries.\n",
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
required_fields <- c("validation_data", "output_dir")
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
required_files <- unname(paths[setdiff(names(paths), "output_dir")])
out_dir <- unname(paths["output_dir"])


# Protect source inputs and existing results; check mode creates no directories.
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
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
full_validation_data <- load_named_object(paths["validation_data"], "full_validation_data")
# 明确截取Fig. 1C Acid_Score；旧五模块及Stage/Grade不属于此入口。
score_cols <- "Acid_Score"

tvn_stats_list <- list()

cancer_list <- unique(full_validation_data$Cancer_Type)

for (cancer in cancer_list) {
  sub_df <- full_validation_data %>% filter(Cancer_Type == cancer)
  
  # 检查是否同时有 T 和 N
  if (length(unique(sub_df$Sample_Type)) < 2) next
  
  # 对每个模块 (包括 Acid Score) 计算差异
  for (mod in score_cols) {
    # Wilcoxon Test
    res <- wilcox.test(sub_df[[mod]] ~ sub_df$Sample_Type) # Normal vs Tumor (default factor order)
    
    # 计算均值差 (Tumor - Normal)
    mean_t <- mean(sub_df[[mod]][sub_df$Sample_Type == "Tumor"], na.rm = TRUE)
    mean_n <- mean(sub_df[[mod]][sub_df$Sample_Type == "Normal"], na.rm = TRUE)
    diff_val <- mean_t - mean_n # 正值表示 Tumor 高
    
    # 存入列表
    tvn_stats_list[[paste(cancer, mod, sep = "_")]] <- data.frame(
      Cancer = cancer,
      Module = mod,
      Diff_Mean = diff_val,
      P_Value = res$p.value
    )
  }
}
# 合并列表
tvn_stats <- bind_rows(tvn_stats_list)

# 标记显著性星号
tvn_stats$Significance <- case_when(
  tvn_stats$P_Value < 0.001 ~ "***",
  tvn_stats$P_Value < 0.01 ~ "**",
  tvn_stats$P_Value < 0.05 ~ "*",
  TRUE ~ ""
)

write.csv(tvn_stats, file.path(out_dir, "Stats_Tumor_vs_Normal.csv"), row.names = FALSE)

tvn_acid <- tvn_stats %>%
  filter(Module == "Acid_Score") %>%
  select(Cancer, Diff_Mean, P_Value, Significance)

# 中位数按源panel数据摘要另表导出；Diff_Mean仍是均值差。
bottom_summary <- full_validation_data %>%
  select(Cancer = Cancer_Type, Sample_Type, Acid_Score) %>%
  group_by(Cancer, Sample_Type) %>%
  summarise(Median_Acid_Score = median(Acid_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample_Type, values_from = Median_Acid_Score) %>%
  left_join(tvn_acid, by = "Cancer")
write.csv(bottom_summary, file.path(out_dir, "Fig1C_bottom_tumor_normal_medians.csv"), row.names = FALSE)
