#!/usr/bin/env Rscript
# Continue source AG1 versus AG3 overall-survival analyses from prepared inputs.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 00_prepare_ag_km.R --config PATH [--check-config]\n",
      "       Rscript 00_prepare_ag_km.R --help\n\n",
      "DCF fields: projects_data, clinical_groups, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation executes the documented AG continuation.\n",
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
required_fields <- c("projects_data", "clinical_groups", "output_dir")
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
projects_file <- unname(paths["projects_data"])
clinical_file <- unname(paths["clinical_groups"])
out_dir <- unname(paths["output_dir"])
required_files <- unname(paths[setdiff(required_fields, "output_dir")])
# Resolve each input's final target before checking output independence.
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)

# Protect inputs and existing results; checks create no directories.
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

suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2); library(dplyr); library(ggpubr)
})
load_named <- function(path, name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) stop("Missing named object: ", name)
  get(name, envir = env, inherits = FALSE)
}
projects <- load_named(projects_file, "projects")
clinical_data_cluster_list <- load_named(clinical_file, "clinical_data_cluster_list")
dir.create(file.path(out_dir, "KM"), recursive = TRUE, showWarnings = FALSE)

# 初始化列表以存储每个癌症类型的生存数据和显著的癌症类型
survival_data_list <- list()
significant_projects <- c()

# 主循环进行数据准备和生存分析
for (project in projects) {
  # 筛选数据，只保留 High 和 Low 组
  survival_data <- clinical_data_cluster_list[[project]]
  survival_data <- survival_data[survival_data$Group %in% c("High", "Low"), ]
  
  # 存储生存数据到列表中
  survival_data_list[[project]] <- survival_data
  
  # 生存分析
  surv_object <- Surv(time = survival_data$OS.time, event = survival_data$OS)
  fit <- survfit(surv_object ~ Group, data = survival_data)
  
  # 计算 p 值并检查显著性
  pval <- surv_pvalue(fit)$pval
  if (pval < 0.05) {
    significant_projects <- c(significant_projects, project)
  }
  
  # 去除 "TCGA-"前缀
  project <- sub("TCGA-", "", project)
  
  # 绘制生存曲线
  p <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,                         # 显示 p 值
    legend.title = "Group",               # 图例标题
    legend.labs = c("High", "Low"),      # 图例标签
    palette = c("red", "blue"), 
    title = paste(project),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    ggtheme = theme_classic()            # 使用白色背景
  )
  
  
  # 确保保存生存曲线的文件夹存在
  km_dir <- file.path(out_dir, "KM")
  
  # 保存生存曲线
  ggsave(filename = file.path(km_dir, paste0(project, "_KM_plot.png")), plot = p$plot, width = 8, height = 6, dpi = 300)
}

# 将生存数据列表和显著的癌症类型保存为 Rdata 文件
save(survival_data_list, file = file.path(out_dir, "tow_group_clinical_data_list.rdata"))
save(significant_projects, file = file.path(out_dir, "significant_projects.rdata"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))


