#!/usr/bin/env Rscript
# Continue source AG1 versus AG3 overall-survival analyses from prepared inputs.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 01_run_ag_cox.R --config PATH [--check-config]\n",
      "       Rscript 01_run_ag_cox.R --help\n\n",
      "DCF fields: projects_data, survival_data, output_dir.\n",
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
required_fields <- c("projects_data", "survival_data", "output_dir")
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
survival_file <- unname(paths["survival_data"])
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

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival) })
load_named <- function(path, name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) stop("Missing named object: ", name)
  get(name, envir = env, inherits = FALSE)
}
projects <- load_named(projects_file, "projects")
survival_data_list <- load_named(survival_file, "survival_data_list")
# 原函数保留未使用的gene_list参数；不加载与模型无关的基因对象。
gene_list <- NULL
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 定义用于 Cox 回归分析的函数
run_cox_analysis <- function(project, gene_list) {
  # 准备临床数据
  clinical_data <- survival_data_list[[project]]
  clinical_data$Group <- as.factor(clinical_data$Group)
  clinical_data$Group <- relevel(clinical_data$Group, ref = "Low")
  
  # 单因素 Cox 回归
  univ_model <- coxph(Surv(OS.time, OS) ~ Group, data = clinical_data, control = coxph.control(iter.max = 300))
  univ_summary <- summary(univ_model)
  
  # 多因素 Cox 回归数据处理
  # 处理 gender 变量，并转化为因子
  if ("gender" %in% colnames(clinical_data)) {
    unique_genders <- unique(clinical_data$gender)
    if (!("MALE" %in% unique_genders && "FEMALE" %in% unique_genders)) {
      clinical_data <- clinical_data %>% select(-gender)
    } else {
      clinical_data$gender <- as.factor(clinical_data$gender)
      #设置female为参考值
      clinical_data$gender <- relevel(factor(clinical_data$gender), ref = "FEMALE")
    }
  }
  
  # 处理 race 变量
  if ("race" %in% colnames(clinical_data)) {
    clinical_data <- clinical_data %>%
      mutate(
        race = case_when(
          race == "WHITE" ~ "WHITE",
          TRUE ~ "other"
        )
      )
    # 转换为因子并设置 "white" 为参考值
    clinical_data$race <- factor(clinical_data$race, levels = c("WHITE", "other"))
    clinical_data$race <- relevel(clinical_data$race, ref = "WHITE")
    # 检查 race 是否只包含一个水平
    if (length(unique(clinical_data$race)) == 1) {
      clinical_data <- clinical_data %>% select(-race)
    }
  }
  
  
  # 删除 age_at_initial_pathologic_diagnosis 中包含 NA 的行
  if ("age_at_initial_pathologic_diagnosis" %in% colnames(clinical_data)) {
    clinical_data <- clinical_data %>%
      filter(!is.na(age_at_initial_pathologic_diagnosis))
  }
  
  summary <- summary(clinical_data)
  
  # 定义协变量
  available_vars <- c("Group", "race", "gender", "age_at_initial_pathologic_diagnosis")
  available_covariates <- intersect(available_vars, colnames(clinical_data))
  
  # 动态构建 Cox 回归公式
  cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(available_covariates, collapse = " + ")))
  
  # 拟合模型
  multiv_model <- coxph(cox_formula, data = clinical_data, control = coxph.control(iter.max = 300))
  multiv_summary <- summary(multiv_model)
  
  # 返回结果
  list(
    univariate = univ_summary,
    multivariate = multiv_summary,
    data_summary = summary
  )
}


# 保存所有结果的列表
results_list <- list()
# 保存每个项目的警告信息
warnings_list <- list()

for (project in projects) {
  warnings <- NULL
  result <- withCallingHandlers(
    tryCatch(
      run_cox_analysis(project, gene_list),
      error = function(e) {
        message(paste("Error in project", project, ":", e$message))
        NULL
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")  # 取消后续显示同一警告
    }
  )
  
  # 将结果存入结果列表
  results_list[[project]] <- result
  # 将警告信息存入警告列表
  if (!is.null(warnings)) {
    warnings_list[[project]] <- warnings
  }
}


# 保存结果到文件中
save(results_list, file = file.path(out_dir, "Cox_analysis_results.RData"))
save(warnings_list, file = file.path(out_dir, "Cox_warnings.RData"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))


