#!/usr/bin/env Rscript
# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Task: ag_survival; step: s02
# Inputs (DCF fields): cox_results,cancer_order
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task ag_survival --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 02_finalize_ag_cox.R --config PATH [--check-config]\n",
      "       Rscript 02_finalize_ag_cox.R --help\n\n",
      "DCF fields: cox_results, cancer_order, output_dir.\n",
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
required_fields <- c("cox_results", "cancer_order", "output_dir")
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
cox_file <- unname(paths["cox_results"])
cancer_order_file <- unname(paths["cancer_order"])
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

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(readr) })
load_named <- function(path, name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) stop("Missing named object: ", name)
  get(name, envir = env, inherits = FALSE)
}
results_list <- load_named(cox_file, "results_list")
panel_dir <- file.path(out_dir, "panel_exports")
source_data_dir <- file.path(out_dir, "source_data")
source_dir <- source_data_dir  # 兼容源诊断提取段的输出目录变量
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
# 保留完整22癌/两模型家族；图中KM癌种选择不进入Cox/BH。
if (length(results_list) != 22L || any(vapply(results_list, is.null, logical(1L)))) {
  stop("Formal finalization requires all 22 source cancer model results; do not shorten the BH family.")
}

result_cancers <- sub("^TCGA-", "", names(results_list))
shared_order_df <- read_csv(cancer_order_file, show_col_types = FALSE)
cancer_order <- shared_order_df$Cancer[shared_order_df$Cancer %in% result_cancers]
# 仅验证完整输入家族身份/显示顺序；不改变源BH范围。
if (is.null(names(results_list)) || anyNA(names(results_list)) ||
    any(!nzchar(names(results_list))) || anyDuplicated(result_cancers) ||
    anyNA(shared_order_df$Cancer) || anyDuplicated(shared_order_df$Cancer) ||
    length(cancer_order) != 22L || !setequal(cancer_order, result_cancers)) {
  stop("Result cancer names must be unique and fully covered by a unique frozen cancer order.")
}


cox_df <- bind_rows(lapply(names(results_list), function(project) {
  result <- results_list[[project]]
  if (is.null(result)) {
    return(NULL)
  }

  univ_result <- result$univariate
  multiv_result <- result$multivariate

  if ("GroupHigh" %in% rownames(multiv_result$coefficients)) {
    multiv_hr <- exp(multiv_result$coefficients["GroupHigh", "coef"])
    multiv_p <- multiv_result$coefficients["GroupHigh", "Pr(>|z|)"]
  } else {
    multiv_hr <- NA_real_
    multiv_p <- NA_real_
  }

  tibble(
    Cancer_Type = sub("^TCGA-", "", project),
    Model_Type = c("Univariate", "Multivariate"),
    Hazard_Ratio = c(
      univ_result$coefficients[1, "exp(coef)"],
      multiv_hr
    ),
    P_Value = c(
      univ_result$coefficients[1, "Pr(>|z|)"],
      multiv_p
    )
  )
}))

cox_df <- cox_df %>%
  mutate(
    FDR = ave(P_Value, Model_Type, FUN = function(x) p.adjust(x, method = "fdr")),
    Significance = ifelse(FDR <= 0.05, "FDR <= 0.05", "FDR > 0.05"),
    Hazard_Ratio_Discrete = cut(
      Hazard_Ratio,
      breaks = c(-Inf, 0.15, 1.0, 2.0, 3.0, Inf),
      labels = c("<0.15", "0.15-1.0", "1.0-2.0", "2.0-3.0", ">3.0"),
      include.lowest = TRUE
    ),
    Cancer_Type = factor(Cancer_Type, levels = cancer_order),
    Model_Type = factor(Model_Type, levels = c("Multivariate", "Univariate"))
  )

write_csv(
  cox_df,
  file.path(source_data_dir, "Fig1E_cox_summary_values.csv")
)

hr_colors <- c(
  "<0.15" = "#5A88C9",
  "0.15-1.0" = "#B7D8E8",
  "1.0-2.0" = "#F4B6C2",
  "2.0-3.0" = "#F47C7C",
  ">3.0" = "#D7301F"
)

p <- ggplot(
  cox_df,
  aes(x = Cancer_Type, y = Model_Type)
) +
  geom_point(
    aes(fill = Hazard_Ratio_Discrete, size = Significance),
    shape = 22,
    colour = "black",
    stroke = 0.45,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = hr_colors,
    drop = FALSE,
    name = "Hazard ratio"
  ) +
  scale_size_manual(
    values = c("FDR > 0.05" = 4.1, "FDR <= 0.05" = 7.0),
    name = "FDR",
    labels = c("FDR > 0.05" = ">0.05", "FDR <= 0.05" = "<=0.05"),
    guide = guide_legend(
      override.aes = list(shape = 22, fill = "white", colour = "black")
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.6, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(3, 4, 1, 4)
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1, override.aes = list(size = 5.5)),
    size = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(shape = 22, fill = "white", colour = "black"),
      label.hjust = 0.5
    )
  )

ggsave(
  file.path(panel_dir, "Fig1D_cox_summary_v1.pdf"),
  p,
  width = 178 / 25.4,
  height = 55 / 25.4,
  units = "in",
  dpi = 300
)

ggsave(
  file.path(panel_dir, "Fig1D_cox_summary_v1.png"),
  p,
  width = 178 / 25.4,
  height = 55 / 25.4,
  units = "in",
  dpi = 300
)

message("Fig1D Cox summary exported successfully.")

extract_model <- function(project, result, model_type) {
  fit_summary <- if (model_type == "Univariate") result$univariate else result$multivariate
  coef_rows <- rownames(fit_summary$coefficients)
  group_row <- if ("GroupHigh" %in% coef_rows) "GroupHigh" else coef_rows[[1]]
  ci <- fit_summary$conf.int[group_row, , drop = FALSE]
  lower <- unname(ci[1, "lower .95"])
  upper <- unname(ci[1, "upper .95"])
  events <- unname(fit_summary$nevent)
  stability <- case_when(
    !is.finite(lower) | !is.finite(upper) | lower <= 0 ~ "NONFINITE_CI_OR_SEPARATION",
    events < 10 ~ "LOW_EVENT_COUNT_LT10",
    upper / lower > 100 ~ "VERY_WIDE_CI",
    TRUE ~ "OK"
  )
  tibble(
    Cancer_Type = sub("^TCGA-", "", project),
    Model_Type = model_type,
    Extracted_Hazard_Ratio = exp(fit_summary$coefficients[group_row, "coef"]),
    CI95_Lower = lower,
    CI95_Upper = upper,
    Extracted_P_Value = fit_summary$coefficients[group_row, "Pr(>|z|)"],
    Model_N = unname(fit_summary$n),
    Event_N = events,
    Model_Terms = paste(coef_rows, collapse = ";"),
    Stability_Flag = stability
  )
}

d_extract <- bind_rows(lapply(names(results_list), function(project) {
  result <- results_list[[project]]
  if (is.null(result)) return(NULL)
  bind_rows(
    extract_model(project, result, "Univariate"),
    extract_model(project, result, "Multivariate")
  )
}))

d_manifest <- cox_df %>%
  left_join(d_extract, by = c("Cancer_Type", "Model_Type")) %>%
  mutate(
    HR_Match = abs(Hazard_Ratio - Extracted_Hazard_Ratio) < 1e-8 * pmax(1, abs(Hazard_Ratio)),
    P_Match = abs(P_Value - Extracted_P_Value) < 1e-10 * pmax(1, abs(P_Value)),
    Contrast = "High Acid versus Low Acid",
    Endpoint = "overall survival"
  )
stopifnot(nrow(d_manifest) == 44L, all(d_manifest$HR_Match), all(d_manifest$P_Match))
write_tsv(d_manifest, file.path(source_dir, "Fig1E_cox_model_manifest.tsv"))


writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
