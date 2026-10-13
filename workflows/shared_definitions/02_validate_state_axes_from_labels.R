#!/usr/bin/env Rscript
# Continue state-axis validation or canonical-label export from upstream tables.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 02_validate_state_axes_from_labels.R --config PATH [--check-config]\n",
      "       Rscript 02_validate_state_axes_from_labels.R --help\n\n",
      "DCF fields: assignment_labels, feature_metadata, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation calculates state-axis statistics and writes outputs.\n",
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

# 相对路径以 DCF 所在目录为基准；既有祖先和符号链接均解析为实际路径。
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
required_fields <- c("assignment_labels", "feature_metadata", "output_dir")
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
label_file <- unname(paths["assignment_labels"])
feature_metadata_file <- unname(paths["feature_metadata"])
out_dir <- unname(paths["output_dir"])
required_files <- c(label_file, feature_metadata_file)

# 防止输出覆盖输入、配置或既有结果；检查模式不会创建目录。
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
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

c_levels <- c("C1", "C2", "C3", "C4")
comparisons <- tibble(
  comparison = c("C1_vs_C2", "C3_vs_C4", "C2_vs_C3", "C1_vs_C4", "C1_vs_C3", "C2_vs_C4"),
  group_1 = c("C1", "C3", "C2", "C1", "C1", "C2"),
  group_2 = c("C2", "C4", "C3", "C4", "C3", "C4"),
  comparison_class = c(
    "high_acid_inflammatory_axis",
    "low_acid_inflammatory_axis",
    "decoupled_intermediate_states",
    "extreme_endpoint",
    "reference",
    "reference"
  )
)

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

safe_quantile <- function(x, p) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

kruskal_one <- function(df, score_col) {
  df <- df %>% filter(!is.na(.data[[score_col]]), !is.na(C_subtype))
  n <- nrow(df)
  n_groups <- n_distinct(df$C_subtype)
  min_group_n <- if (n_groups > 0) min(as.integer(table(df$C_subtype))) else NA_integer_
  if (n_groups < 2) {
    return(tibble(
      n = n,
      n_groups = n_groups,
      min_group_n = min_group_n,
      kruskal_statistic = NA_real_,
      kruskal_df = NA_real_,
      kruskal_p = NA_real_,
      epsilon2 = NA_real_
    ))
  }
  kt <- kruskal.test(reformulate("C_subtype", response = score_col), data = df)
  h <- as.numeric(kt$statistic)
  eps2 <- (h - n_groups + 1) / (n - n_groups)
  tibble(
    n = n,
    n_groups = n_groups,
    min_group_n = min_group_n,
    kruskal_statistic = h,
    kruskal_df = as.numeric(kt$parameter),
    kruskal_p = kt$p.value,
    epsilon2 = max(0, min(1, eps2))
  )
}

pairwise_one <- function(df, score_col) {
  bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
    comp <- comparisons[i, ]
    g1 <- comp$group_1
    g2 <- comp$group_2
    x <- df %>% filter(C_subtype == g1) %>% pull(all_of(score_col))
    y <- df %>% filter(C_subtype == g2) %>% pull(all_of(score_col))
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    p <- if (length(x) > 0 && length(y) > 0) {
      tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    tibble(
      comparison = comp$comparison,
      comparison_class = comp$comparison_class,
      group_1 = g1,
      group_2 = g2,
      n_group_1 = length(x),
      n_group_2 = length(y),
      median_group_1 = median(x, na.rm = TRUE),
      median_group_2 = median(y, na.rm = TRUE),
      median_diff_group_1_minus_group_2 = median_group_1 - median_group_2,
      wilcox_p = p
    )
  }))
}

summarise_score_by_cancer <- function(df, score_col) {
  df %>%
    group_by(Cancer, C_subtype, C_label) %>%
    summarise(
      n = sum(!is.na(.data[[score_col]])),
      mean_value = mean(.data[[score_col]], na.rm = TRUE),
      median_value = median(.data[[score_col]], na.rm = TRUE),
      q25_value = safe_quantile(.data[[score_col]], 0.25),
      q75_value = safe_quantile(.data[[score_col]], 0.75),
      .groups = "drop"
    ) %>%
    mutate(score = score_col, .before = 1) %>%
    arrange(score, Cancer, C_subtype)
}

write_metric_set <- function(df, score_col, prefix) {
  summary_by_cancer <- summarise_score_by_cancer(df, score_col)
  write_tsv(summary_by_cancer, file.path(out_dir, paste0(prefix, "_summary_by_cancer.tsv")))

  equalweight_summary <- summary_by_cancer %>%
    group_by(score, C_subtype, C_label) %>%
    summarise(
      n_cancers = n_distinct(Cancer[n > 0]),
      total_n = sum(n),
      mean_of_cancer_medians = mean(median_value, na.rm = TRUE),
      median_of_cancer_medians = median(median_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(score, C_subtype)
  write_tsv(equalweight_summary, file.path(out_dir, paste0(prefix, "_summary_equalweight.tsv")))

  kruskal_by_cancer <- df %>%
    group_by(Cancer) %>%
    group_modify(~ kruskal_one(.x, score_col)) %>%
    ungroup() %>%
    mutate(
      score = score_col,
      kruskal_p_adj_BH = p.adjust(kruskal_p, method = "BH"),
      .before = 1
    ) %>%
    arrange(kruskal_p)
  write_tsv(kruskal_by_cancer, file.path(out_dir, paste0(prefix, "_kruskal_by_cancer.tsv")))

  pairwise_by_cancer <- df %>%
    group_by(Cancer) %>%
    group_modify(~ pairwise_one(.x, score_col)) %>%
    ungroup() %>%
    group_by(Cancer) %>%
    mutate(wilcox_p_adj_within_cancer_BH = p.adjust(wilcox_p, method = "BH")) %>%
    ungroup() %>%
    mutate(score = score_col, .before = 1) %>%
    arrange(Cancer, comparison)
  write_tsv(pairwise_by_cancer, file.path(out_dir, paste0(prefix, "_pairwise_by_cancer.tsv")))

  pairwise_summary <- pairwise_by_cancer %>%
    group_by(score, comparison, comparison_class, group_1, group_2) %>%
    summarise(
      n_cancers_tested = sum(!is.na(wilcox_p)),
      n_p_lt_0.05 = sum(wilcox_p < 0.05, na.rm = TRUE),
      n_p_adj_lt_0.05 = sum(wilcox_p_adj_within_cancer_BH < 0.05, na.rm = TRUE),
      n_group_1_gt_group_2 = sum(median_diff_group_1_minus_group_2 > 0, na.rm = TRUE),
      n_group_2_gt_group_1 = sum(median_diff_group_1_minus_group_2 < 0, na.rm = TRUE),
      median_of_cancer_median_diff = median(median_diff_group_1_minus_group_2, na.rm = TRUE),
      paired_wilcox_p = {
        diff <- median_diff_group_1_minus_group_2[!is.na(median_diff_group_1_minus_group_2)]
        if (length(diff) > 0) {
          tryCatch(wilcox.test(diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(paired_wilcox_p_adj_BH = p.adjust(paired_wilcox_p, method = "BH")) %>%
    arrange(score, comparison)
  write_tsv(pairwise_summary, file.path(out_dir, paste0(prefix, "_pairwise_summary_by_contrast.tsv")))

  wide <- summary_by_cancer %>%
    select(Cancer, C_subtype, median_value) %>%
    pivot_wider(names_from = C_subtype, values_from = median_value) %>%
    filter(if_all(all_of(c_levels), ~ !is.na(.x))) %>%
    arrange(Cancer)

  paired_pairwise <- bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
    comp <- comparisons[i, ]
    g1 <- comp$group_1
    g2 <- comp$group_2
    diff <- wide[[g1]] - wide[[g2]]
    p <- if (length(diff[!is.na(diff)]) > 0) {
      tryCatch(wilcox.test(diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    tibble(
      score = score_col,
      comparison = comp$comparison,
      comparison_class = comp$comparison_class,
      group_1 = g1,
      group_2 = g2,
      n_cancers = length(diff),
      n_group_1_gt_group_2 = sum(diff > 0, na.rm = TRUE),
      n_group_2_gt_group_1 = sum(diff < 0, na.rm = TRUE),
      median_diff_group_1_minus_group_2 = median(diff, na.rm = TRUE),
      paired_wilcox_p = p
    )
  })) %>%
    mutate(paired_wilcox_p_adj_BH = p.adjust(paired_wilcox_p, method = "BH")) %>%
    arrange(score, comparison)
  write_tsv(paired_pairwise, file.path(out_dir, paste0(prefix, "_pan_cancer_paired_pairwise.tsv")))

  invisible(list(
    summary_by_cancer = summary_by_cancer,
    equalweight_summary = equalweight_summary,
    kruskal_by_cancer = kruskal_by_cancer,
    pairwise_by_cancer = pairwise_by_cancer,
    pairwise_summary = pairwise_summary,
    paired_pairwise = paired_pairwise
  ))
}

message("=== 6.0 C1-C4 axis validation ===")
message("Input labels: ", label_file)
message("Output dir: ", out_dir)

feature_metadata <- read_tsv(feature_metadata_file, show_col_types = FALSE) %>%
  arrange(display_order)
feature_ids <- feature_metadata$feature_id
inflammatory_features <- setdiff(feature_ids, "LM12_Antimicrobial_Peptides")

labels <- read_csv(label_file, show_col_types = FALSE)
required_cols <- c("Sample", "Cancer", "Acid_Score", "AG", "IG", "AGIG", "C_subtype", "C_label", feature_ids)
missing_cols <- setdiff(required_cols, names(labels))
if (length(missing_cols) > 0) {
  stop("标签表缺失必要列: ", paste(missing_cols, collapse = ", "))
}

labels <- labels %>%
  mutate(
    C_subtype = factor(as.character(C_subtype), levels = c_levels),
    AG = factor(as.character(AG), levels = c("AG1", "AG2", "AG3")),
    IG = factor(as.character(IG), levels = c("IG1", "IG2"))
  ) %>%
  group_by(Cancer) %>%
  mutate(
    Acid_Score_z = safe_z(Acid_Score),
    LM12_inflammatory_axis = rowMeans(across(all_of(inflammatory_features)), na.rm = TRUE),
    LM12_all_feature_axis = rowMeans(across(all_of(feature_ids)), na.rm = TRUE),
    LM12_inflammatory_axis_z = safe_z(LM12_inflammatory_axis),
    LM12_all_feature_axis_z = safe_z(LM12_all_feature_axis)
  ) %>%
  ungroup()

feature_z_df <- labels %>%
  select(Sample, Cancer, C_subtype, all_of(feature_ids)) %>%
  group_by(Cancer) %>%
  mutate(across(all_of(feature_ids), safe_z, .names = "{.col}_z")) %>%
  ungroup() %>%
  select(Sample, Cancer, C_subtype, ends_with("_z"))

axis_sample_scores <- labels %>%
  select(
    Sample, Cancer, AG, IG, AGIG, C_subtype, C_label,
    Acid_Score, Acid_Score_z,
    LM12_inflammatory_axis, LM12_inflammatory_axis_z,
    LM12_all_feature_axis, LM12_all_feature_axis_z,
    all_of(feature_ids)
  ) %>%
  left_join(feature_z_df, by = c("Sample", "Cancer", "C_subtype"))
write_tsv(axis_sample_scores, file.path(out_dir, "c1c4_axis_sample_scores.tsv.gz"))

message("[1/4] 样本量与状态中心")
counts_by_cancer <- labels %>%
  count(Cancer, C_subtype, name = "n") %>%
  complete(Cancer, C_subtype = factor(c_levels, levels = c_levels), fill = list(n = 0)) %>%
  group_by(Cancer) %>%
  mutate(
    total_samples = sum(n),
    all_groups_present = all(n > 0),
    all_groups_n_ge_10 = all(n >= 10),
    all_groups_n_ge_20 = all(n >= 20),
    min_group_n = min(n)
  ) %>%
  ungroup() %>%
  arrange(Cancer, C_subtype)
write_tsv(counts_by_cancer, file.path(out_dir, "c1c4_counts_by_cancer_long.tsv"))

centroids_by_cancer <- labels %>%
  group_by(Cancer, C_subtype, C_label) %>%
  summarise(
    n = n(),
    median_acid_z = median(Acid_Score_z, na.rm = TRUE),
    median_inflammatory_axis_z = median(LM12_inflammatory_axis_z, na.rm = TRUE),
    median_all_lm12_axis_z = median(LM12_all_feature_axis_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cancer, C_subtype)
write_tsv(centroids_by_cancer, file.path(out_dir, "c1c4_state_centroids_by_cancer.tsv"))

centroids_equalweight <- centroids_by_cancer %>%
  group_by(C_subtype, C_label) %>%
  summarise(
    n_cancers = n_distinct(Cancer[n > 0]),
    total_n = sum(n),
    mean_median_acid_z = mean(median_acid_z, na.rm = TRUE),
    median_median_acid_z = median(median_acid_z, na.rm = TRUE),
    mean_median_inflammatory_axis_z = mean(median_inflammatory_axis_z, na.rm = TRUE),
    median_median_inflammatory_axis_z = median(median_inflammatory_axis_z, na.rm = TRUE),
    mean_median_all_lm12_axis_z = mean(median_all_lm12_axis_z, na.rm = TRUE),
    median_median_all_lm12_axis_z = median(median_all_lm12_axis_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(C_subtype)
write_tsv(centroids_equalweight, file.path(out_dir, "c1c4_state_centroids_equalweight.tsv"))

message("[2/4] Acid Score 统计")
acid_results <- write_metric_set(labels, "Acid_Score_z", "acid_score_z")

message("[3/4] LM12 inflammatory-remodeling axis 统计")
inflammatory_results <- write_metric_set(labels, "LM12_inflammatory_axis_z", "lm12_inflammatory_axis_z")

message("[4/4] 12 个 LM-Imm feature 统计")
feature_long <- labels %>%
  select(Sample, Cancer, C_subtype, C_label, all_of(feature_ids)) %>%
  pivot_longer(cols = all_of(feature_ids), names_to = "feature_id", values_to = "value") %>%
  group_by(Cancer, feature_id) %>%
  mutate(value_z = safe_z(value)) %>%
  ungroup() %>%
  left_join(feature_metadata, by = "feature_id")

lm12_feature_summary_by_cancer <- feature_long %>%
  group_by(feature_id, short_label, display_order, Cancer, C_subtype, C_label) %>%
  summarise(
    n = sum(!is.na(value_z)),
    mean_value_z = mean(value_z, na.rm = TRUE),
    median_value_z = median(value_z, na.rm = TRUE),
    q25_value_z = safe_quantile(value_z, 0.25),
    q75_value_z = safe_quantile(value_z, 0.75),
    .groups = "drop"
  ) %>%
  arrange(display_order, Cancer, C_subtype)
write_tsv(lm12_feature_summary_by_cancer, file.path(out_dir, "lm12_feature_summary_by_cancer.tsv"))

lm12_feature_summary_equalweight <- lm12_feature_summary_by_cancer %>%
  group_by(feature_id, short_label, display_order, C_subtype, C_label) %>%
  summarise(
    n_cancers = n_distinct(Cancer[n > 0]),
    total_n = sum(n),
    mean_of_cancer_medians_z = mean(median_value_z, na.rm = TRUE),
    median_of_cancer_medians_z = median(median_value_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(display_order, C_subtype)
write_tsv(lm12_feature_summary_equalweight, file.path(out_dir, "lm12_feature_summary_equalweight.tsv"))

lm12_feature_pairwise_by_cancer <- feature_long %>%
  group_by(feature_id, short_label, display_order, Cancer) %>%
  group_modify(~ pairwise_one(.x, "value_z")) %>%
  ungroup() %>%
  group_by(feature_id, Cancer) %>%
  mutate(wilcox_p_adj_within_feature_cancer_BH = p.adjust(wilcox_p, method = "BH")) %>%
  ungroup() %>%
  arrange(display_order, Cancer, comparison)
write_tsv(lm12_feature_pairwise_by_cancer, file.path(out_dir, "lm12_feature_pairwise_by_cancer.tsv"))

lm12_feature_pairwise_summary <- lm12_feature_pairwise_by_cancer %>%
  group_by(feature_id, short_label, display_order, comparison, comparison_class, group_1, group_2) %>%
  summarise(
    n_cancers_tested = sum(!is.na(wilcox_p)),
    n_p_lt_0.05 = sum(wilcox_p < 0.05, na.rm = TRUE),
    n_p_adj_lt_0.05 = sum(wilcox_p_adj_within_feature_cancer_BH < 0.05, na.rm = TRUE),
    n_group_1_gt_group_2 = sum(median_diff_group_1_minus_group_2 > 0, na.rm = TRUE),
    n_group_2_gt_group_1 = sum(median_diff_group_1_minus_group_2 < 0, na.rm = TRUE),
    median_of_cancer_median_diff = median(median_diff_group_1_minus_group_2, na.rm = TRUE),
    paired_wilcox_p = {
      diff <- median_diff_group_1_minus_group_2[!is.na(median_diff_group_1_minus_group_2)]
      if (length(diff) > 0) {
        tryCatch(wilcox.test(diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  group_by(feature_id) %>%
  mutate(paired_wilcox_p_adj_within_feature_BH = p.adjust(paired_wilcox_p, method = "BH")) %>%
  ungroup() %>%
  arrange(display_order, comparison)
write_tsv(lm12_feature_pairwise_summary, file.path(out_dir, "lm12_feature_pairwise_summary_by_contrast.tsv"))

axis_focus_summary <- bind_rows(
  acid_results$pairwise_summary %>%
    mutate(axis = "Acid_Score_z", feature_id = NA_character_, short_label = "Acid Score"),
  inflammatory_results$pairwise_summary %>%
    mutate(axis = "LM12_inflammatory_axis_z", feature_id = NA_character_, short_label = "LM12 inflammatory axis")
) %>%
  select(axis, feature_id, short_label, everything()) %>%
  arrange(axis, comparison)
write_tsv(axis_focus_summary, file.path(out_dir, "c1c4_axis_pairwise_summary_by_contrast.tsv"))

run_summary_lines <- c(
  "analysis\tC1C4_axis_validation",
  paste0("created_at\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("n_cancers\t", n_distinct(labels$Cancer)),
  paste0("n_samples\t", nrow(labels)),
  paste0("n_features\t", length(feature_ids)),
  paste0("feature_ids\t", paste(feature_ids, collapse = ";")),
  paste0("comparisons\t", paste(comparisons$comparison, collapse = ";")),
  paste0("acid_C1_vs_C4_n_group_1_gt_group_2\t", acid_results$pairwise_summary$n_group_1_gt_group_2[acid_results$pairwise_summary$comparison == "C1_vs_C4"]),
  paste0("acid_C1_vs_C4_n_group_2_gt_group_1\t", acid_results$pairwise_summary$n_group_2_gt_group_1[acid_results$pairwise_summary$comparison == "C1_vs_C4"]),
  paste0("lm12_axis_C1_vs_C2_n_group_1_gt_group_2\t", inflammatory_results$pairwise_summary$n_group_1_gt_group_2[inflammatory_results$pairwise_summary$comparison == "C1_vs_C2"]),
  paste0("lm12_axis_C3_vs_C4_n_group_1_gt_group_2\t", inflammatory_results$pairwise_summary$n_group_1_gt_group_2[inflammatory_results$pairwise_summary$comparison == "C3_vs_C4"])
)
writeLines(run_summary_lines, file.path(out_dir, "run_summary.txt"))

message("=== 02_validate_c1c4_subtype_axes.R complete ===")
message("Saved outputs to: ", out_dir)
