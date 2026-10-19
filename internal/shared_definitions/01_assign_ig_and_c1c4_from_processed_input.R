#!/usr/bin/env Rscript
# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Task: c1c4_state_assignment; step: s01
# Inputs (DCF fields): processed_data,lm12_manifest,lm12_genes
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_state_assignment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 01_assign_ig_and_c1c4_from_processed_input.R --config PATH [--check-config] [--test] [--cancers=BRCA,COAD]\n",
      "       Rscript 01_assign_ig_and_c1c4_from_processed_input.R --help\n\n",
      "DCF fields: processed_data, lm12_manifest, lm12_genes, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "--test runs analysis for the first three input cancers unless --cancers is supplied.\n",
      "--cancers selects named cancers for actual analysis; unknown names are intersected with input names.\n",
      "A configured invocation without --check-config runs analysis and writes to an independent empty output.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (length(args) < 2L || args[1L] != "--config" || !nzchar(args[2L])) {
  usage()
  stop("Expected --config PATH followed by optional workflow flags.", call. = FALSE)
}
flags <- if (length(args) > 2L) args[-c(1L, 2L)] else character()
if (any(!(flags %in% c("--check-config", "--test") | grepl("^--cancers=.+$", flags))) ||
    sum(flags == "--check-config") > 1L || sum(flags == "--test") > 1L) {
  stop("Unknown or repeated flag. Use --help.", call. = FALSE)
}
is_test <- "--test" %in% args
cancers_arg <- grep("^--cancers=", args, value = TRUE)

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
required_fields <- c("processed_data", "lm12_manifest", "lm12_genes", "output_dir")
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
input_data_file <- unname(paths["processed_data"])
manifest_file <- unname(paths["lm12_manifest"])
gene_long_file <- unname(paths["lm12_genes"])
result_dir <- unname(paths["output_dir"])
required_files <- c(input_data_file, manifest_file, gene_long_file)

# 防止输出覆盖输入、配置或既有结果；检查模式不会创建目录。
protected_paths <- c(dirname(required_files), config_path)
if (result_dir == "/" || result_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, result_dir) || is_within(result_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(result_dir) && !dir.exists(result_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(result_dir)) {
  if (file.access(result_dir, mode = 7L) != 0L ||
      length(list.files(result_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All three inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- result_dir
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
  library(GSVA)
})
if (!dir.create(result_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(result_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

message("=== 6.0 AG3 x IG2 -> C1-C4 subtype build ===")
message("Input data: ", input_data_file)
message("LM12 manifest: ", manifest_file)
message("Output dir: ", result_dir)

load(input_data_file)
if (!exists("pancancer_data_list")) {
  stop("Object `pancancer_data_list` not found in ", input_data_file)
}

manifest <- read_tsv(manifest_file, show_col_types = FALSE) %>%
  arrange(display_order)
gene_long <- read_tsv(gene_long_file, show_col_types = FALSE)

feature_metadata <- manifest %>%
  transmute(
    display_order,
    feature_id,
    short_label,
    pathway,
    functional_group,
    direction_label,
    stable_gene_count
  )

feature_ids <- feature_metadata$feature_id
immune_high_features <- setdiff(feature_ids, "LM12_Antimicrobial_Peptides")
row_label_map <- setNames(feature_metadata$short_label, feature_metadata$feature_id)

gene_sets <- lapply(feature_ids, function(fid) {
  gene_long %>%
    filter(feature_id == fid) %>%
    pull(gene) %>%
    unique() %>%
    sort()
})
names(gene_sets) <- feature_ids

selected_cancers <- names(pancancer_data_list)
if (length(cancers_arg) > 0) {
  selected_cancers <- strsplit(sub("^--cancers=", "", cancers_arg[1]), ",")[[1]]
}
if (is_test && length(cancers_arg) == 0) {
  selected_cancers <- head(selected_cancers, 3)
}
selected_cancers <- intersect(selected_cancers, names(pancancer_data_list))
if (length(selected_cancers) == 0) {
  stop("没有可处理的癌种。")
}

message("Cancers to process: ", paste(selected_cancers, collapse = ", "))

ag_levels <- c("AG1", "AG2", "AG3")
ig_levels <- c("IG1", "IG2")
c_levels <- c("C1", "C2", "C3", "C4")
agig_levels <- as.vector(outer(ag_levels, ig_levels, paste, sep = "_"))

ag_map <- c(High = "AG1", Medium = "AG2", Low = "AG3")
ag_label_map <- c(
  AG1 = "acid_high",
  AG2 = "acid_intermediate",
  AG3 = "acid_low"
)
ig_label_map <- c(
  IG1 = "lm12_inflammatory_remodeling_high",
  IG2 = "lm12_inflammatory_remodeling_low"
)
c_map <- c(
  AG1_IG1 = "C1",
  AG1_IG2 = "C2",
  AG2_IG1 = "C3",
  AG2_IG2 = "C2",
  AG3_IG1 = "C3",
  AG3_IG2 = "C4"
)
c_label_map <- c(
  C1 = "high_acid__high_inflammatory_remodeling",
  C2 = "high_or_intermediate_acid__low_inflammatory_remodeling",
  C3 = "intermediate_or_low_acid__high_inflammatory_remodeling",
  C4 = "low_acid__low_inflammatory_remodeling"
)

score_list <- list()
coverage_list <- list()
ig_label_list <- list()
ig_manifest_list <- list()
ig_center_list <- list()
sample_count_records <- list()

for (cancer in selected_cancers) {
  message("\n--- Processing ", cancer, " ---")
  current_data <- pancancer_data_list[[cancer]]
  expr_matrix <- current_data$expression_data
  acid_scores <- current_data$acid_scores
  clinical_data <- current_data$clinical_data

  if (is.null(expr_matrix) || is.null(acid_scores) || is.null(clinical_data)) {
    stop("癌种 ", cancer, " 缺失 expression_data、acid_scores 或 clinical_data。")
  }
  if (!"sample" %in% colnames(clinical_data)) {
    stop("癌种 ", cancer, " 的 clinical_data 缺失 sample 列。")
  }
  if (!"Group" %in% colnames(clinical_data)) {
    stop("癌种 ", cancer, " 的 clinical_data 缺失 Group 列，无法定义 AG。")
  }

  common_samples <- Reduce(
    intersect,
    list(colnames(expr_matrix), colnames(acid_scores), clinical_data$sample)
  )
  if (length(common_samples) < 10) {
    stop("癌种 ", cancer, " 共同样本数 < 10。")
  }

  clinical_use <- clinical_data %>%
    filter(sample %in% common_samples) %>%
    distinct(sample, .keep_all = TRUE)
  ordered_samples <- clinical_use$sample

  expr_use <- expr_matrix[, ordered_samples, drop = FALSE]
  acid_values <- as.numeric(acid_scores[1, ordered_samples])
  names(acid_values) <- ordered_samples

  if (any(is.na(acid_values))) {
    stop("癌种 ", cancer, " 的 Acid_Score 存在 NA。")
  }

  expr_genes <- rownames(expr_use)
  coverage_df <- bind_rows(lapply(feature_ids, function(feature_id) {
    genes_total <- gene_sets[[feature_id]]
    genes_present <- intersect(genes_total, expr_genes)
    tibble(
      Cancer = cancer,
      feature_id = feature_id,
      short_label = row_label_map[[feature_id]],
      genes_total = length(genes_total),
      genes_present = length(genes_present),
      coverage = length(genes_present) / length(genes_total)
    )
  }))
  coverage_list[[cancer]] <- coverage_df

  gene_sets_present <- lapply(feature_ids, function(feature_id) {
    intersect(gene_sets[[feature_id]], expr_genes)
  })
  names(gene_sets_present) <- feature_ids

  zero_coverage_features <- names(gene_sets_present)[lengths(gene_sets_present) == 0]
  if (length(zero_coverage_features) > 0) {
    stop(
      "癌种 ", cancer, " 存在 gene set 完全无覆盖: ",
      paste(zero_coverage_features, collapse = ", ")
    )
  }

  expr_log2 <- log2(expr_use + 1)
  ssgsea_param <- ssgseaParam(
    exprData = expr_log2,
    geneSets = gene_sets_present,
    normalize = TRUE
  )
  score_matrix <- gsva(ssgsea_param, verbose = FALSE)

  score_df <- as.data.frame(t(score_matrix)) %>%
    rownames_to_column("Sample") %>%
    mutate(
      Cancer = cancer,
      Acid_Score = acid_values[Sample],
      Acid_Group = clinical_use$Group[match(Sample, clinical_use$sample)],
      .after = Sample
    )

  unknown_groups <- setdiff(unique(score_df$Acid_Group), names(ag_map))
  if (length(unknown_groups) > 0) {
    stop("癌种 ", cancer, " 出现未知 Acid_Group: ", paste(unknown_groups, collapse = ", "))
  }

  score_list[[cancer]] <- score_df
  sample_count_records[[cancer]] <- tibble(Cancer = cancer, n_samples = nrow(score_df))

  immune_axis <- rowMeans(score_df[, immune_high_features, drop = FALSE], na.rm = TRUE)
  cluster_input <- as.matrix(score_df[, feature_ids, drop = FALSE])
  rownames(cluster_input) <- score_df$Sample

  dist_matrix <- dist(cluster_input, method = "euclidean")
  hclust_res <- hclust(dist_matrix, method = "ward.D")
  raw_labels <- cutree(hclust_res, k = 2)

  tmp_label <- tibble(
    Sample = names(raw_labels),
    Cancer = cancer,
    IG2_raw_cluster = as.integer(raw_labels),
    Acid_Score = score_df$Acid_Score[match(names(raw_labels), score_df$Sample)],
    immune_inflammatory_axis = immune_axis[match(names(raw_labels), score_df$Sample)]
  )

  cluster_stats <- tmp_label %>%
    group_by(IG2_raw_cluster) %>%
    summarise(
      n_samples = n(),
      mean_immune_inflammatory_axis = mean(immune_inflammatory_axis, na.rm = TRUE),
      median_immune_inflammatory_axis = median(immune_inflammatory_axis, na.rm = TRUE),
      mean_acid_score = mean(Acid_Score, na.rm = TRUE),
      median_acid_score = median(Acid_Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_immune_inflammatory_axis)) %>%
    mutate(
      IG = factor(c("IG1", "IG2"), levels = ig_levels),
      IG_Label = unname(ig_label_map[as.character(IG)]),
      immune_rank_within_cancer = c(1L, 2L),
      Cancer = cancer,
      .before = 1
    )

  ig_label_list[[cancer]] <- tmp_label %>%
    select(Sample, Cancer, IG2_raw_cluster) %>%
    left_join(
      cluster_stats %>%
        select(Cancer, IG2_raw_cluster, IG, IG_Label, immune_rank_within_cancer),
      by = c("Cancer", "IG2_raw_cluster")
    )

  ig_manifest_list[[cancer]] <- cluster_stats %>%
    select(Cancer, IG, IG_Label, IG2_raw_cluster, immune_rank_within_cancer, everything())

  ig_center_list[[cancer]] <- score_df %>%
    select(Sample, all_of(feature_ids)) %>%
    inner_join(ig_label_list[[cancer]] %>% select(Sample, IG), by = "Sample") %>%
    group_by(IG) %>%
    summarise(across(all_of(feature_ids), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(Cancer = cancer, .before = 1)

  message("  samples=", nrow(score_df), "; IG cluster sizes=", paste(cluster_stats$n_samples, collapse = "/"))
}

score_df_all <- bind_rows(score_list) %>%
  arrange(Cancer, Sample)
coverage_summary <- bind_rows(coverage_list) %>%
  left_join(feature_metadata, by = "feature_id") %>%
  select(
    Cancer,
    display_order,
    feature_id,
    short_label = short_label.x,
    pathway,
    functional_group,
    genes_total,
    genes_present,
    coverage
  ) %>%
  arrange(Cancer, display_order)

ig_labels <- bind_rows(ig_label_list)
ig_manifest <- bind_rows(ig_manifest_list) %>%
  arrange(Cancer, immune_rank_within_cancer)
ig_centers <- bind_rows(ig_center_list) %>%
  arrange(Cancer, IG)
sample_count_summary <- bind_rows(sample_count_records) %>%
  arrange(Cancer)

label_df <- score_df_all %>%
  left_join(ig_labels, by = c("Sample", "Cancer")) %>%
  mutate(
    AG = unname(ag_map[Acid_Group]),
    AG_Label = unname(ag_label_map[AG]),
    AGIG = paste(AG, as.character(IG), sep = "_"),
    C_subtype = unname(c_map[AGIG]),
    C_label = unname(c_label_map[C_subtype]),
    immune_inflammatory_axis = rowMeans(across(all_of(immune_high_features)), na.rm = TRUE)
  )

if (any(is.na(label_df$AG))) {
  stop("AG 映射失败，存在未分配样本。")
}
if (any(is.na(label_df$IG))) {
  stop("IG 映射失败，存在未分配样本。")
}
if (any(is.na(label_df$C_subtype))) {
  stop("C_subtype 映射失败，存在未分配 AGIG。")
}

label_df <- label_df %>%
  mutate(
    AG = factor(AG, levels = ag_levels),
    IG = factor(IG, levels = ig_levels),
    AGIG = factor(AGIG, levels = agig_levels),
    C_subtype = factor(C_subtype, levels = c_levels)
  ) %>%
  arrange(Cancer, C_subtype, AG, IG, Sample)

atomic_label_df <- label_df %>%
  select(
    Sample,
    Cancer,
    Acid_Score,
    Acid_Group,
    AG,
    AG_Label,
    IG2_raw_cluster,
    IG,
    IG_Label,
    immune_rank_within_cancer,
    immune_inflammatory_axis,
    AGIG,
    C_subtype,
    C_label,
    all_of(feature_ids)
  )

counts_overall <- atomic_label_df %>%
  count(C_subtype, AGIG, AG, IG, C_label, name = "n_samples") %>%
  arrange(C_subtype, AG, IG)

counts_by_cancer_long <- atomic_label_df %>%
  count(Cancer, C_subtype, name = "n_samples") %>%
  complete(
    Cancer = sort(unique(atomic_label_df$Cancer)),
    C_subtype = factor(c_levels, levels = c_levels),
    fill = list(n_samples = 0L)
  ) %>%
  left_join(
    tibble(C_subtype = factor(c_levels, levels = c_levels), C_label = unname(c_label_map[c_levels])),
    by = "C_subtype"
  ) %>%
  arrange(Cancer, C_subtype)

counts_by_cancer_wide <- counts_by_cancer_long %>%
  select(Cancer, C_subtype, n_samples) %>%
  pivot_wider(names_from = C_subtype, values_from = n_samples, values_fill = 0) %>%
  mutate(
    total_samples = rowSums(across(all_of(c_levels))),
    all_groups_present = if_all(all_of(c_levels), ~ .x > 0),
    all_groups_n_ge_10 = if_all(all_of(c_levels), ~ .x >= 10),
    all_groups_n_ge_20 = if_all(all_of(c_levels), ~ .x >= 20),
    min_group_n = pmin(C1, C2, C3, C4)
  ) %>%
  arrange(Cancer)

agig_counts_by_cancer <- atomic_label_df %>%
  count(Cancer, AGIG, AG, IG, name = "n_samples") %>%
  complete(
    Cancer = sort(unique(atomic_label_df$Cancer)),
    AGIG = factor(agig_levels, levels = agig_levels),
    fill = list(n_samples = 0L)
  ) %>%
  mutate(
    AG = factor(sub("_IG[12]$", "", as.character(AGIG)), levels = ag_levels),
    IG = factor(sub("^AG[123]_", "", as.character(AGIG)), levels = ig_levels)
  ) %>%
  arrange(Cancer, AG, IG)

write_tsv(feature_metadata, file.path(result_dir, "lm12_feature_metadata.tsv"))
write_tsv(coverage_summary, file.path(result_dir, "lm12_gene_coverage.tsv"))
write_tsv(score_df_all, file.path(result_dir, "pancancer_lm12_scores_by_sample.tsv.gz"))
write_tsv(ig_manifest, file.path(result_dir, "ig2_cluster_manifest_by_cancer.tsv"))
write_tsv(ig_centers, file.path(result_dir, "ig2_cluster_feature_medians_by_cancer.tsv"))
write_csv(atomic_label_df, file.path(result_dir, "ag3ig2_c1c4_sample_labels.csv"))
write_tsv(counts_overall, file.path(result_dir, "ag3ig2_c1c4_counts_overall.tsv"))
write_tsv(counts_by_cancer_long, file.path(result_dir, "ag3ig2_c1c4_counts_by_cancer_long.tsv"))
write_tsv(counts_by_cancer_wide, file.path(result_dir, "ag3ig2_c1c4_counts_by_cancer_wide.tsv"))
write_tsv(agig_counts_by_cancer, file.path(result_dir, "ag3ig2_atomic_counts_by_cancer.tsv"))
write_tsv(sample_count_summary, file.path(result_dir, "sample_count_summary_by_cancer.tsv"))

save(
  score_list,
  score_df_all,
  feature_metadata,
  gene_sets,
  coverage_summary,
  ig_labels,
  ig_manifest,
  ig_centers,
  atomic_label_df,
  counts_overall,
  counts_by_cancer_long,
  counts_by_cancer_wide,
  agig_counts_by_cancer,
  file = file.path(result_dir, "ag3ig2_c1c4_subtype_build_objects.Rdata")
)

summary_lines <- c(
  paste("analysis", "AG3xIG2_C1C4_formal_subtype_build", sep = "\t"),
  paste("created_at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sep = "\t"),
  paste("test_mode", is_test, sep = "\t"),
  paste("n_cancers", length(selected_cancers), sep = "\t"),
  paste("n_samples", nrow(atomic_label_df), sep = "\t"),
  paste("n_features", length(feature_ids), sep = "\t"),
  paste("feature_ids", paste(feature_ids, collapse = ";"), sep = "\t"),
  paste("n_cancers_all_groups_present", sum(counts_by_cancer_wide$all_groups_present), sep = "\t"),
  paste("n_cancers_all_groups_n_ge_10", sum(counts_by_cancer_wide$all_groups_n_ge_10), sep = "\t"),
  paste("n_cancers_all_groups_n_ge_20", sum(counts_by_cancer_wide$all_groups_n_ge_20), sep = "\t"),
  paste("counts_overall", paste(paste(counts_overall$C_subtype, counts_overall$AGIG, counts_overall$n_samples, sep = ":"), collapse = ";"), sep = "\t")
)
writeLines(summary_lines, file.path(result_dir, "run_summary.txt"))

message("\nSaved outputs to: ", result_dir)
message("=== 01_build_ag3ig2_c1c4_subtypes.R complete ===")
