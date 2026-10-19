#!/usr/bin/env Rscript
# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Task: c1c4_state_assignment; step: s03
# Inputs (DCF fields): assignment_labels,axis_scores,feature_metadata
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_state_assignment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 03_export_canonical_labels.R --config PATH [--check-config]\n",
      "       Rscript 03_export_canonical_labels.R --help\n\n",
      "DCF fields: assignment_labels, axis_scores, feature_metadata, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation exports an independent copy in the canonical label format.\n",
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
required_fields <- c("assignment_labels", "axis_scores", "feature_metadata", "output_dir")
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
label_01_file <- unname(paths["assignment_labels"])
axis_score_file <- unname(paths["axis_scores"])
feature_metadata_file <- unname(paths["feature_metadata"])
out_dir <- unname(paths["output_dir"])
required_files <- c(label_01_file, axis_score_file, feature_metadata_file)

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

canonical_version <- "AG3IG2_C1C4_v1_20260615"
created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

required_files <- c(label_01_file, axis_score_file, feature_metadata_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("缺失输入文件:\n", paste(missing_files, collapse = "\n"))
}

c_levels <- c("C1", "C2", "C3", "C4")
ag_levels <- c("AG1", "AG2", "AG3")
ig_levels <- c("IG1", "IG2")
c_map <- c(
  AG1_IG1 = "C1",
  AG1_IG2 = "C2",
  AG2_IG2 = "C2",
  AG2_IG1 = "C3",
  AG3_IG1 = "C3",
  AG3_IG2 = "C4"
)

message("=== 6.0 freeze canonical AG3 x IG2 -> C1-C4 labels ===")
message("Input 01 labels: ", label_01_file)
message("Input 02 axis scores: ", axis_score_file)
message("Output dir: ", out_dir)

feature_metadata <- read_tsv(feature_metadata_file, show_col_types = FALSE) %>%
  arrange(display_order)
feature_ids <- feature_metadata$feature_id
feature_z_cols <- paste0(feature_ids, "_z")

labels_01 <- read_csv(label_01_file, show_col_types = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    Cancer = as.character(Cancer),
    AG = as.character(AG),
    IG = as.character(IG),
    AGIG = as.character(AGIG),
    C_subtype = as.character(C_subtype)
  )

axis_scores <- read_tsv(axis_score_file, show_col_types = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    Cancer = as.character(Cancer),
    AG = as.character(AG),
    IG = as.character(IG),
    AGIG = as.character(AGIG),
    C_subtype = as.character(C_subtype)
  )

required_01_cols <- c(
  "Sample", "Cancer", "Acid_Group", "AG", "AG_Label", "IG2_raw_cluster",
  "IG", "IG_Label", "immune_rank_within_cancer", "AGIG", "C_subtype",
  "C_label", feature_ids
)
required_axis_cols <- c(
  "Sample", "Cancer", "AG", "IG", "AGIG", "C_subtype", "C_label",
  "Acid_Score", "Acid_Score_z", "LM12_inflammatory_axis",
  "LM12_inflammatory_axis_z", "LM12_all_feature_axis",
  "LM12_all_feature_axis_z", feature_ids, feature_z_cols
)

missing_01_cols <- setdiff(required_01_cols, names(labels_01))
missing_axis_cols <- setdiff(required_axis_cols, names(axis_scores))
if (length(missing_01_cols) > 0) {
  stop("01 标签表缺失必要列: ", paste(missing_01_cols, collapse = ", "))
}
if (length(missing_axis_cols) > 0) {
  stop("02 轴得分表缺失必要列: ", paste(missing_axis_cols, collapse = ", "))
}

message("[1/5] 输入一致性校验")
dup_01 <- labels_01 %>% count(Sample, Cancer) %>% filter(n > 1)
dup_axis <- axis_scores %>% count(Sample, Cancer) %>% filter(n > 1)
if (nrow(dup_01) > 0 || nrow(dup_axis) > 0) {
  stop("存在重复 Sample-Cancer 键，无法冻结 canonical label。")
}

compare_cols <- c("AG", "IG", "AGIG", "C_subtype")
joined_check <- labels_01 %>%
  select(Sample, Cancer, all_of(compare_cols)) %>%
  inner_join(
    axis_scores %>% select(Sample, Cancer, all_of(compare_cols)),
    by = c("Sample", "Cancer"),
    suffix = c("_01", "_axis")
  )

if (nrow(joined_check) != nrow(labels_01) || nrow(joined_check) != nrow(axis_scores)) {
  stop("01 标签表和 02 轴得分表的 Sample-Cancer 键不完全一致。")
}

for (col in compare_cols) {
  mismatch_n <- sum(joined_check[[paste0(col, "_01")]] != joined_check[[paste0(col, "_axis")]], na.rm = TRUE)
  if (mismatch_n > 0) {
    stop("01 与 02 的 ", col, " 不一致，mismatch_n = ", mismatch_n)
  }
}

mapped_c <- unname(c_map[labels_01$AGIG])
if (any(is.na(mapped_c)) || any(mapped_c != labels_01$C_subtype)) {
  bad <- labels_01 %>% filter(is.na(mapped_c) | mapped_c != C_subtype) %>% distinct(AGIG, C_subtype)
  stop("AGIG -> C_subtype 确定性映射校验失败:\n", paste(capture.output(print(bad)), collapse = "\n"))
}

message("[2/5] 构建 canonical full/minimal label 表")
metadata_cols <- c(
  "Sample", "Cancer", "Acid_Group", "AG_Label", "IG2_raw_cluster",
  "IG_Label", "immune_rank_within_cancer"
)
labels_meta <- labels_01 %>%
  select(all_of(metadata_cols))

canonical_labels <- axis_scores %>%
  left_join(labels_meta, by = c("Sample", "Cancer")) %>%
  mutate(
    canonical_version = canonical_version,
    C_subtype = factor(C_subtype, levels = c_levels),
    AG = factor(AG, levels = ag_levels),
    IG = factor(IG, levels = ig_levels),
    AGIG = factor(AGIG, levels = names(c_map)),
    freeze_date = created_at
  ) %>%
  arrange(Cancer, C_subtype, AG, IG, Sample) %>%
  mutate(canonical_row_id = row_number()) %>%
  select(
    canonical_version,
    canonical_row_id,
    freeze_date,
    Sample,
    Cancer,
    C_subtype,
    C_label,
    AG,
    AG_Label,
    Acid_Group,
    IG,
    IG_Label,
    AGIG,
    IG2_raw_cluster,
    immune_rank_within_cancer,
    Acid_Score,
    Acid_Score_z,
    LM12_inflammatory_axis,
    LM12_inflammatory_axis_z,
    LM12_all_feature_axis,
    LM12_all_feature_axis_z,
    all_of(feature_ids),
    all_of(feature_z_cols)
  )

canonical_minimal <- canonical_labels %>%
  select(
    canonical_version,
    canonical_row_id,
    freeze_date,
    Sample,
    Cancer,
    C_subtype,
    C_label,
    AG,
    AG_Label,
    Acid_Group,
    IG,
    IG_Label,
    AGIG,
    Acid_Score,
    Acid_Score_z,
    LM12_inflammatory_axis,
    LM12_inflammatory_axis_z
  )

message("[3/5] 生成样本量和字典")
counts_long <- canonical_labels %>%
  count(Cancer, C_subtype, C_label, name = "n_samples") %>%
  complete(
    Cancer,
    C_subtype = factor(c_levels, levels = c_levels),
    fill = list(n_samples = 0)
  ) %>%
  group_by(Cancer) %>%
  mutate(
    total_samples = sum(n_samples),
    all_groups_present = all(n_samples > 0),
    all_groups_n_ge_10 = all(n_samples >= 10),
    all_groups_n_ge_20 = all(n_samples >= 20),
    min_group_n = min(n_samples)
  ) %>%
  ungroup() %>%
  arrange(Cancer, C_subtype)

label_lookup <- canonical_labels %>%
  distinct(C_subtype, C_label)
counts_long <- counts_long %>%
  select(-C_label) %>%
  left_join(label_lookup, by = "C_subtype") %>%
  select(Cancer, C_subtype, C_label, everything())

counts_wide <- counts_long %>%
  select(Cancer, C_subtype, n_samples) %>%
  pivot_wider(names_from = C_subtype, values_from = n_samples) %>%
  mutate(
    total_samples = rowSums(across(all_of(c_levels))),
    min_group_n = do.call(pmin, c(across(all_of(c_levels)), list(na.rm = TRUE))),
    all_groups_present = if_all(all_of(c_levels), ~ .x > 0),
    all_groups_n_ge_10 = if_all(all_of(c_levels), ~ .x >= 10),
    all_groups_n_ge_20 = if_all(all_of(c_levels), ~ .x >= 20)
  ) %>%
  arrange(Cancer)

counts_overall <- canonical_labels %>%
  count(C_subtype, C_label, name = "n_samples") %>%
  complete(C_subtype = factor(c_levels, levels = c_levels), fill = list(n_samples = 0)) %>%
  select(C_subtype, n_samples) %>%
  left_join(label_lookup, by = "C_subtype") %>%
  select(C_subtype, C_label, n_samples)

subtype_dictionary <- tribble(
  ~field, ~level, ~label, ~definition, ~construction_rule, ~interpretation_note,
  "AG", "AG1", "acid_high", "Within-cancer high-acid group", "clinical_data$Group == High", "The acid axis is fixed as a within-cancer three-level stratification.",
  "AG", "AG2", "acid_intermediate", "Within-cancer intermediate-acid group", "clinical_data$Group == Medium", "Intermediate acid state is retained instead of forcing a binary contrast.",
  "AG", "AG3", "acid_low", "Within-cancer low-acid group", "clinical_data$Group == Low", "Low-acid endpoint.",
  "IG", "IG1", "lm12_inflammatory_remodeling_high", "Higher LM12 inflammatory-remodeling program", "Within-cancer hierarchical clustering k=2 on 12 LM-Imm scores; clusters are oriented by the LM12 inflammatory axis mean.", "IG is not interpreted as immune infiltration abundance or anti-tumor immune strength.",
  "IG", "IG2", "lm12_inflammatory_remodeling_low", "Lower LM12 inflammatory-remodeling program", "Within-cancer hierarchical clustering k=2 on 12 LM-Imm scores; clusters are oriented by the LM12 inflammatory axis mean.", "IG is not interpreted as immune infiltration abundance or anti-tumor immune strength.",
  "C_subtype", "C1", "high_acid__high_inflammatory_remodeling", "High acid plus high inflammatory-remodeling program", "C1 = AG1_IG1", "Primary endpoint contrast with C4.",
  "C_subtype", "C2", "high_or_intermediate_acid__low_inflammatory_remodeling", "High/intermediate acid plus low inflammatory-remodeling program", "C2 = AG1_IG2 + AG2_IG2", "Useful for disentangling acid and inflammatory-remodeling axes against C1 or C3.",
  "C_subtype", "C3", "intermediate_or_low_acid__high_inflammatory_remodeling", "Intermediate/low acid plus high inflammatory-remodeling program", "C3 = AG2_IG1 + AG3_IG1", "Useful for disentangling acid and inflammatory-remodeling axes against C2 or C4.",
  "C_subtype", "C4", "low_acid__low_inflammatory_remodeling", "Low acid plus low inflammatory-remodeling program", "C4 = AG3_IG2", "Primary endpoint contrast with C1."
)

message("[4/5] 写出冻结文件")
full_csv <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_labels.csv")
minimal_csv <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_labels_minimal.csv")
rdata_file <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_labels.Rdata")
dictionary_file <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_dictionary.tsv")
counts_long_file <- file.path(out_dir, "final_ag3ig2_c1c4_counts_by_cancer_long.tsv")
counts_wide_file <- file.path(out_dir, "final_ag3ig2_c1c4_counts_by_cancer_wide.tsv")
counts_overall_file <- file.path(out_dir, "final_ag3ig2_c1c4_counts_overall.tsv")
manifest_file <- file.path(out_dir, "final_ag3ig2_c1c4_freeze_manifest.tsv")
checksums_file <- file.path(out_dir, "final_ag3ig2_c1c4_file_checksums.tsv")

write_csv(canonical_labels, full_csv)
write_csv(canonical_minimal, minimal_csv)
write_tsv(subtype_dictionary, dictionary_file)
write_tsv(counts_long, counts_long_file)
write_tsv(counts_wide, counts_wide_file)
write_tsv(counts_overall, counts_overall_file)

freeze_metadata <- list(
  canonical_version = canonical_version,
  created_at = created_at,
  n_samples = nrow(canonical_labels),
  n_cancers = n_distinct(canonical_labels$Cancer),
  c_levels = c_levels,
  c_map = c_map,
  input_label_01_file = label_01_file,
  input_axis_score_file = axis_score_file,
  feature_metadata_file = feature_metadata_file,
  no_write_to_data_dir = TRUE
)

subtype_labels <- canonical_labels
subtype_labels_minimal <- canonical_minimal
save(
  subtype_labels,
  subtype_labels_minimal,
  subtype_dictionary,
  counts_long,
  counts_wide,
  counts_overall,
  freeze_metadata,
  file = rdata_file
)

manifest <- tibble(
  key = c(
    "canonical_version",
    "created_at",
    "n_samples",
    "n_cancers",
    "n_features",
    "n_cancers_all_groups_present",
    "n_cancers_all_groups_n_ge_10",
    "n_cancers_all_groups_n_ge_20",
    "input_label_01_file",
    "input_label_01_md5",
    "input_axis_score_file",
    "input_axis_score_md5",
    "feature_metadata_file",
    "feature_metadata_md5",
    "write_to_data_dir",
    "legacy_data_label_note"
  ),
  value = c(
    canonical_version,
    created_at,
    as.character(nrow(canonical_labels)),
    as.character(n_distinct(canonical_labels$Cancer)),
    as.character(length(feature_ids)),
    as.character(sum(counts_wide$all_groups_present)),
    as.character(sum(counts_wide$all_groups_n_ge_10)),
    as.character(sum(counts_wide$all_groups_n_ge_20)),
    label_01_file,
    unname(tools::md5sum(label_01_file)),
    axis_score_file,
    unname(tools::md5sum(axis_score_file)),
    feature_metadata_file,
    unname(tools::md5sum(feature_metadata_file)),
    "FALSE",
    "data/final_pancancer_subtype_labels.csv is legacy LM-Imm C1/C2/C3 and is not used by the current C1-C4 mainline."
  )
)
write_tsv(manifest, manifest_file)

message("[5/5] 输出校验和运行摘要")
output_files <- c(
  full_csv,
  minimal_csv,
  rdata_file,
  dictionary_file,
  counts_long_file,
  counts_wide_file,
  counts_overall_file,
  manifest_file
)
checksums <- tibble(
  file = output_files,
  md5 = unname(tools::md5sum(output_files)),
  size_bytes = file.info(output_files)$size
)
write_tsv(checksums, checksums_file)

summary_lines <- c(
  "analysis\tFreeze_AG3IG2_C1C4_canonical_labels",
  paste0("created_at\t", created_at),
  paste0("canonical_version\t", canonical_version),
  paste0("n_samples\t", nrow(canonical_labels)),
  paste0("n_cancers\t", n_distinct(canonical_labels$Cancer)),
  paste0("n_features\t", length(feature_ids)),
  paste0("counts_overall\t", paste(paste0(counts_overall$C_subtype, ":", counts_overall$n_samples), collapse = ";")),
  paste0("n_cancers_all_groups_present\t", sum(counts_wide$all_groups_present)),
  paste0("n_cancers_all_groups_n_ge_10\t", sum(counts_wide$all_groups_n_ge_10)),
  paste0("n_cancers_all_groups_n_ge_20\t", sum(counts_wide$all_groups_n_ge_20)),
  paste0("full_csv\t", full_csv),
  paste0("minimal_csv\t", minimal_csv),
  paste0("rdata_file\t", rdata_file),
  "write_to_data_dir\tFALSE"
)
writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))

message("=== 04_freeze_c1c4_labels.R complete ===")
message("Saved outputs to: ", out_dir)
