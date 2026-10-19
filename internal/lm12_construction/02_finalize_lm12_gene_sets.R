#!/usr/bin/env Rscript
# Purpose: Construct the adopted twelve immune pathway definitions from cross-cancer evidence.
# Task: lm12_definition; step: s02
# Inputs (DCF fields): all_votes,gene_summary,gene_long,candidate_terms,family_dictionary
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lm12_definition --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 02_finalize_lm12_gene_sets.R --config PATH [--check-config]\n",
      "       Rscript 02_finalize_lm12_gene_sets.R --help\n\n",
      "DCF fields: all_votes, gene_summary, gene_long, candidate_terms, family_dictionary, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation finalizes the study-specific expected 12 LM12 features.\n",
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
required_fields <- c("all_votes", "gene_summary", "gene_long", "candidate_terms", "family_dictionary", "output_dir")
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
all_vote_file <- unname(paths["all_votes"])
gene_summary_file <- unname(paths["gene_summary"])
gene_long_file <- unname(paths["gene_long"])
candidate_terms_file <- unname(paths["candidate_terms"])
family_dict_file <- unname(paths["family_dictionary"])
out_dir <- unname(paths["output_dir"])
required_files <- c(all_vote_file, gene_summary_file, gene_long_file, candidate_terms_file, family_dict_file)

# Protect inputs, configuration and existing results; check mode creates no directories.
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
  library(readr)
  library(stringr)
})
gene_set_dir <- file.path(out_dir, "gene_sets")
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
if (!dir.create(gene_set_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(gene_set_dir)) {
  stop("Cannot create gene_sets output directory.", call. = FALSE)
}

required_files <- c(
  all_vote_file,
  gene_summary_file,
  gene_long_file,
  candidate_terms_file,
  family_dict_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("缺失输入文件:\n", paste(missing_files, collapse = "\n"))
}

strict_recurrent_cutoff <- 0.80
strict_core_cutoff <- 0.90
original_like_core_cutoff <- 0.60
original_like_minor_cutoff <- 0.40

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_.-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

annotation_df <- data.frame(
  pathway = c(
    "Antiviral mechanism by IFN-stimulated genes",
    "MHC class II antigen presentation",
    "Fcgamma receptor (FCGR) dependent phagocytosis",
    "Interleukin-1 family signaling",
    "Class I MHC mediated antigen processing & presentation",
    "Neutrophil degranulation",
    "Interleukin-17 signaling",
    "Cytosolic sensors of pathogen-associated DNA",
    "Toll-like Receptor Cascades",
    "TCR signaling",
    "C-type lectin receptors (CLRs)",
    "Antimicrobial peptides"
  ),
  feature_id = c(
    "LM12_IFN_ISG",
    "LM12_MHC_II",
    "LM12_FCGR_Phagocytosis",
    "LM12_IL1",
    "LM12_MHC_I",
    "LM12_Neutrophil_Degranulation",
    "LM12_IL17",
    "LM12_Cytosolic_DNA_Sensing",
    "LM12_TLR_Cascades",
    "LM12_TCR_Signaling",
    "LM12_CLR",
    "LM12_Antimicrobial_Peptides"
  ),
  short_label = c(
    "IFN_ISG",
    "MHC_II",
    "FCGR_Phagocytosis",
    "IL1",
    "MHC_I",
    "Neutrophil_Degranulation",
    "IL17",
    "Cytosolic_DNA_Sensing",
    "TLR_Cascades",
    "TCR_Signaling",
    "CLR",
    "Antimicrobial_Peptides"
  ),
  functional_group = c(
    "ifn_antigen_presentation",
    "ifn_antigen_presentation",
    "myeloid_inflammation",
    "myeloid_inflammation",
    "ifn_antigen_presentation",
    "myeloid_inflammation",
    "myeloid_inflammation",
    "innate_sensing",
    "innate_sensing",
    "t_cell_signaling",
    "innate_sensing",
    "low_acid_barrier_defense"
  ),
  display_order = seq_len(12),
  stringsAsFactors = FALSE
)

expected_pathways <- annotation_df$pathway

all_votes <- read_tsv(all_vote_file, show_col_types = FALSE)
gene_summary <- read_tsv(gene_summary_file, show_col_types = FALSE) %>%
  select(direction, pathway, stable_gene_count, leading_edge_threshold_count)
gene_long <- read_tsv(gene_long_file, show_col_types = FALSE)
candidate_terms <- read_tsv(candidate_terms_file, show_col_types = FALSE)
family_dict <- read_tsv(family_dict_file, show_col_types = FALSE)

primary_candidates <- candidate_terms %>%
  filter(analysis_tier == "primary") %>%
  distinct(candidate_set, scope_layer, family_id, local_term, match_type, note)

candidate_votes <- all_votes %>%
  inner_join(primary_candidates, by = c("pathway" = "local_term")) %>%
  left_join(gene_summary, by = c("direction", "pathway")) %>%
  mutate(
    stable_gene_count = replace_na(stable_gene_count, 0L),
    original_like_class = case_when(
      vote_freq >= original_like_core_cutoff ~ "Core",
      vote_freq >= original_like_minor_cutoff ~ "Minor",
      TRUE ~ "Sporadic"
    ),
    strict_class = case_when(
      vote_freq >= strict_core_cutoff ~ "Core",
      vote_freq >= strict_recurrent_cutoff ~ "Minor",
      TRUE ~ "Below_Strict"
    )
  ) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)

write_tsv(candidate_votes, file.path(out_dir, "predefined_primary_candidate_all_votes.tsv"))

selected_original_like <- candidate_votes %>%
  filter(vote_freq >= original_like_minor_cutoff) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)
write_tsv(
  selected_original_like,
  file.path(out_dir, "selected_primary_candidates_original_like_60_40.tsv")
)

selected_strict <- candidate_votes %>%
  filter(vote_freq >= strict_recurrent_cutoff) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), desc(stable_gene_count), pathway)
write_tsv(selected_strict, file.path(out_dir, "selected_primary_candidates_strict_80.tsv"))

unexpected_pathways <- setdiff(selected_strict$pathway, expected_pathways)
missing_expected_pathways <- setdiff(expected_pathways, selected_strict$pathway)
if (length(unexpected_pathways) > 0 || length(missing_expected_pathways) > 0) {
  stop(
    paste0(
      "严格 0.80 结果与预期 12 条 pathway 不一致。\nUnexpected: ",
      paste(unexpected_pathways, collapse = "; "),
      "\nMissing: ",
      paste(missing_expected_pathways, collapse = "; ")
    )
  )
}

lm_imm_12_manifest <- selected_strict %>%
  left_join(annotation_df, by = "pathway") %>%
  left_join(
    family_dict %>%
      distinct(
        scope_layer,
        family_id,
        family_label = fig1_or_project_basis,
        family_description = note
      ),
    by = c("scope_layer", "family_id")
  ) %>%
  mutate(
    selected_rule = "strict_primary_candidate_vote_freq_ge_0.80",
    feature_class = strict_class,
    direction_label = direction
  ) %>%
  arrange(display_order) %>%
  select(
    display_order,
    feature_id,
    short_label,
    pathway,
    direction_label,
    feature_class,
    vote_count,
    vote_freq,
    mean_NES,
    median_NES,
    stable_gene_count,
    leading_edge_threshold_count,
    scope_layer,
    family_id,
    family_label,
    functional_group,
    selected_rule
  )

lm_imm_12_gene_sets_long <- gene_long %>%
  semi_join(lm_imm_12_manifest, by = c("pathway")) %>%
  left_join(
    lm_imm_12_manifest %>%
      select(feature_id, short_label, pathway, functional_group, direction_label, feature_class),
    by = "pathway"
  ) %>%
  arrange(match(feature_id, lm_imm_12_manifest$feature_id), desc(leading_edge_count), gene) %>%
  select(
    feature_id,
    short_label,
    pathway,
    functional_group,
    direction = direction_label,
    feature_class,
    gene,
    leading_edge_count,
    vote_count,
    leading_edge_threshold_count
  )

for (feature_name in lm_imm_12_manifest$feature_id) {
  sub_df <- lm_imm_12_gene_sets_long %>%
    filter(feature_id == feature_name)
  pathway_name <- unique(sub_df$pathway)
  filename <- paste0(feature_name, "__", safe_filename(pathway_name), ".txt")
  writeLines(unique(sub_df$gene), file.path(gene_set_dir, filename))
}

function_annotation_df <- lm_imm_12_manifest %>%
  mutate(
    functional_group_label = case_when(
      functional_group == "ifn_antigen_presentation" ~ "IFN_and_antigen_presentation",
      functional_group == "t_cell_signaling" ~ "T_cell_signaling",
      functional_group == "innate_sensing" ~ "Innate_sensing",
      functional_group == "myeloid_inflammation" ~ "Myeloid_inflammation",
      functional_group == "low_acid_barrier_defense" ~ "Low_acid_barrier_defense",
      TRUE ~ functional_group
    )
  ) %>%
  select(
    feature_id,
    pathway,
    functional_group,
    functional_group_label,
    family_id,
    family_label,
    direction_label
  )

write_tsv(lm_imm_12_manifest, file.path(out_dir, "lm_imm_12_manifest.tsv"))
write_tsv(lm_imm_12_manifest, file.path(out_dir, "lm_imm_12_pathway_vote_summary.tsv"))
write_tsv(lm_imm_12_gene_sets_long, file.path(out_dir, "lm_imm_12_gene_sets_long.tsv"))
write_tsv(function_annotation_df, file.path(out_dir, "lm_imm_12_functional_annotation.tsv"))

summary_lines <- c(
  "LM-Imm 12 pathway-level gene set finalization",
  paste0("Primary candidate term count: ", n_distinct(primary_candidates$local_term)),
  paste0("Strict recurrence cutoff: ", strict_recurrent_cutoff),
  paste0("Strict Core cutoff: ", strict_core_cutoff),
  paste0(
    "Strict selected features: ", nrow(lm_imm_12_manifest),
    " (High_Acid_Positive=", sum(lm_imm_12_manifest$direction_label == "High_Acid_Positive"),
    ", Low_Acid_Negative=", sum(lm_imm_12_manifest$direction_label == "Low_Acid_Negative"), ")"
  ),
  paste0(
    "Functional groups: ",
    paste(unique(function_annotation_df$functional_group_label), collapse = "; ")
  ),
  paste0("Output directory: ", out_dir)
)

writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
