#!/usr/bin/env Rscript
# Purpose: Construct the adopted twelve immune pathway definitions from cross-cancer evidence.
# Task: lm12_definition; step: s01
# Inputs (DCF fields): gsea_results_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lm12_definition --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 01_build_votes_and_stable_genes.R --config PATH [--check-config]\n",
      "       Rscript 01_build_votes_and_stable_genes.R --help\n\n",
      "DCF fields: gsea_results_dir, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation calculates recurrence votes and stable leading-edge members.\n",
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
required_fields <- c("gsea_results_dir", "output_dir")
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
full_result_dir <- unname(paths["gsea_results_dir"])
out_dir <- unname(paths["output_dir"])
# Check only file metadata using the source's filename rule; resolve symlink targets.
if (!dir.exists(full_result_dir) || file.access(full_result_dir, mode = 5L) != 0L) {
  stop("gsea_results_dir must be a readable searchable directory.", call. = FALSE)
}
gsea_metadata_files <- list.files(full_result_dir,
  pattern = "_GSEA_full_results\\.csv$", full.names = TRUE)
if (length(gsea_metadata_files) == 0L ||
    any(!file_test("-f", gsea_metadata_files)) ||
    any(file.access(gsea_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(gsea_metadata_files)$size) | file.info(gsea_metadata_files)$size <= 0)) {
  stop("Expected at least one readable nonempty regular GSEA result CSV.", call. = FALSE)
}
gsea_real_paths <- vapply(gsea_metadata_files, normalizePath, character(1L),
                         winslash = "/", mustWork = TRUE)
required_files <- gsea_real_paths


# Protect inputs, configuration and existing results; check mode creates no directories.
protected_paths <- c(full_result_dir, dirname(required_files), config_path)
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
  library(stringr)
  library(purrr)
  library(readr)
})
gene_set_dir <- file.path(out_dir, "gene_sets")
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
if (!dir.create(gene_set_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(gene_set_dir)) {
  stop("Cannot create gene_sets output directory.", call. = FALSE)
}

padj_cutoff <- 0.05
core_freq_cutoff <- 0.60
minor_freq_cutoff <- 0.40
leading_edge_freq_cutoff <- 0.70

parse_leading_edge <- function(x) {
  x <- ifelse(is.na(x), "", x)
  genes <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes[genes != ""]
}

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_.-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

gsea_files <- list.files(
  full_result_dir,
  pattern = "_GSEA_full_results\\.csv$",
  full.names = TRUE
)
if (length(gsea_files) == 0) {
  stop("未找到每癌种 GSEA 结果文件。")
}

message("读取每癌种 GSEA 结果...")
gsea_df <- map_dfr(gsea_files, function(f) {
  cancer <- sub("_GSEA_full_results\\.csv$", "", basename(f))
  read.csv(f, check.names = FALSE) %>%
    mutate(cancer = cancer)
})

n_cancers <- length(unique(gsea_df$cancer))
message("GSEA 癌种数: ", n_cancers)

select_pathways_one_direction <- function(direction_label, direction_filter) {
  sig_df <- gsea_df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    filter(direction_filter(NES)) %>%
    mutate(direction = direction_label)

  vote_df <- sig_df %>%
    distinct(direction, pathway, cancer, .keep_all = TRUE) %>%
    group_by(direction, pathway) %>%
    summarise(
      vote_count = n_distinct(cancer),
      vote_freq = vote_count / n_cancers,
      mean_NES = mean(NES, na.rm = TRUE),
      median_NES = median(NES, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      recurrence_class = case_when(
        vote_freq >= core_freq_cutoff ~ "Core",
        vote_freq >= minor_freq_cutoff ~ "Minor",
        TRUE ~ "Sporadic"
      )
    ) %>%
    arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)

  selected_df <- vote_df %>%
    filter(vote_freq >= minor_freq_cutoff)

  list(sig_df = sig_df, vote_df = vote_df, selected_df = selected_df)
}

high_res <- select_pathways_one_direction("High_Acid_Positive", function(x) x > 0)
low_res <- select_pathways_one_direction("Low_Acid_Negative", function(x) x < 0)

all_vote_df <- bind_rows(high_res$vote_df, low_res$vote_df) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)
selected_pathways <- bind_rows(high_res$selected_df, low_res$selected_df) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)

write_tsv(all_vote_df, file.path(out_dir, "all_pathway_votes_original_like.tsv"))
write_tsv(selected_pathways, file.path(out_dir, "selected_core_minor_pathways_original_like.tsv"))

message("构建 pathway-level stable leading-edge gene set...")
sig_all <- bind_rows(high_res$sig_df, low_res$sig_df)
gene_set_rows <- list()
gene_set_summary <- selected_pathways %>%
  mutate(
    stable_gene_count = 0L,
    leading_edge_threshold_count = ceiling(vote_count * leading_edge_freq_cutoff)
  )

for (i in seq_len(nrow(selected_pathways))) {
  row <- selected_pathways[i, ]
  sub_df <- sig_all %>%
    filter(direction == row$direction, pathway == row$pathway) %>%
    distinct(cancer, .keep_all = TRUE)

  leading_genes <- unlist(lapply(sub_df$leadingEdge, parse_leading_edge), use.names = FALSE)
  gene_counts <- sort(table(leading_genes), decreasing = TRUE)
  threshold_count <- ceiling(row$vote_count * leading_edge_freq_cutoff)
  stable_genes <- names(gene_counts)[as.integer(gene_counts) >= threshold_count]

  gene_set_summary$stable_gene_count[i] <- length(stable_genes)

  if (length(stable_genes) > 0) {
    gene_set_rows[[length(gene_set_rows) + 1]] <- data.frame(
      direction = row$direction,
      pathway = row$pathway,
      recurrence_class = row$recurrence_class,
      gene = stable_genes,
      leading_edge_count = as.integer(gene_counts[stable_genes]),
      vote_count = row$vote_count,
      leading_edge_threshold_count = threshold_count,
      stringsAsFactors = FALSE
    )

    prefix <- paste(
      row$direction,
      row$recurrence_class,
      safe_filename(row$pathway),
      sep = "__"
    )
    writeLines(unique(stable_genes), file.path(gene_set_dir, paste0(prefix, ".txt")))
  }
}

gene_set_long <- bind_rows(gene_set_rows) %>%
  arrange(direction, desc(vote_count), pathway, desc(leading_edge_count), gene)

write_tsv(gene_set_summary, file.path(out_dir, "original_like_gene_set_summary.tsv"))
write_tsv(gene_set_long, file.path(out_dir, "original_like_gene_sets_long.tsv"))

summary_lines <- c(
  "Original-like pathway audit",
  paste0("Cancer types: ", n_cancers),
  paste0("padj cutoff: ", padj_cutoff),
  paste0("Recurrence rule: Core >= ", core_freq_cutoff, ", Minor >= ", minor_freq_cutoff),
  paste0("Stable leading-edge rule: freq >= ", leading_edge_freq_cutoff),
  paste0(
    "Selected pathways: ", nrow(selected_pathways),
    " (High_Acid_Positive=", sum(selected_pathways$direction == "High_Acid_Positive"),
    ", Low_Acid_Negative=", sum(selected_pathways$direction == "Low_Acid_Negative"), ")"
  ),
  paste0("Non-empty stable gene sets: ", n_distinct(gene_set_long$pathway)),
  paste0("Output directory: ", out_dir)
)

writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
