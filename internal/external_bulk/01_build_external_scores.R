#!/usr/bin/env Rscript
# Purpose: Score external bulk cohorts and test the fixed LM12 pathway enrichment.
# Task: external_bulk_enrichment; step: s01
# Inputs (DCF fields): expression_dir,metadata_dir,final_genes,lm12_manifest,lm12_genes
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task external_bulk_enrichment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 01_build_external_scores.R --config PATH [--check-config]\n",
      "       Rscript 01_build_external_scores.R --help\n\n",
      "DCF fields: expression_dir, metadata_dir, final_genes, lm12_manifest, lm12_genes, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation prepares expression and calculates Acid plus LM12 ssGSEA scores.\n",
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
required_fields <- c("expression_dir", "metadata_dir", "final_genes", "lm12_manifest", "lm12_genes", "output_dir")
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
geo_source_dir <- unname(paths["expression_dir"])
metadata_dir <- unname(paths["metadata_dir"])
acid_gene_file <- unname(paths["final_genes"])
lm_imm_12_manifest_file <- unname(paths["lm12_manifest"])
lm_imm_12_gene_long_file <- unname(paths["lm12_genes"])
out_dir <- unname(paths["output_dir"])
manifest_file <- file.path(metadata_dir, "external_cohort_manifest.tsv")
gse40419_annot_file <- file.path(metadata_dir, "GSE40419_sample_annotation.tsv")
gse127165_annot_file <- file.path(metadata_dir, "GSE127165_sample_annotation.tsv")
gse96058_annot_file <- file.path(metadata_dir, "GSE96058_expression_metadata_manifest.tsv")
expr_file_map <- c(
  GSE40419 = file.path(geo_source_dir, "GSE40419_LC-87_RPKM_expression.txt.gz"),
  GSE127165 = file.path(geo_source_dir, "GSE127165_114_samples_RNA_seq_FPKM.txt.gz"),
  GSE96058 = file.path(geo_source_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
)
gse40419_file <- unname(expr_file_map["GSE40419"])
gse127165_file <- unname(expr_file_map["GSE127165"])
gse96058_file <- unname(expr_file_map["GSE96058"])
required_files <- c(acid_gene_file, manifest_file, gse40419_annot_file,
  gse127165_annot_file, gse96058_annot_file, lm_imm_12_manifest_file,
  lm_imm_12_gene_long_file, expr_file_map)

# Inspect required fixed filenames without reading or decompressing their contents.
input_dirs <- c(geo_source_dir, metadata_dir)
if (any(!dir.exists(input_dirs)) || any(file.access(input_dirs, mode = 5L) != 0L)) {
  stop("Input directories must be readable and searchable.", call. = FALSE)
}
if (any(!file_test("-f", required_files)) || any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All fixed-name inputs must be readable nonempty regular files.", call. = FALSE)
}
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)


# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(input_dirs, dirname(required_files), config_path)
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

# Isolate the named signature object so other serialized objects cannot overwrite paths.
load_named_object <- function(path, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(path, envir = input_env)
  if (!exists(object_name, envir = input_env, inherits = FALSE)) {
    stop("Missing input object: ", object_name, call. = FALSE)
  }
  get(object_name, envir = input_env, inherits = FALSE)
}
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GSVA)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readr)
})
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

gene_alias_map <- c(
  "PKM2" = "PKM"
)

normalize_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

apply_gene_alias_map <- function(gene_symbols) {
  alias_hit <- gene_symbols %in% names(gene_alias_map)
  gene_symbols[alias_hit] <- unname(gene_alias_map[gene_symbols[alias_hit]])
  list(
    gene_symbols = gene_symbols,
    alias_applied_n = sum(alias_hit, na.rm = TRUE)
  )
}

aggregate_duplicate_rows_mean <- function(expr_mat) {
  duplicate_n <- sum(duplicated(rownames(expr_mat)))
  if (duplicate_n == 0) {
    return(list(expr_mat = expr_mat, duplicate_n = 0L))
  }
  summed <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(
    matrix(1, nrow = nrow(expr_mat), ncol = ncol(expr_mat)),
    group = rownames(expr_mat),
    reorder = FALSE
  )
  list(expr_mat = summed / counts, duplicate_n = duplicate_n)
}

build_expression_qc <- function(cohort, gene_n_raw, gene_n_nonempty, gene_n_unique, duplicate_gene_n, sample_n, alias_applied_n) {
  data.frame(
    Cohort = cohort,
    Gene_N_Raw = gene_n_raw,
    Gene_N_NonEmpty = gene_n_nonempty,
    Gene_N_Unique = gene_n_unique,
    Duplicate_Gene_N = duplicate_gene_n,
    Alias_Applied_N = alias_applied_n,
    Sample_N = sample_n,
    stringsAsFactors = FALSE
  )
}

prepare_gse40419_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  gene_symbols <- normalize_gene_symbol(raw_df$gene)
  alias_res <- apply_gene_alias_map(gene_symbols)
  expr_df <- raw_df[, -(1:6), drop = FALSE]
  expr_mat <- as.matrix(data.frame(lapply(expr_df, as.numeric), check.names = FALSE))
  rownames(expr_mat) <- alias_res$gene_symbols
  keep <- nzchar(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  agg <- aggregate_duplicate_rows_mean(expr_mat)
  list(
    expr_mat = agg$expr_mat,
    qc_df = build_expression_qc(
      cohort = "GSE40419",
      gene_n_raw = nrow(raw_df),
      gene_n_nonempty = sum(keep),
      gene_n_unique = nrow(agg$expr_mat),
      duplicate_gene_n = agg$duplicate_n,
      sample_n = ncol(agg$expr_mat),
      alias_applied_n = alias_res$alias_applied_n
    )
  )
}

prepare_gse127165_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  ensembl_ids <- normalize_gene_symbol(sub("\\..*$", "", raw_df[[1]]))
  symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  gene_symbols <- normalize_gene_symbol(unname(symbol_map[ensembl_ids]))
  alias_res <- apply_gene_alias_map(gene_symbols)
  expr_df <- raw_df[, -1, drop = FALSE]
  expr_mat <- as.matrix(data.frame(lapply(expr_df, as.numeric), check.names = FALSE))
  rownames(expr_mat) <- alias_res$gene_symbols
  keep <- nzchar(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  agg <- aggregate_duplicate_rows_mean(expr_mat)
  list(
    expr_mat = agg$expr_mat,
    qc_df = build_expression_qc(
      cohort = "GSE127165",
      gene_n_raw = nrow(raw_df),
      gene_n_nonempty = sum(keep),
      gene_n_unique = nrow(agg$expr_mat),
      duplicate_gene_n = agg$duplicate_n,
      sample_n = ncol(agg$expr_mat),
      alias_applied_n = alias_res$alias_applied_n
    )
  )
}

prepare_gse96058_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = ",", data.table = FALSE)
  gene_symbols <- normalize_gene_symbol(raw_df[[1]])
  alias_res <- apply_gene_alias_map(gene_symbols)
  expr_df <- raw_df[, -1, drop = FALSE]
  expr_mat <- as.matrix(data.frame(lapply(expr_df, as.numeric), check.names = FALSE))
  rownames(expr_mat) <- alias_res$gene_symbols
  keep <- nzchar(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  agg <- aggregate_duplicate_rows_mean(expr_mat)
  list(
    expr_mat = agg$expr_mat,
    qc_df = build_expression_qc(
      cohort = "GSE96058",
      gene_n_raw = nrow(raw_df),
      gene_n_nonempty = sum(keep),
      gene_n_unique = nrow(agg$expr_mat),
      duplicate_gene_n = agg$duplicate_n,
      sample_n = ncol(agg$expr_mat),
      alias_applied_n = alias_res$alias_applied_n
    )
  )
}

calc_coverage_tables <- function(cohort, expr_mat, acid_genes, lm_imm_12_manifest, lm_imm_12_gene_long) {
  expr_genes <- rownames(expr_mat)

  acid_detail_df <- data.frame(
    Cohort = cohort,
    Gene = acid_genes,
    Signature = "Acid_Score",
    Present = acid_genes %in% expr_genes,
    stringsAsFactors = FALSE
  )

  acid_summary_df <- data.frame(
    Cohort = cohort,
    Signature = "Acid_Score",
    Signature_Name = "Acid_Score",
    Gene_Set_Size = length(acid_genes),
    Covered_Gene_N = sum(acid_detail_df$Present),
    Coverage_Ratio = sum(acid_detail_df$Present) / length(acid_genes),
    stringsAsFactors = FALSE
  )

  lm_imm_12_summary_df <- lm_imm_12_gene_long %>%
    group_by(feature_id, pathway) %>%
    summarise(
      Gene_Set_Size = n_distinct(gene),
      Covered_Gene_N = sum(unique(gene) %in% expr_genes),
      .groups = "drop"
    ) %>%
    left_join(
      lm_imm_12_manifest %>% dplyr::select(feature_id, short_label, pathway),
      by = c("feature_id", "pathway")
    ) %>%
    mutate(
      Cohort = cohort,
      Signature = feature_id,
      Signature_Name = pathway,
      Coverage_Ratio = Covered_Gene_N / Gene_Set_Size
    ) %>%
    dplyr::select(Cohort, Signature, Signature_Name, Gene_Set_Size, Covered_Gene_N, Coverage_Ratio)

  list(
    acid_detail_df = acid_detail_df,
    signature_summary_df = bind_rows(acid_summary_df, lm_imm_12_summary_df)
  )
}

final_selected_genes <- load_named_object(acid_gene_file, "final_selected_genes")
if (!exists("final_selected_genes")) {
  stop("final_selected_genes.rdata 中未找到 final_selected_genes。")
}
acid_genes <- unique(as.character(final_selected_genes$Gene))

cohort_manifest <- fread(manifest_file, sep = "\t", data.table = FALSE)
gse40419_annot <- fread(gse40419_annot_file, sep = "\t", data.table = FALSE)
gse127165_annot <- fread(gse127165_annot_file, sep = "\t", data.table = FALSE)
gse96058_annot <- fread(gse96058_annot_file, sep = "\t", data.table = FALSE)
lm_imm_12_manifest <- read_tsv(lm_imm_12_manifest_file, show_col_types = FALSE)
lm_imm_12_gene_long <- read_tsv(lm_imm_12_gene_long_file, show_col_types = FALSE)

lm_imm_12_gene_sets <- split(
  lm_imm_12_gene_long$gene,
  lm_imm_12_gene_long$feature_id
)
score_col_order <- c("Acid_Score", lm_imm_12_manifest$feature_id)
signatures_for_scoring <- c(list(Acid_Score = acid_genes), lm_imm_12_gene_sets)

annotation_map <- list(
  GSE40419 = gse40419_annot,
  GSE127165 = gse127165_annot,
  GSE96058 = gse96058_annot
)

expr_list <- list(
  GSE40419 = prepare_gse40419_expr(gse40419_file),
  GSE127165 = prepare_gse127165_expr(gse127165_file),
  GSE96058 = prepare_gse96058_expr(gse96058_file)
)

all_score_rows <- list()
expression_qc_list <- list()
coverage_detail_list <- list()
coverage_summary_list <- list()

for (cohort in names(expr_list)) {
  expr_mat <- expr_list[[cohort]]$expr_mat
  expression_qc_list[[cohort]] <- expr_list[[cohort]]$qc_df

  coverage_res <- calc_coverage_tables(
    cohort = cohort,
    expr_mat = expr_mat,
    acid_genes = acid_genes,
    lm_imm_12_manifest = lm_imm_12_manifest,
    lm_imm_12_gene_long = lm_imm_12_gene_long
  )
  coverage_detail_list[[cohort]] <- coverage_res$acid_detail_df
  coverage_summary_list[[cohort]] <- coverage_res$signature_summary_df

  expr_log2 <- log2(pmax(expr_mat, 0) + 1)
  ssgsea_param <- ssgseaParam(
    exprData = expr_log2,
    geneSets = signatures_for_scoring,
    normalize = TRUE
  )
  score_mat <- gsva(ssgsea_param, verbose = FALSE)
  score_df <- as.data.frame(t(score_mat), check.names = FALSE, stringsAsFactors = FALSE)
  score_df$Sample_ID <- rownames(score_df)

  annot_df <- annotation_map[[cohort]]
  score_df <- annot_df %>%
    left_join(score_df, by = "Sample_ID") %>%
    dplyr::select(
      any_of(c(
        "Cohort", "Cancer", "Sample_ID", "Sample_Type", "Pair_ID",
        "Is_Replicate", "Primary_Sample_ID", "Has_Clinical_Metadata"
      )),
      all_of(score_col_order)
    )

  all_score_rows[[cohort]] <- score_df

  write.table(
    score_df,
    file = file.path(out_dir, paste0(cohort, "_score_matrix.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

external_score_matrix <- bind_rows(all_score_rows) %>%
  arrange(match(Cohort, c("GSE40419", "GSE127165", "GSE96058")), Sample_ID)

write.table(
  external_score_matrix,
  file = file.path(out_dir, "external_score_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bind_rows(expression_qc_list),
  file = file.path(out_dir, "external_expression_preparation_qc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bind_rows(coverage_detail_list),
  file = file.path(out_dir, "acid_gene_coverage_detail.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bind_rows(coverage_summary_list),
  file = file.path(out_dir, "signature_coverage_qc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("外部 GEO score matrix 已保存到: ", out_dir)
