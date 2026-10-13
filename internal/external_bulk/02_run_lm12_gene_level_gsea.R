#!/usr/bin/env Rscript
# Continue cohort-specific external Acid-linked LM12 direction validation.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 02_run_lm12_gene_level_gsea.R --config PATH [--check-config]\n",
      "       Rscript 02_run_lm12_gene_level_gsea.R --help\n\n",
      "DCF fields: expression_dir, score_dir, lm12_manifest, lm12_genes, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation runs cohort-specific High versus Low gene-level LM12 GSEA.\n",
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
required_fields <- c("expression_dir", "score_dir", "lm12_manifest", "lm12_genes", "output_dir")
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
score_dir <- unname(paths["score_dir"])
manifest_file <- unname(paths["lm12_manifest"])
gene_long_file <- unname(paths["lm12_genes"])
out_dir <- unname(paths["output_dir"])
expr_file_map <- c(
  GSE40419 = file.path(geo_source_dir, "GSE40419_LC-87_RPKM_expression.txt.gz"),
  GSE127165 = file.path(geo_source_dir, "GSE127165_114_samples_RNA_seq_FPKM.txt.gz"),
  GSE96058 = file.path(geo_source_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
)
score_file_map <- setNames(
  file.path(score_dir, paste0(names(expr_file_map), "_score_matrix.tsv")),
  names(expr_file_map)
)
required_files <- c(manifest_file, gene_long_file, expr_file_map, score_file_map)

# Inspect required fixed filenames without reading or decompressing their contents.
input_dirs <- c(geo_source_dir, score_dir)
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

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(limma)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
set.seed(20260804)
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

gene_alias_map <- c("PKM2" = "PKM")

normalize_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  hit <- x %in% names(gene_alias_map)
  x[hit] <- unname(gene_alias_map[x[hit]])
  x
}

aggregate_duplicate_rows_mean <- function(expr_mat) {
  if (!anyDuplicated(rownames(expr_mat))) return(expr_mat)
  summed <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(
    matrix(1, nrow = nrow(expr_mat), ncol = ncol(expr_mat)),
    group = rownames(expr_mat), reorder = FALSE
  )
  summed / counts
}

prepare_gse40419_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  gene_symbols <- normalize_gene_symbol(raw_df$gene)
  expr_mat <- as.matrix(data.frame(lapply(raw_df[, -(1:6), drop = FALSE], as.numeric), check.names = FALSE))
  rownames(expr_mat) <- gene_symbols
  expr_mat <- expr_mat[nzchar(rownames(expr_mat)), , drop = FALSE]
  aggregate_duplicate_rows_mean(expr_mat)
}

prepare_gse127165_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  ensembl_ids <- sub("\\..*$", "", as.character(raw_df[[1]]))
  symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
  )
  gene_symbols <- normalize_gene_symbol(unname(symbol_map[ensembl_ids]))
  expr_mat <- as.matrix(data.frame(lapply(raw_df[, -1, drop = FALSE], as.numeric), check.names = FALSE))
  rownames(expr_mat) <- gene_symbols
  expr_mat <- expr_mat[nzchar(rownames(expr_mat)), , drop = FALSE]
  aggregate_duplicate_rows_mean(expr_mat)
}

prepare_gse96058_expr <- function(path) {
  # 该矩阵较大，保持data.table数值列并尽早释放原对象，避免data.frame/lapply产生多份全矩阵副本。
  raw_dt <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = ",", data.table = TRUE)
  gene_symbols <- normalize_gene_symbol(raw_dt[[1]])
  raw_dt[[1]] <- NULL
  expr_mat <- as.matrix(raw_dt)
  rm(raw_dt)
  gc(verbose = FALSE)
  rownames(expr_mat) <- gene_symbols
  expr_mat <- expr_mat[nzchar(rownames(expr_mat)), , drop = FALSE]
  aggregate_duplicate_rows_mean(expr_mat)
}

expr_loader_map <- list(
  GSE40419 = prepare_gse40419_expr,
  GSE127165 = prepare_gse127165_expr,
  GSE96058 = prepare_gse96058_expr
)

manifest_df <- read_tsv(manifest_file, show_col_types = FALSE) %>% arrange(display_order)
gene_long_df <- read_tsv(gene_long_file, show_col_types = FALSE)
if (
  nrow(manifest_df) != 12 ||
    anyDuplicated(manifest_df$feature_id) > 0 ||
    !all(manifest_df$display_order == seq_len(12))
) {
  stop("LM12 manifest必须为12条唯一feature且display_order严格为1..12。")
}

lm12_gene_sets <- split(gene_long_df$gene, gene_long_df$feature_id)
lm12_gene_sets <- lapply(lm12_gene_sets[manifest_df$feature_id], function(x) {
  unique(normalize_gene_symbol(x[nzchar(x)]))
})

assign_high_low <- function(score_df) {
  score_df <- score_df %>% arrange(desc(Acid_Score), Sample_ID)
  high_n <- ceiling(nrow(score_df) / 2)
  score_df$Acid_Group <- c(rep("High_Acid", high_n), rep("Low_Acid", nrow(score_df) - high_n))
  score_df
}

run_one_analysis <- function(cohort, expr_log2, score_df, analysis_set) {
  common_samples <- intersect(score_df$Sample_ID, colnames(expr_log2))
  score_work <- score_df %>%
    filter(Sample_ID %in% common_samples, Sample_Type == "Tumor", is.finite(Acid_Score)) %>%
    distinct(Sample_ID, .keep_all = TRUE)
  score_work <- assign_high_low(score_work)

  expr_work <- expr_log2[, score_work$Sample_ID, drop = FALSE]
  group_factor <- factor(score_work$Acid_Group, levels = c("Low_Acid", "High_Acid"))
  design <- model.matrix(~ group_factor)
  fit <- eBayes(lmFit(expr_work, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$Gene <- rownames(tt)

  rank_df <- tt %>%
    transmute(Gene, logFC, AveExpr, moderated_t = t, P_Value = P.Value, FDR = adj.P.Val, B) %>%
    filter(nzchar(Gene), is.finite(moderated_t)) %>%
    arrange(desc(moderated_t), Gene)
  if (anyDuplicated(rank_df$Gene) > 0) stop(cohort, ": ranked genes存在重复。")
  ranks <- rank_df$moderated_t
  names(ranks) <- rank_df$Gene
  ranks <- sort(ranks, decreasing = TRUE)

  coverage_df <- manifest_df %>%
    transmute(
      display_order, feature_id, short_label, pathway, direction_label, functional_group,
      gene_set_size = vapply(lm12_gene_sets[feature_id], length, integer(1)),
      covered_gene_n = vapply(lm12_gene_sets[feature_id], function(gs) sum(gs %in% names(ranks)), integer(1)),
      coverage_ratio = covered_gene_n / gene_set_size
    )

  fgsea_res <- fgseaMultilevel(
    pathways = lm12_gene_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  ) %>%
    as.data.frame()
  fgsea_res$leadingEdge <- vapply(fgsea_res$leadingEdge, paste, collapse = ",", character(1))

  result_df <- manifest_df %>%
    dplyr::select(display_order, feature_id, short_label, pathway, direction_label, functional_group) %>%
    left_join(fgsea_res, by = c("feature_id" = "pathway")) %>%
    left_join(coverage_df %>% dplyr::select(feature_id, gene_set_size, covered_gene_n, coverage_ratio), by = "feature_id") %>%
    mutate(
      Cohort = cohort,
      Cancer = unique(score_work$Cancer),
      analysis_set = analysis_set,
      High_Acid_N = sum(score_work$Acid_Group == "High_Acid"),
      Low_Acid_N = sum(score_work$Acid_Group == "Low_Acid"),
      expected_nes_sign = ifelse(direction_label == "High_Acid_Positive", 1, -1),
      direction_consistent = sign(NES) == expected_nes_sign
    ) %>%
    dplyr::select(
      Cohort, Cancer, analysis_set, display_order, feature_id, short_label, pathway,
      direction_label, functional_group, High_Acid_N, Low_Acid_N,
      gene_set_size, covered_gene_n, coverage_ratio,
      size, ES, NES, pval, padj, log2err, leadingEdge,
      expected_nes_sign, direction_consistent
    ) %>%
    arrange(display_order)

  if (nrow(result_df) != 12 || any(!is.finite(result_df$NES))) {
    stop(cohort, " / ", analysis_set, ": 未产生12条有限NES结果。")
  }

  prefix <- paste0(cohort, "_lm12_gene_level_gsea", ifelse(analysis_set == "primary", "", paste0("_", analysis_set)))
  write_tsv(score_work %>% dplyr::select(Cohort, Cancer, Sample_ID, Acid_Score, Acid_Group, any_of(c("Is_Replicate", "Primary_Sample_ID"))), file.path(out_dir, paste0(prefix, "_sample_groups.tsv")))
  write_tsv(coverage_df, file.path(out_dir, paste0(prefix, "_coverage.tsv")))
  write_tsv(result_df, file.path(out_dir, paste0(prefix, "_results.tsv")))
  rank_con <- gzfile(file.path(out_dir, paste0(prefix, "_ranked_t.tsv.gz")), open = "wt")
  write.table(rank_df, rank_con, sep = "\t", row.names = FALSE, quote = FALSE)
  close(rank_con)
  saveRDS(
    list(
      cohort = cohort,
      cancer = unique(score_work$Cancer),
      analysis_set = analysis_set,
      sample_groups = score_work,
      pathway_order = manifest_df$feature_id,
      pathway_gene_sets = lm12_gene_sets,
      fgsea_results = result_df,
      ranked_stats = ranks
    ),
    file.path(out_dir, paste0(prefix, "_objects.rds"))
  )

  list(results = result_df, coverage = coverage_df, sample_groups = score_work, rank_n = length(ranks))
}

primary_results <- list()
qa_rows <- list()
sensitivity_results <- list()

for (cohort in names(expr_file_map)) {
  message("读取并处理: ", cohort)
  expr_mat <- expr_loader_map[[cohort]](expr_file_map[[cohort]])
  expr_log2 <- log2(pmax(expr_mat, 0) + 1)
  rm(expr_mat)
  gc(verbose = FALSE)

  score_df <- read_tsv(score_file_map[[cohort]], show_col_types = FALSE)
  score_primary <- score_df
  if (cohort == "GSE96058") {
    score_primary <- score_df %>% filter(!as.logical(Is_Replicate))
  }

  primary <- run_one_analysis(cohort, expr_log2, score_primary, "primary")
  primary_results[[cohort]] <- primary$results
  qa_rows[[cohort]] <- tibble(
    Cohort = cohort,
    Cancer = unique(primary$results$Cancer),
    analysis_set = "primary",
    Sample_N = nrow(primary$sample_groups),
    High_Acid_N = sum(primary$sample_groups$Acid_Group == "High_Acid"),
    Low_Acid_N = sum(primary$sample_groups$Acid_Group == "Low_Acid"),
    Ranked_Gene_N = primary$rank_n,
    Pathway_N = nrow(primary$results),
    Direction_Consistent_N = sum(primary$results$direction_consistent),
    FDR_lt_0_05_N = sum(primary$results$padj < 0.05)
  )

  if (cohort == "GSE96058") {
    sensitivity <- run_one_analysis(cohort, expr_log2, score_df, "all_samples_sensitivity")
    sensitivity_results[[cohort]] <- sensitivity$results
    qa_rows[[paste0(cohort, "_all")]] <- tibble(
      Cohort = cohort,
      Cancer = unique(sensitivity$results$Cancer),
      analysis_set = "all_samples_sensitivity",
      Sample_N = nrow(sensitivity$sample_groups),
      High_Acid_N = sum(sensitivity$sample_groups$Acid_Group == "High_Acid"),
      Low_Acid_N = sum(sensitivity$sample_groups$Acid_Group == "Low_Acid"),
      Ranked_Gene_N = sensitivity$rank_n,
      Pathway_N = nrow(sensitivity$results),
      Direction_Consistent_N = sum(sensitivity$results$direction_consistent),
      FDR_lt_0_05_N = sum(sensitivity$results$padj < 0.05)
    )
    rm(sensitivity)
    gc(verbose = FALSE)
  }

  rm(expr_log2, score_df, score_primary)
  gc(verbose = FALSE)
}

combined_primary <- bind_rows(primary_results) %>%
  arrange(match(Cohort, names(expr_file_map)), display_order)
combined_qa <- bind_rows(qa_rows)

if (
  nrow(combined_primary) != 36 ||
    anyDuplicated(combined_primary[c("Cohort", "feature_id")]) > 0 ||
    !all(combined_qa$Pathway_N == 12)
) {
  stop("三队列LM12 GSEA最终QA失败。")
}

write_tsv(combined_primary, file.path(out_dir, "external_lm_imm_12_gene_level_gsea_results.tsv"))
write_tsv(combined_qa, file.path(out_dir, "external_lm_imm_12_gene_level_gsea_qa.tsv"))
if (length(sensitivity_results) > 0) {
  write_tsv(bind_rows(sensitivity_results), file.path(out_dir, "external_lm_imm_12_gene_level_gsea_sensitivity_results.tsv"))
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "external_lm_imm_12_gene_level_gsea_sessionInfo.txt"))

message("外部LM12 gene-level GSEA完成。")
message("主分析结果: ", nrow(combined_primary), " rows")
message(paste0(combined_qa$Cohort, "/", combined_qa$analysis_set, ": n=", combined_qa$Sample_N, ", direction=", combined_qa$Direction_Consistent_N, "/12, FDR<0.05=", combined_qa$FDR_lt_0_05_N, collapse = "\n"))
message("输出目录: ", out_dir)
