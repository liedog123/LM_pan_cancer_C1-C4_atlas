#!/usr/bin/env Rscript
# Purpose: Calculate the reported tumour-programme and within-signature gene correlations.
# Task: signature_correlations; step: s02
# Inputs (DCF fields): final_genes,tumor_scores,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task signature_correlations --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 02_correlate_final_nine_tumor_program.R --config PATH [--check-config]\n",
      "       Rscript 02_correlate_final_nine_tumor_program.R --help\n\n",
      "DCF fields: final_genes, tumor_scores, tpm_root, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary invocation computes final-nine within-cancer correlations and BH adjustment; no plotting.\n",
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
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link)) stop("Unresolvable symlink: ", path, call. = FALSE)
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
required_fields <- c("final_genes", "tumor_scores", "tpm_root", "output_dir")
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
final_file <- unname(paths["final_genes"])
score_file <- unname(paths["tumor_scores"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
# Inspect path metadata only; actual project keys are validated during an ordinary run.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
tpm_metadata_files <- list.files(tpm_root, pattern = "_mrna_expr_tpm\\.rdata$",
                                 recursive = TRUE, full.names = TRUE)
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L),
                         winslash = "/", mustWork = TRUE)
required_files <- c(final_file, score_file, tpm_real_paths)

# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
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

# Preserve the source raw TPM and complete TPM/score sample intersection.
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
load_object <- function(path, expected) {
  e <- new.env(parent = emptyenv())
  loaded <- load(path, envir = e)
  if (!expected %in% loaded) stop("Missing object ", expected, " in ", path)
  e[[expected]]
}
safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(c(n = sum(ok), r = NA, p = NA))
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  c(n = sum(ok), r = unname(z$estimate), p = z$p.value)
}
formal <- read.delim(final_file, check.names = FALSE)
gene_col <- if ("Gene" %in% names(formal)) "Gene" else names(formal)[1]
genes <- unique(as.character(formal[[gene_col]]))
stopifnot(length(genes) == 9L)

scores <- load_object(score_file, "tumor_gsva_list")
projects <- sort(names(scores))
stopifnot(length(projects) == 21L)

# Validate identifiers without changing the source sample definition.
if (anyNA(genes) || any(!nzchar(genes))) stop("Final gene names must be nonmissing and nonempty")
if (anyNA(projects) || any(!nzchar(projects)) || anyDuplicated(names(scores))) stop("Invalid project names")
cor_rows <- list()
for (project in projects) {
  tpm <- load_object(file.path(tpm_root, project, paste0(project, "_mrna_expr_tpm.rdata")), "mrna_expr_tpm")
  tpm <- as.matrix(tpm)
  storage.mode(tpm) <- "double"
  score_mat <- as.matrix(scores[[project]])
  score_vec <- as.numeric(score_mat[1, ])
  names(score_vec) <- colnames(score_mat)
  samples <- intersect(colnames(tpm), names(score_vec))
  if (is.null(colnames(tpm)) || is.null(rownames(tpm)) || anyDuplicated(colnames(tpm)) ||
      anyDuplicated(rownames(tpm)) || is.null(colnames(score_mat)) || anyDuplicated(colnames(score_mat)) ||
      nrow(score_mat) < 1L) stop("Input matrices require unique sample/gene identifiers in ", project)
  if (!all(genes %in% rownames(tpm))) stop("Final gene absent in ", project)

  project_cor <- do.call(rbind, lapply(genes, function(gene) {
    z <- safe_pearson(as.numeric(tpm[gene, samples]), score_vec[samples])
    data.frame(project = project, gene = gene, n_complete = as.integer(z["n"]),
               pearson_r = z["r"], pvalue = z["p"])
  }))
  project_cor$fdr <- p.adjust(project_cor$pvalue, method = "BH")
  project_cor$significance <- ifelse(project_cor$fdr < 0.001, "***",
    ifelse(project_cor$fdr < 0.01, "**", ifelse(project_cor$fdr < 0.05, "*", "")))
  cor_rows[[project]] <- project_cor

}
cor_long <- do.call(rbind, cor_rows)
cor_long$gene_order <- match(cor_long$gene, genes)
cor_long$cancer_order <- match(cor_long$project, projects)
cor_long <- cor_long[order(cor_long$gene_order, cor_long$cancer_order), ]

write_tsv(cor_long, file.path(out_dir, "panel_a_final_nine_complete_pearson.tsv"))
checks <- data.frame(check = c("genes_9", "cancers_21", "complete_a_189", "no_missing_a"),
  pass = c(length(genes) == 9L, length(projects) == 21L, nrow(cor_long) == 189L, !anyNA(cor_long$pearson_r)))
write_tsv(checks, file.path(out_dir, "machine_checks.tsv"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
if (!all(checks$pass)) stop("Final-nine correlation QA failed")
message("Final-nine correlation table exported.")
