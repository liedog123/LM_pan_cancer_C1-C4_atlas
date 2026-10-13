#!/usr/bin/env Rscript
# Compute final-nine correlations from existing study intermediates.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 03_correlate_final_nine_gene_pairs.R --config PATH [--check-config]\n",
      "       Rscript 03_correlate_final_nine_gene_pairs.R --help\n\n",
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
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(c(n = n, rho = NA, p = NA))
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = n, rho = unname(z$estimate), p = z$p.value)
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
pairs <- combn(genes, 2, simplify = FALSE)
pair_levels <- vapply(pairs, function(gp) paste(gp, collapse = " / "), character(1))

rows <- list()
for (project in projects) {
  tpm <- load_object(file.path(tpm_root, project, paste0(project, "_mrna_expr_tpm.rdata")), "mrna_expr_tpm")
  tpm <- as.matrix(tpm)
  storage.mode(tpm) <- "double"
  score_mat <- as.matrix(scores[[project]])
  samples <- intersect(colnames(tpm), colnames(score_mat))
  if (is.null(colnames(tpm)) || is.null(rownames(tpm)) || anyDuplicated(colnames(tpm)) ||
      anyDuplicated(rownames(tpm)) || is.null(colnames(score_mat)) || anyDuplicated(colnames(score_mat)))
    stop("Input matrices require unique sample/gene identifiers in ", project)
  if (!all(genes %in% rownames(tpm))) stop("Final gene absent in ", project)
  d <- do.call(rbind, lapply(pairs, function(gp) {
    z <- safe_spearman(as.numeric(tpm[gp[1], samples]), as.numeric(tpm[gp[2], samples]))
    data.frame(project = project, gene1 = gp[1], gene2 = gp[2], gene_pair = paste(gp, collapse = " / "),
               n_complete = as.integer(z["n"]), spearman_rho = z["rho"], pvalue = z["p"])
  }))
  d$fdr <- p.adjust(d$pvalue, method = "BH")
  d$significance <- ifelse(d$fdr < 0.001, "***", ifelse(d$fdr < 0.01, "**", ifelse(d$fdr < 0.05, "*", "")))
  rows[[project]] <- d
}
pair_long <- do.call(rbind, rows)
pair_long$pair_order <- match(pair_long$gene_pair, pair_levels)
pair_long$cancer_order <- match(pair_long$project, projects)
pair_long <- pair_long[order(pair_long$pair_order, pair_long$cancer_order), ]
write_tsv(pair_long, file.path(out_dir, "panel_b_final_nine_pairwise_spearman_complete.tsv"))

checks <- data.frame(check = c("genes_9", "gene_pairs_36", "cancers_21", "complete_pair_grid_756", "no_missing_rho"),
  pass = c(length(genes) == 9L, length(pair_levels) == 36L, length(projects) == 21L,
           nrow(pair_long) == 756L, !anyNA(pair_long$spearman_rho)))
write_tsv(checks, file.path(out_dir, "machine_checks.tsv"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
if (!all(checks$pass)) stop("Final-nine correlation QA failed")
message("Final-nine correlation table exported.")
