#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 14_run_melanoma_pseudobulk_programs.R --config PATH [--check-config]\n       Rscript 14_run_melanoma_pseudobulk_programs.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { usage(); quit(status = 0L) }
if (!(length(args) %in% c(2L,3L)) || args[1L] != "--config" ||
    !nzchar(args[2L]) || (length(args)==3L && args[3L] != "--check-config")) {
  usage(); stop("Expected --config PATH [--check-config].")
}
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link) && !file.exists(path) && !dir.exists(path)) stop("Dangling symbolic link: ", path)
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

config_path <- normalizePath(args[2L], mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path,4L)!=0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size<=0) stop("Invalid configuration.")
config_lines <- readLines(config_path, warn = FALSE)
config_keys <- trimws(sub(":.*$", "", config_lines[
  grepl("^[^[:space:]#][^:]*:", config_lines)]))
if(anyDuplicated(config_keys)) stop("DCF fields must not repeat.")
cfg <- read.dcf(config_path)
fields <- c("contract", "stage3_gate", "source_object", "frozen_cell_metadata", "sample_labels", "readiness", "catalog", "membership_full", "membership_nondefinition", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("contract", "stage3_gate", "source_object", "frozen_cell_metadata", "sample_labels", "readiness", "catalog", "membership_full", "membership_nondefinition")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),input_dirs,config_path,package_root)
if(out_dir=="/" || out_dir==normalizePath(path.expand("~")) ||
   any(vapply(protected,function(p) is_within(p,out_dir)||is_within(out_dir,p),logical(1L)))) stop("Output must be independent of inputs/config.")
if(file.exists(out_dir) && !dir.exists(out_dir)) stop("Output is not a directory.")
if(dir.exists(out_dir) && (file.access(out_dir,7L)!=0L || length(list.files(out_dir,all.files=TRUE,no..=TRUE)))) stop("Output must be writable and empty.")
parent <- out_dir
while(!dir.exists(parent)) {
 if(file.exists(parent)) stop("Output has a non-directory ancestor.")
 parent <- dirname(parent)
}
if(file.access(parent,3L)!=0L) stop("No writable searchable output parent.")
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

# ---- Retained source body; only documented path/read adapters ----
cfg_paths <- paths


Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(edgeR)
  library(GSVA)
})

options(stringsAsFactors = FALSE)

SCRIPT_NAME <- "7.3scRNA_c1c4_state_closure/scripts/04a_build_melanoma_pseudobulk_and_program_scores.R"
CONTRACT_VERSION <- "STAGE4_MELANOMA_C3_VALIDATION_FROZEN_V0.1_20260722"
CELLTYPE_TO_COMPARTMENT <- c(
  "CD8 T" = "T/NK", "CD4 T" = "T/NK",
  "Monocyte/Macrophage" = "Myeloid", "cDC" = "Myeloid",
  "Fibroblast" = "Stromal/CAF", "Tumor" = "Malignant"
)

log_event <- function(level, msg, fields = list()) {
  clean <- function(x) gsub('"', "'", as.character(x), fixed = TRUE)
  extras <- if (length(fields) == 0L) "" else paste0(
    " ", paste(sprintf('%s="%s"', names(fields), vapply(fields, clean, character(1L))), collapse = " ")
  )
  cat(sprintf('[%s] %s script="%s" msg="%s"%s\n',
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, SCRIPT_NAME, clean(msg), extras))
  flush.console()
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = TRUE, na = "NA")

safe_md5 <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

file_record <- function(role, path, critical = TRUE) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  data.table(
    role = role, path = path, critical = critical, exists = exists,
    bytes = if (exists) as.numeric(info$size) else NA_real_,
    nonempty = exists && !dir.exists(path) && isTRUE(info$size > 0),
    modified_time = if (exists) format(info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
    md5 = safe_md5(path)
  )
}

get_counts <- function(object) {
  mat <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", layer = "counts"),
    error = function(e) tryCatch(
      suppressWarnings(Seurat::GetAssayData(object, assay = "RNA", slot = "counts")),
      error = function(e2) NULL
    )
  )
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) stop("RNA counts matrix is unavailable")
  if (!inherits(mat, "sparseMatrix")) mat <- as(mat, "dgCMatrix")
  mat
}

zscore_rows <- function(mat) {
  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, sd, na.rm = TRUE)
  out <- sweep(mat, 1L, means, "-")
  out <- sweep(out, 1L, sds, "/")
  out[!is.finite(out)] <- NA_real_
  out
}

aggregate_celltype_counts <- function(counts, meta, cell_type) {
  cell_type_value <- as.character(cell_type)[1L]
  md <- meta[get("cell_type") == cell_type_value]
  if (nrow(md) == 0L) stop("No cells for cell type: ", cell_type)
  idx <- match(md$cell_id, colnames(counts))
  if (anyNA(idx)) stop("Frozen cell metadata does not match counts for: ", cell_type)
  sample_levels <- sort(unique(md$sample_ID_long))
  sample_factor <- factor(md$sample_ID_long, levels = sample_levels)
  design <- sparse.model.matrix(~ 0 + sample_factor)
  colnames(design) <- sample_levels
  pb <- counts[, idx, drop = FALSE] %*% design
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sample_levels
  pb
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  count_dir <- file.path(outdir, "counts")
  dir.create(count_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- as.list(cfg_paths[c("contract", "stage3_gate", "source_object", "frozen_cell_metadata", "sample_labels", "readiness", "catalog", "membership_full", "membership_nondefinition")])

  log_event("START", "build melanoma sample-by-native-cell-type pseudobulk and frozen program scores", list(output = outdir))
  input_manifest <- rbindlist(Map(function(role, path) file_record(role, path, TRUE), names(paths), unname(unlist(paths))))
  write_tsv(input_manifest, file.path(outdir, "input_file_manifest.tsv"))
  if (nrow(input_manifest[!exists | !nonempty]) > 0L) stop("Stage 4 score inputs are incomplete")

  contract <- fread(paths$contract)
  if (!any(contract$contract_item == "contract_version" & contract$frozen_value == CONTRACT_VERSION)) {
    stop("Stage 4 contract version mismatch")
  }
  stage3 <- fread(paths$stage3_gate)
  if (stage3$status[[1L]] != "PASS" || !isTRUE(stage3$stage4_authorized[[1L]])) {
    stop("Stage 3 did not authorize Stage 4")
  }

  cell_meta <- fread(paths$frozen_cell_metadata)
  cell_meta <- cell_meta[flag_keep_for_scoring %in% TRUE & cell_type %in% names(CELLTYPE_TO_COMPARTMENT)]
  labels <- fread(paths$sample_labels)
  label_cols <- c("sample_ID_long", "patient", "sample_C_subtype", "state", "response_short",
                  "treatment_group", "subtype", "tissue", "day")
  if (!all(label_cols %in% names(labels))) stop("Melanoma label columns are incomplete")
  labels <- unique(labels[, ..label_cols], by = "sample_ID_long")

  log_event("CHECKPOINT", "load frozen melanoma Seurat object", list(path = paths$source_object))
  object <- readRDS(paths$source_object)
  counts <- get_counts(object)
  rm(object)
  gc(verbose = FALSE)
  if (!all(cell_meta$cell_id %in% colnames(counts))) stop("Frozen metadata cells are missing from source object")

  catalog <- fread(paths$catalog)
  membership_full <- fread(paths$membership_full)
  membership_nondefinition <- fread(paths$membership_nondefinition)
  score_rows <- list()
  coverage_rows <- list()
  sample_rows <- list()
  count_manifest_rows <- list()

  for (cell_type in names(CELLTYPE_TO_COMPARTMENT)) {
    cell_type_value <- as.character(cell_type)[1L]
    compartment_value <- unname(CELLTYPE_TO_COMPARTMENT[[cell_type_value]])
    log_event("CHECKPOINT", "aggregate and score native cell type", list(cell_type = cell_type_value, compartment = compartment_value))
    pb <- aggregate_celltype_counts(counts, cell_meta, cell_type_value)
    cell_counts <- cell_meta[get("cell_type") == cell_type_value, .(n_cells = .N), by = sample_ID_long]
    sample_sheet <- merge(data.table(sample_ID_long = colnames(pb)), cell_counts,
                          by = "sample_ID_long", all.x = TRUE)
    sample_sheet <- merge(sample_sheet, labels, by = "sample_ID_long", all.x = TRUE)
    sample_sheet[, `:=`(
      cell_type = cell_type_value,
      compartment = compartment_value,
      unit_id = paste(safe_name(cell_type_value), sample_ID_long, sep = "||"),
      library_size = as.numeric(colSums(pb)),
      pre_post = fifelse(state == "PRE", "Pre", "Post"),
      eligible_ge20 = n_cells >= 20L,
      eligible_ge30 = n_cells >= 30L,
      eligible_ge50 = n_cells >= 50L
    )]
    if (anyNA(sample_sheet$patient) || anyNA(sample_sheet$sample_C_subtype)) {
      stop("Sample labels failed to join for: ", cell_type_value)
    }
    idx_all <- match(sample_sheet$sample_ID_long, colnames(pb))
    pb <- pb[, idx_all, drop = FALSE]
    colnames(pb) <- sample_sheet$unit_id

    count_path <- file.path(count_dir, paste0(safe_name(cell_type_value), "_pseudobulk_counts.rds"))
    saveRDS(list(counts = pb, sample_sheet = sample_sheet), count_path, compress = TRUE)
    count_manifest_rows[[length(count_manifest_rows) + 1L]] <- data.table(
      cell_type = cell_type_value, compartment = compartment_value, file = count_path,
      n_genes = nrow(pb), n_units = ncol(pb), n_units_ge20 = sum(sample_sheet$eligible_ge20)
    )
    sample_rows[[length(sample_rows) + 1L]] <- sample_sheet

    eligible <- sample_sheet[eligible_ge20 == TRUE]
    idx <- match(eligible$unit_id, colnames(pb))
    y <- DGEList(counts = pb[, idx, drop = FALSE])
    y <- calcNormFactors(y)
    logcpm <- cpm(y, log = TRUE, prior.count = 2)

    program_ids <- catalog[get("compartment") == compartment_value, program_id]
    full_sets <- lapply(program_ids, function(id) intersect(membership_full[program_id == id, unique(gene)], rownames(logcpm)))
    nondef_sets <- lapply(program_ids, function(id) intersect(membership_nondefinition[program_id == id, unique(gene)], rownames(logcpm)))
    names(full_sets) <- paste(program_ids, "FULL", sep = "||")
    names(nondef_sets) <- paste(program_ids, "NONDEFINITION", sep = "||")
    gene_sets <- c(full_sets, nondef_sets)
    if (any(lengths(gene_sets) < 3L)) stop("Fewer than 3 scored genes for a frozen program in: ", cell_type_value)

    score_mat <- gsva(
      ssgseaParam(exprData = logcpm, geneSets = gene_sets, minSize = 1L,
                  alpha = 0.25, normalize = TRUE),
      verbose = FALSE
    )
    score_z <- zscore_rows(score_mat)
    raw_long <- melt(as.data.table(t(score_mat), keep.rownames = "unit_id"), id.vars = "unit_id",
                     variable.name = "score_id", value.name = "score_raw")
    z_long <- melt(as.data.table(t(score_z), keep.rownames = "unit_id"), id.vars = "unit_id",
                   variable.name = "score_id", value.name = "score_z_global")
    long <- merge(raw_long, z_long, by = c("unit_id", "score_id"), all = TRUE)
    long[, c("program_id", "score_layer") := tstrsplit(score_id, "||", fixed = TRUE)]
    long[, score_id := NULL]
    long <- merge(long, eligible, by = "unit_id", all.x = TRUE, sort = FALSE)
    long <- merge(long, catalog[, .(program_id, program_order, biological_family, evidence_use,
                                    higher_score_interpretation, allowed_claim)],
                  by = "program_id", all.x = TRUE, sort = FALSE)
    long[, contract_version := CONTRACT_VERSION]
    setorder(long, program_order, score_layer, sample_ID_long)
    score_rows[[length(score_rows) + 1L]] <- long

    coverage_rows[[length(coverage_rows) + 1L]] <- rbindlist(lapply(names(gene_sets), function(id) {
      parts <- strsplit(id, "||", fixed = TRUE)[[1L]]
      data.table(cell_type = cell_type_value, compartment = compartment_value,
                 program_id = parts[[1L]], score_layer = parts[[2L]],
                 n_genes_scored = length(gene_sets[[id]]), n_units_scored = ncol(logcpm))
    }))
    rm(pb, y, logcpm, score_mat, score_z)
    gc(verbose = FALSE)
  }
  rm(counts)
  gc(verbose = FALSE)

  scores <- rbindlist(score_rows, fill = TRUE)
  coverage <- rbindlist(coverage_rows, fill = TRUE)
  sample_sheet_all <- rbindlist(sample_rows, fill = TRUE)
  count_manifest <- rbindlist(count_manifest_rows, fill = TRUE)
  fwrite(scores, file.path(outdir, "melanoma_pseudobulk_program_scores.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  write_tsv(coverage, file.path(outdir, "melanoma_program_score_coverage.tsv"))
  write_tsv(sample_sheet_all, file.path(outdir, "melanoma_pseudobulk_sample_sheet.tsv"))
  write_tsv(count_manifest, file.path(outdir, "melanoma_pseudobulk_counts_manifest.tsv"))

  concordance <- dcast(scores[, .(unit_id, cell_type, program_id, score_layer, score_raw)],
                       unit_id + cell_type + program_id ~ score_layer, value.var = "score_raw")
  concordance <- concordance[, .(
    full_nondefinition_spearman = cor(FULL, NONDEFINITION, method = "spearman", use = "pairwise.complete.obs"),
    max_abs_difference = max(abs(FULL - NONDEFINITION), na.rm = TRUE)
  ), by = .(cell_type, program_id)]
  concordance <- merge(catalog[, .(program_id, program_order, compartment, definition_overlap_fraction)],
                       concordance, by = "program_id", all.y = TRUE)
  setorder(concordance, program_order, cell_type)
  write_tsv(concordance, file.path(outdir, "melanoma_full_nondefinition_score_concordance.tsv"))

  expected_units <- sum(vapply(names(CELLTYPE_TO_COMPARTMENT), function(ct) {
    catalog[compartment == CELLTYPE_TO_COMPARTMENT[[ct]], .N]
  }, integer(1L)))
  checks <- data.table(
    check_id = c(
      "stage3_authorized", "six_predefined_celltypes", "all_source_cells_match",
      "all_sample_labels_joined", "all_expected_celltype_program_units_two_layers",
      "all_scores_finite", "all_count_bundles_nonempty", "qc_program_restricted"
    ),
    passed = c(
      TRUE, uniqueN(sample_sheet_all$cell_type) == 6L, TRUE,
      !anyNA(sample_sheet_all$patient) && !anyNA(sample_sheet_all$sample_C_subtype),
      uniqueN(coverage[, paste(cell_type, program_id)]) == expected_units &&
        all(coverage[, .N, by = .(cell_type, program_id)]$N == 2L),
      all(is.finite(scores$score_raw)) && all(is.finite(scores$score_z_global)),
      nrow(count_manifest) == 6L && all(file.exists(count_manifest$file)) && all(file.info(count_manifest$file)$size > 0),
      nrow(catalog[program_id == "MALIGNANT_CELL_CYCLE_QC" & evidence_use == "QC_ONLY"]) == 1L
    )
  )
  gate_pass <- all(checks$passed)
  gate <- data.table(
    gate = "STAGE4_MELANOMA_PSEUDOBULK_AND_SCORE_FREEZE",
    status = if (gate_pass) "PASS" else "FAIL",
    contract_version = CONTRACT_VERSION,
    n_celltypes = uniqueN(sample_sheet_all$cell_type),
    n_sample_celltype_units = nrow(sample_sheet_all),
    n_units_ge20 = sum(sample_sheet_all$eligible_ge20),
    n_celltype_program_units = expected_units,
    n_score_rows = nrow(scores),
    program_modeling_authorized = gate_pass,
    boundary = "Label-independent melanoma pseudobulk scores only; no C3 comparison has been tested at this step."
  )
  write_tsv(checks, file.path(outdir, "stage4_score_checks.tsv"))
  write_tsv(gate, file.path(outdir, "stage4_score_gate.tsv"))

  output_files <- list.files(outdir, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[file.info(output_files)$isdir == FALSE]
  output_manifest <- rbindlist(lapply(output_files, function(path) file_record(basename(path), path, TRUE)))
  write_tsv(output_manifest, file.path(outdir, "output_manifest.tsv"))
  log_event("DONE", "melanoma pseudobulk and frozen program scores completed", list(
    status = gate$status, n_score_rows = nrow(scores),
    elapsed_minutes = round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2)
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
