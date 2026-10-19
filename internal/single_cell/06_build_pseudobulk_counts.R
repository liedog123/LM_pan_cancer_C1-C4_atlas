#!/usr/bin/env Rscript
# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Task: single_cell_pseudobulk; step: s06
# Inputs (DCF fields): seurat_object,cell_metadata,sample_mapping,exclusion_genes,g1_summary
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 06_build_pseudobulk_counts.R --config PATH [--check-config]\n       Rscript 06_build_pseudobulk_counts.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("seurat_object", "cell_metadata", "sample_mapping", "exclusion_genes", "g1_summary", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "cell_metadata", "sample_mapping", "exclusion_genes", "g1_summary")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character()
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


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(SummarizedExperiment)
  library(S4Vectors)
})

options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02a_build_pseudobulk_counts.R"

log_event <- function(level, msg, ...) {
  extras <- list(...)
  extra_text <- if (length(extras) == 0L) "" else paste(
    vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)),
    collapse = " "
  )
  cat(sprintf("[%s] %s script=%s msg=%s%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg),
              if (nzchar(extra_text)) " " else "", extra_text))
  flush.console()
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  counts_dir <- file.path(outdir, "counts")
  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(seurat_rds = unname(paths["seurat_object"]), cell_metadata = unname(paths["cell_metadata"]),
                sample_mapping = unname(paths["sample_mapping"]), exclusion_genes = unname(paths["exclusion_genes"]),
                g1_summary = unname(paths["g1_summary"]))
  missing <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]
  if (length(missing) > 0L) stop("Missing inputs: ", paste(missing, collapse = ", "))

  log_event("START", "build sample-by-compartment pseudobulk raw counts",
            input = paths$seurat_rds, output = outdir)

  g1 <- fread(paths$g1_summary)
  if (!identical(g1$status[[1]], "PASS")) stop("G1 status is not PASS")
  cell_meta <- fread(paths$cell_metadata)
  mapping <- fread(paths$sample_mapping)
  exclusion <- fread(paths$exclusion_genes)

  log_event("CHECKPOINT", "load Seurat object and extract RNA counts")
  obj <- readRDS(paths$seurat_rds)
  if (!"RNA" %in% names(obj@assays) || !"counts" %in% Layers(obj[["RNA"]])) {
    stop("RNA/counts layer is unavailable")
  }
  counts <- LayerData(obj[["RNA"]], layer = "counts")
  genes <- rownames(counts)
  cells <- colnames(counts)
  if (!inherits(counts, "sparseMatrix")) stop("RNA/counts is not sparse")
  if (anyDuplicated(genes)) stop("Duplicated feature names in RNA/counts")
  ord <- match(cells, cell_meta$cell_id)
  if (anyNA(ord) || length(ord) != ncol(counts)) stop("Cell metadata does not match counts columns")
  cell_meta <- cell_meta[ord]
  if (!identical(cell_meta$cell_id, cells)) stop("Cell order verification failed")
  cell_library <- as.numeric(Matrix::colSums(counts))
  rm(obj)
  invisible(gc())

  non_other <- cell_meta$project_broad_compartment != "Other"
  if (!any(non_other)) stop("No non-Other cells")
  cell_meta[, cell_library_size := cell_library]
  cell_meta[, unit_id := paste(sample_unit_key, project_broad_compartment, sep = "||")]

  units <- cell_meta[non_other, .(
    n_cells = .N,
    expected_library_size = sum(cell_library_size),
    study_from_cells = unique(study)[1],
    cancer_type_from_cells = unique(cancer_type)[1]
  ), by = .(unit_id, sample_unit_key, project_broad_compartment)]
  setorder(units, project_broad_compartment, sample_unit_key)
  units[, unit_index := .I]

  eligible_idx <- which(non_other)
  group_index <- match(cell_meta$unit_id[eligible_idx], units$unit_id)
  if (anyNA(group_index)) stop("Failed to map cells to pseudobulk units")

  log_event("CHECKPOINT", "aggregate sparse counts", n_cells = length(eligible_idx), n_units = nrow(units))
  aggregation_matrix <- sparseMatrix(
    i = eligible_idx,
    j = group_index,
    x = 1,
    dims = c(ncol(counts), nrow(units))
  )
  pb <- counts %*% aggregation_matrix
  rownames(pb) <- genes
  colnames(pb) <- units$unit_id
  rm(aggregation_matrix)
  invisible(gc())

  observed_library <- as.numeric(Matrix::colSums(pb))
  detected_genes <- as.integer(Matrix::colSums(pb > 0))
  units[, `:=`(
    library_size = observed_library,
    detected_genes = detected_genes,
    library_size_match = observed_library == expected_library_size,
    eligible_ge20 = n_cells >= 20L,
    eligible_ge30 = n_cells >= 30L,
    eligible_ge50 = n_cells >= 50L
  )]
  units[, log10_library_size := log10(library_size + 1)]
  units[, log10_n_cells := log10(n_cells + 1)]

  mapping_keep <- setdiff(names(mapping), intersect(names(mapping), c("patient_key")))
  sample_sheet <- merge(units, mapping[, c("sample_unit_key", "patient_key", setdiff(mapping_keep, "sample_unit_key")), with = FALSE],
                        by = "sample_unit_key", all.x = TRUE, sort = FALSE)
  sample_sheet[, cohort_block := paste(study, cancer_type, sep = "||")]
  setorder(sample_sheet, unit_index)
  if (anyNA(sample_sheet$patient_key)) stop("Missing patient mapping after aggregation")
  if (!identical(sample_sheet$unit_id, colnames(pb))) stop("Pseudobulk sample order mismatch")

  gene_meta <- data.table(
    gene = genes,
    is_mitochondrial = grepl("^MT-", genes),
    is_ribosomal = grepl("^RP[SL][0-9]", genes),
    is_definition_exclusion = genes %in% exclusion$gene
  )
  fwrite(gene_meta, file.path(outdir, "pseudobulk_gene_metadata.tsv"), sep = "\t", na = "NA")

  log_event("CHECKPOINT", "save compartment-level SummarizedExperiment objects")
  manifest_rows <- list()
  compartments <- unique(sample_sheet$project_broad_compartment)
  for (compartment in compartments) {
    idx <- which(sample_sheet$project_broad_compartment == compartment)
    comp_counts <- pb[, idx, drop = FALSE]
    comp_meta <- copy(sample_sheet[idx])
    rownames(comp_meta) <- comp_meta$unit_id
    se <- SummarizedExperiment(
      assays = list(counts = comp_counts),
      rowData = DataFrame(gene_meta, row.names = gene_meta$gene),
      colData = DataFrame(comp_meta, row.names = comp_meta$unit_id)
    )
    file_path <- file.path(counts_dir, paste0("pseudobulk_", safe_name(compartment), ".rds"))
    saveRDS(se, file_path, compress = "gzip")
    vals <- comp_counts@x
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
      compartment = compartment,
      file = file_path,
      n_genes = nrow(comp_counts),
      n_units = ncol(comp_counts),
      n_nonzero = length(vals),
      total_library_size = sum(comp_counts),
      nonnegative = all(vals >= 0),
      integer_like = all(abs(vals - round(vals)) < 1e-8),
      file_bytes = as.numeric(file.info(file_path)$size),
      md5 = unname(tools::md5sum(file_path))
    )
    rm(se, comp_counts)
    invisible(gc())
  }
  manifest <- rbindlist(manifest_rows)
  fwrite(manifest, file.path(outdir, "pseudobulk_counts_manifest.tsv"), sep = "\t", na = "NA")
  fwrite(sample_sheet, file.path(outdir, "pseudobulk_sample_sheet.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")

  qc <- sample_sheet[, .(
    unit_id, sample_unit_key, patient_key, study, cancer_type,
    project_broad_compartment, n_cells, expected_library_size,
    observed_library_size = library_size,
    library_size_difference = library_size - expected_library_size,
    library_size_match, detected_genes, eligible_ge20, eligible_ge30, eligible_ge50
  )]
  fwrite(qc, file.path(outdir, "pseudobulk_aggregation_qc.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")

  global_expected <- sum(cell_library[non_other])
  global_observed <- sum(pb)
  counts_pass <- all(qc$library_size_match) && global_expected == global_observed &&
    all(manifest$nonnegative) && all(manifest$integer_like)
  gate <- data.table(
    gate = "G2A_pseudobulk_counts",
    status = if (counts_pass) "PASS" else "FAIL",
    n_cells_aggregated = sum(non_other),
    n_units = nrow(sample_sheet),
    n_compartments = uniqueN(sample_sheet$project_broad_compartment),
    n_genes = nrow(pb),
    global_expected_library = global_expected,
    global_observed_library = global_observed,
    n_unit_mismatches = sum(!qc$library_size_match),
    note = "Technical count aggregation gate only; no differential-expression inference was performed."
  )
  fwrite(gate, file.path(outdir, "g2a_counts_gate_summary.tsv"), sep = "\t", na = "NA")

  run_manifest <- data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version",
             "Matrix_version", "SummarizedExperiment_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string,
              as.character(packageVersion("Matrix")), as.character(packageVersion("SummarizedExperiment")), outdir)
  )
  fwrite(run_manifest, file.path(outdir, "02a_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pseudobulk count aggregation completed",
            runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
            n_units = nrow(sample_sheet), status = gate$status[[1]], output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
