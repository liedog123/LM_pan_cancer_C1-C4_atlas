#!/usr/bin/env Rscript
# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Task: single_cell_states; step: s04
# Inputs (DCF fields): seurat_object,cell_metadata,sample_metadata,celltype_to_broad,compartment_manifest,sample_inclusion_rule,treatment_map,lm12_manifest,lm12_gene_long
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_states --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 04_score_initial_lm12.R --config PATH [--check-config]\n       Rscript 04_score_initial_lm12.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("seurat_object", "cell_metadata", "sample_metadata", "celltype_to_broad", "compartment_manifest", "sample_inclusion_rule", "treatment_map", "lm12_manifest", "lm12_gene_long", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "cell_metadata", "sample_metadata", "celltype_to_broad", "compartment_manifest", "sample_inclusion_rule", "treatment_map", "lm12_manifest", "lm12_gene_long")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),config_path,package_root)
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

# Source-preserved defaults; no study filter is applied.
options(stringsAsFactors = FALSE)
outdir <- out_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
required_pkgs <- c("AUCell", "Matrix")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0L) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

input_rds <- unname(paths["seurat_object"])
input_cell_metadata <- unname(paths["cell_metadata"])
input_sample_metadata <- unname(paths["sample_metadata"])
celltype_to_broad_path <- unname(paths["celltype_to_broad"])
compartment_manifest_path <- unname(paths["compartment_manifest"])
sample_inclusion_rule_path <- unname(paths["sample_inclusion_rule"])
treatment_map_path <- unname(paths["treatment_map"])
lm12_manifest_path <- unname(paths["lm12_manifest"])
lm12_gene_long_path <- unname(paths["lm12_gene_long"])
score_assay <- "RNA"
rank_plot_stats <- FALSE
rank_split_by_blocks <- TRUE
rank_keep_zeroes_as_na <- FALSE
auc_norm <- TRUE
auc_rank_fraction <- 0.05
exclude_compartments <- c("Other")
study_filter_arg <- ""

write_tsv <- function(df, path, compress = FALSE) {
  con <- if (compress) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

safe_char <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  trimws(x)
}

round_or_na <- function(x, digits = 6L) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  round(x, digits = digits)
}

collapse_or_none <- function(x, n_max = 10L) {
  x <- unique(safe_char(x))
  x <- x[x != "" & x != "NA"]
  if (length(x) == 0L) {
    return("none")
  }
  paste(head(sort(x), n_max), collapse = ",")
}

parse_filter_values <- function(x) {
  vals <- safe_char(unlist(strsplit(x, ",")))
  vals[vals != "" & vals != "NA"]
}

`%notin%` <- function(x, table) !(x %in% table)

extract_scoring_matrix <- function(seurat_obj, assay_name) {
  assay_obj <- seurat_obj@assays[[assay_name]]
  if (is.null(assay_obj)) {
    stop("Assay not found: ", assay_name)
  }

  # AUCell 只依赖单细胞内的表达排序，这里优先使用 counts，避免依赖不稳定的 data 层。
  counts_mat <- assay_obj@counts
  if (!is.null(counts_mat) && nrow(counts_mat) > 0L && ncol(counts_mat) > 0L) {
    return(list(matrix = counts_mat, slot_used = "counts"))
  }

  data_mat <- assay_obj@data
  if (!is.null(data_mat) && nrow(data_mat) > 0L && ncol(data_mat) > 0L) {
    return(list(matrix = data_mat, slot_used = "data"))
  }

  stop("No usable expression matrix found in assay: ", assay_name)
}

required_files <- c(
  input_rds,
  input_cell_metadata,
  input_sample_metadata,
  celltype_to_broad_path,
  compartment_manifest_path,
  sample_inclusion_rule_path,
  treatment_map_path,
  lm12_manifest_path,
  lm12_gene_long_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

message("Loading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
expr_info <- extract_scoring_matrix(obj, score_assay)
expr_mat <- expr_info$matrix
score_slot_used <- expr_info$slot_used
obj_cell_ids <- colnames(obj)

cell_metadata <- read.delim(gzfile(input_cell_metadata), check.names = FALSE)
sample_metadata <- read.delim(input_sample_metadata, check.names = FALSE)
celltype_to_broad <- read.delim(celltype_to_broad_path, check.names = FALSE)
compartment_manifest <- read.delim(compartment_manifest_path, check.names = FALSE)
sample_inclusion_rule <- read.delim(sample_inclusion_rule_path, check.names = FALSE)
treatment_map <- read.delim(treatment_map_path, check.names = FALSE)
lm12_manifest <- read.delim(lm12_manifest_path, check.names = FALSE)
lm12_gene_long <- read.delim(lm12_gene_long_path, check.names = FALSE)

required_cell_cols <- c(
  "cell_id", "study", "cancer_type", "donor_id", "Sample", "donor_id_pre_post",
  "sample_unit_key", "pre_post", "outcome", "Combined_outcome", "treatment_raw",
  "treatment_class_broad", "cell_types", "cell_types_update",
  "Cell_type_broad_original", "project_broad_compartment", "nCount_RNA", "nFeature_RNA"
)
missing_cell_cols <- required_cell_cols[!required_cell_cols %in% colnames(cell_metadata)]
if (length(missing_cell_cols) > 0L) {
  stop("analysis_ready_cell_metadata is missing required columns: ", paste(missing_cell_cols, collapse = ", "))
}

if (!setequal(cell_metadata$cell_id, obj_cell_ids)) {
  stop("Cell IDs differ between frozen metadata and Seurat object.")
}
cell_metadata <- cell_metadata[match(obj_cell_ids, cell_metadata$cell_id), , drop = FALSE]
if (any(is.na(cell_metadata$cell_id))) {
  stop("Failed to align frozen cell metadata to Seurat object.")
}

main_compartment_manifest <- compartment_manifest[compartment_manifest$analysis_mode == "main", , drop = FALSE]
if (nrow(main_compartment_manifest) == 0L) {
  stop("No rows with analysis_mode == 'main' found in compartment_sample_set_manifest.")
}
broad_levels <- safe_char(main_compartment_manifest$project_broad_compartment)
broad_levels <- broad_levels[!duplicated(broad_levels)]
scoring_broad_levels <- broad_levels[broad_levels %notin% exclude_compartments]
compartment_role_lookup <- setNames(
  safe_char(main_compartment_manifest$compartment_role),
  safe_char(main_compartment_manifest$project_broad_compartment)
)

feature_metadata <- lm12_manifest[order(lm12_manifest$display_order), c(
  "display_order", "feature_id", "short_label", "pathway", "direction_label",
  "feature_class", "functional_group", "stable_gene_count"
)]
feature_ids <- safe_char(feature_metadata$feature_id)
gene_sets <- lapply(feature_ids, function(fid) {
  genes <- unique(safe_char(lm12_gene_long$gene[lm12_gene_long$feature_id == fid]))
  genes <- genes[genes != "" & genes != "NA"]
  sort(genes)
})
names(gene_sets) <- feature_ids

expr_genes <- rownames(expr_mat)
gene_sets_present <- lapply(gene_sets, function(genes) intersect(genes, expr_genes))
names(gene_sets_present) <- feature_ids

feature_gene_coverage_overall <- data.frame(
  feature_id = feature_ids,
  short_label = feature_metadata$short_label,
  pathway = feature_metadata$pathway,
  genes_total = vapply(gene_sets, length, integer(1)),
  genes_present_in_assay = vapply(gene_sets_present, length, integer(1)),
  coverage_in_assay = round_or_na(vapply(seq_along(gene_sets_present), function(i) {
    denom <- length(gene_sets[[i]])
    if (denom == 0L) {
      return(NA_real_)
    }
    length(gene_sets_present[[i]]) / denom
  }, numeric(1))),
  stringsAsFactors = FALSE
)

if (any(feature_gene_coverage_overall$genes_present_in_assay == 0L)) {
  bad_features <- feature_gene_coverage_overall$feature_id[feature_gene_coverage_overall$genes_present_in_assay == 0L]
  stop("Some LM12 features have zero gene coverage in the selected assay: ", paste(bad_features, collapse = ", "))
}

study_filter <- parse_filter_values(study_filter_arg)
if (length(study_filter) > 0L) {
  cell_metadata <- cell_metadata[safe_char(cell_metadata$study) %in% study_filter, , drop = FALSE]
  if (nrow(cell_metadata) == 0L) {
    stop("No cells remain after applying --study_filter.")
  }
}

cell_metadata <- cell_metadata[cell_metadata$project_broad_compartment %in% broad_levels, , drop = FALSE]
cell_metadata <- cell_metadata[!(cell_metadata$project_broad_compartment %in% exclude_compartments), , drop = FALSE]
if (nrow(cell_metadata) == 0L) {
  stop("No non-excluded cells available for Step 02 scoring.")
}

cell_metadata$project_broad_compartment <- factor(
  safe_char(cell_metadata$project_broad_compartment),
  levels = scoring_broad_levels
)
cell_metadata$study <- safe_char(cell_metadata$study)
cell_metadata$study_compartment_key <- paste(cell_metadata$study, cell_metadata$project_broad_compartment, sep = "||")
cell_metadata$compartment_role <- unname(compartment_role_lookup[as.character(cell_metadata$project_broad_compartment)])

study_compartment_keys <- unique(cell_metadata$study_compartment_key)
study_compartment_keys <- study_compartment_keys[order(
  sub("\\|\\|.*$", "", study_compartment_keys),
  match(sub("^.*\\|\\|", "", study_compartment_keys), levels(cell_metadata$project_broad_compartment))
)]

all_signature_genes <- sort(unique(unlist(gene_sets_present)))
score_rows <- list()
manifest_rows <- list()
coverage_rows <- list()
summary_rows <- list()

for (unit_key in study_compartment_keys) {
  unit_meta <- cell_metadata[cell_metadata$study_compartment_key == unit_key, , drop = FALSE]
  unit_meta <- unit_meta[order(unit_meta$sample_unit_key, unit_meta$cell_id), , drop = FALSE]
  unit_cells <- unit_meta$cell_id
  unit_expr <- expr_mat[, unit_cells, drop = FALSE]

  unit_study <- unit_meta$study[1]
  unit_compartment <- as.character(unit_meta$project_broad_compartment[1])
  unit_role <- safe_char(unit_meta$compartment_role[1])
  unit_cancer_types <- collapse_or_none(unique(unit_meta$cancer_type))
  n_unit_cells <- nrow(unit_meta)
  n_unit_samples <- length(unique(unit_meta$sample_unit_key))
  n_unit_donors <- length(unique(unit_meta$donor_id))
  n_unit_cell_types_update <- length(unique(unit_meta$cell_types_update))
  auc_max_rank <- ceiling(auc_rank_fraction * nrow(unit_expr))

  message("Running AUCell for ", unit_key, " (cells=", n_unit_cells, ").")
  rankings <- AUCell::AUCell_buildRankings(
    exprMat = unit_expr,
    featureType = "genes",
    plotStats = rank_plot_stats,
    splitByBlocks = rank_split_by_blocks,
    keepZeroesAsNA = rank_keep_zeroes_as_na,
    verbose = FALSE
  )
  auc_obj <- AUCell::AUCell_calcAUC(
    geneSets = gene_sets_present,
    rankings = rankings,
    nCores = 1,
    normAUC = auc_norm,
    aucMaxRank = auc_max_rank,
    verbose = FALSE
  )
  auc_mat <- AUCell::getAUC(auc_obj)
  auc_mat <- t(as.matrix(auc_mat))
  auc_df <- as.data.frame(auc_mat, check.names = FALSE, stringsAsFactors = FALSE)
  auc_df$cell_id <- rownames(auc_mat)
  auc_df <- auc_df[, c("cell_id", feature_ids), drop = FALSE]
  colnames(auc_df) <- c("cell_id", paste0("initial__", feature_ids))

  unit_score_table <- merge(
    unit_meta,
    auc_df,
    by = "cell_id",
    all.x = TRUE,
    sort = FALSE
  )
  unit_score_table <- unit_score_table[match(unit_cells, unit_score_table$cell_id), , drop = FALSE]
  unit_score_table$score_method_02 <- "AUCell"
  unit_score_table$score_assay_02 <- score_assay
  unit_score_table$score_slot_02 <- score_slot_used
  unit_score_table$score_status_02 <- "processed"
  score_rows[[length(score_rows) + 1L]] <- unit_score_table

  detected_counts <- Matrix::rowSums(unit_expr[all_signature_genes, , drop = FALSE] > 0)
  detected_fraction <- detected_counts / n_unit_cells

  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    study = unit_study,
    cancer_type_summary = unit_cancer_types,
    project_broad_compartment = unit_compartment,
    compartment_role = unit_role,
    study_compartment_key = unit_key,
    n_cells = n_unit_cells,
    n_sample_units = n_unit_samples,
    n_donors = n_unit_donors,
    n_unique_cell_types_update = n_unit_cell_types_update,
    n_cells_scored = nrow(auc_df),
    score_method = "AUCell",
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    auc_max_rank = auc_max_rank,
    stringsAsFactors = FALSE
  )

  for (fid in feature_ids) {
    genes_present <- gene_sets_present[[fid]]
    gene_detected_fraction <- detected_fraction[genes_present]
    feature_scores <- unit_score_table[[paste0("initial__", fid)]]

    coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(
      study = unit_study,
      cancer_type_summary = unit_cancer_types,
      project_broad_compartment = unit_compartment,
      compartment_role = unit_role,
      study_compartment_key = unit_key,
      feature_id = fid,
      short_label = feature_metadata$short_label[feature_metadata$feature_id == fid][1],
      pathway = feature_metadata$pathway[feature_metadata$feature_id == fid][1],
      genes_total = length(gene_sets[[fid]]),
      genes_present_in_assay = length(genes_present),
      genes_detected_in_unit = sum(detected_counts[genes_present] > 0),
      median_gene_detected_fraction = round_or_na(stats::median(gene_detected_fraction)),
      min_gene_detected_fraction = round_or_na(min(gene_detected_fraction)),
      max_gene_detected_fraction = round_or_na(max(gene_detected_fraction)),
      stringsAsFactors = FALSE
    )

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      study = unit_study,
      cancer_type_summary = unit_cancer_types,
      project_broad_compartment = unit_compartment,
      compartment_role = unit_role,
      study_compartment_key = unit_key,
      feature_id = fid,
      short_label = feature_metadata$short_label[feature_metadata$feature_id == fid][1],
      pathway = feature_metadata$pathway[feature_metadata$feature_id == fid][1],
      genes_total = length(gene_sets[[fid]]),
      genes_present_in_assay = length(genes_present),
      n_cells = length(feature_scores),
      mean_initial_score = round_or_na(mean(feature_scores)),
      median_initial_score = round_or_na(stats::median(feature_scores)),
      sd_initial_score = round_or_na(stats::sd(feature_scores)),
      q05 = round_or_na(stats::quantile(feature_scores, probs = 0.05, names = FALSE)),
      q95 = round_or_na(stats::quantile(feature_scores, probs = 0.95, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }
}

score_table <- do.call(rbind, score_rows)
manifest_df <- do.call(rbind, manifest_rows)
coverage_df <- do.call(rbind, coverage_rows)
summary_df <- do.call(rbind, summary_rows)

score_table <- score_table[order(
  score_table$study,
  match(as.character(score_table$project_broad_compartment), scoring_broad_levels),
  score_table$sample_unit_key,
  score_table$cell_id
), , drop = FALSE]
manifest_df <- manifest_df[order(
  manifest_df$study,
  match(manifest_df$project_broad_compartment, scoring_broad_levels)
), , drop = FALSE]
coverage_df <- coverage_df[order(
  coverage_df$study,
  match(coverage_df$project_broad_compartment, scoring_broad_levels),
  match(coverage_df$feature_id, feature_ids)
), , drop = FALSE]
summary_df <- summary_df[order(
  summary_df$study,
  match(summary_df$project_broad_compartment, scoring_broad_levels),
  match(summary_df$feature_id, feature_ids)
), , drop = FALSE]

write_tsv(manifest_df, file.path(outdir, "study_compartment_manifest.tsv"))
write_tsv(coverage_df, file.path(outdir, "study_compartment_feature_gene_coverage.tsv"))
write_tsv(summary_df, file.path(outdir, "study_compartment_feature_initial_score_summary.tsv"))
write_tsv(score_table, file.path(outdir, "cell_level_lm12_initial_scores.tsv.gz"), compress = TRUE)

bundle <- list(
  parameters = list(
    score_method = "AUCell",
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    auc_rank_fraction = auc_rank_fraction,
    exclude_compartments = exclude_compartments,
    study_filter = study_filter
  ),
  feature_metadata = feature_metadata,
  feature_gene_coverage_overall = feature_gene_coverage_overall,
  gene_sets = gene_sets,
  gene_sets_present = gene_sets_present,
  study_compartment_manifest = manifest_df,
  study_compartment_feature_gene_coverage = coverage_df,
  study_compartment_feature_initial_score_summary = summary_df,
  cell_level_score_table = score_table,
  cell_metadata_columns = colnames(cell_metadata),
  sample_metadata = sample_metadata,
  celltype_to_broad = celltype_to_broad,
  compartment_manifest = compartment_manifest,
  sample_inclusion_rule = sample_inclusion_rule,
  treatment_map = treatment_map
)
saveRDS(bundle, file.path(outdir, "initial_scoring_bundle.rds"))

summary_lines <- c(
  "02_studywise_cell_level_lm12_initial_scoring completed.",
  paste0("Input Seurat object: ", input_rds),
  paste0("Input cell metadata: ", input_cell_metadata),
  paste0("Input sample metadata: ", input_sample_metadata),
  paste0("Score assay: ", score_assay),
  paste0("Score slot used: ", score_slot_used),
  paste0("Score method: AUCell"),
  paste0("AUC rank fraction: ", auc_rank_fraction),
  paste0("Studies processed: ", length(unique(manifest_df$study))),
  paste0("Study x compartment units processed: ", nrow(manifest_df)),
  paste0("Cells scored: ", nrow(score_table)),
  paste0("Output cell-level table: ", file.path(outdir, "cell_level_lm12_initial_scores.tsv.gz")),
  paste0("Output bundle: ", file.path(outdir, "initial_scoring_bundle.rds"))
)
writeLines(summary_lines, con = file.path(outdir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
