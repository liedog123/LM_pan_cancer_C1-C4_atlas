# Purpose: single cell states analysis stages.
# Configuration: config/single_cell.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/04_score_initial_lm12.R
stage_single_cell_04_score_initial_lm12 <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "cell_metadata", "sample_metadata", "celltype_to_broad", "compartment_manifest", "sample_inclusion_rule", "treatment_map", "lm12_manifest", "lm12_gene_long")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`write_tsv` <- bind_stage_helper(shared_write_tsv_006, environment())
`safe_char` <- bind_stage_helper(shared_safe_char_042, environment())
`round_or_na` <- bind_stage_helper(shared_round_or_na_045, environment())

collapse_or_none <- function(x, n_max = 10L) {
  x <- unique(safe_char(x))
  x <- x[x != "" & x != "NA"]
  if (length(x) == 0L) {
    return("none")
  }
  paste(head(sort(x), n_max), collapse = ",")
}
`parse_filter_values` <- bind_stage_helper(shared_parse_filter_values_040, environment())
`%notin%` <- bind_stage_helper(shared__notin__055, environment())
`extract_scoring_matrix` <- bind_stage_helper(shared_extract_scoring_matrix_035, environment())

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
}

# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/05_refine_and_rescore_lm12.R
stage_single_cell_05_refine_and_rescore_lm12 <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "initial_bundle")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
input_bundle <- unname(paths["initial_bundle"])
score_assay <- "RNA"
rank_plot_stats <- FALSE
rank_split_by_blocks <- TRUE
rank_keep_zeroes_as_na <- FALSE
auc_norm <- TRUE
auc_rank_fraction <- 0.05
negative_correlation_threshold <- 0
min_cells_for_correlation <- 25L
low_support_gene_threshold <- 5L
min_retained_genes_for_direct_refined_scoring <- 1L
allow_original_fallback <- TRUE
study_filter_arg <- ""
study_compartment_filter_arg <- ""
`write_tsv` <- bind_stage_helper(shared_write_tsv_006, environment())
`safe_char` <- bind_stage_helper(shared_safe_char_042, environment())
`round_or_na` <- bind_stage_helper(shared_round_or_na_045, environment())

collapse_or_none <- function(x, n_max = 20L) {
  x <- unique(safe_char(x))
  x <- x[x != "" & x != "NA"]
  if (length(x) == 0L) {
    return("none")
  }
  paste(head(sort(x), n_max), collapse = ",")
}
`parse_filter_values` <- bind_stage_helper(shared_parse_filter_values_040, environment())
`%notin%` <- bind_stage_helper(shared__notin__055, environment())
`extract_scoring_matrix` <- bind_stage_helper(shared_extract_scoring_matrix_035, environment())

calc_spearman_safe <- function(x, y, min_cells = 25L) {
  keep <- is.finite(x) & is.finite(y)
  n_keep <- sum(keep)
  if (n_keep < min_cells) {
    return(c(rho = NA_real_, p_value = NA_real_, n = n_keep))
  }
  if (length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(c(rho = NA_real_, p_value = NA_real_, n = n_keep))
  }

  test_res <- suppressWarnings(
    stats::cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
  )
  c(rho = unname(test_res$estimate), p_value = test_res$p.value, n = n_keep)
}

score_gene_sets_auc <- function(rankings, gene_sets, auc_max_rank, auc_norm = TRUE) {
  if (length(gene_sets) == 0L) {
    return(matrix(numeric(0), nrow = ncol(rankings), ncol = 0L))
  }
  auc_obj <- AUCell::AUCell_calcAUC(
    geneSets = gene_sets,
    rankings = rankings,
    nCores = 1,
    normAUC = auc_norm,
    aucMaxRank = auc_max_rank,
    verbose = FALSE
  )
  auc_mat <- AUCell::getAUC(auc_obj)
  auc_mat <- t(as.matrix(auc_mat))
  auc_mat
}

required_files <- c(input_rds, input_bundle)
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

bundle <- readRDS(input_bundle)
required_bundle_names <- c(
  "parameters",
  "feature_metadata",
  "gene_sets_present",
  "study_compartment_manifest",
  "cell_level_score_table",
  "sample_metadata",
  "celltype_to_broad",
  "compartment_manifest",
  "sample_inclusion_rule",
  "treatment_map"
)
missing_bundle_names <- required_bundle_names[required_bundle_names %notin% names(bundle)]
if (length(missing_bundle_names) > 0L) {
  stop("Initial scoring bundle is missing components: ", paste(missing_bundle_names, collapse = ", "))
}

feature_metadata <- bundle$feature_metadata
gene_sets_present <- bundle$gene_sets_present
study_compartment_manifest <- bundle$study_compartment_manifest
cell_score_table <- bundle$cell_level_score_table
sample_metadata <- bundle$sample_metadata
celltype_to_broad <- bundle$celltype_to_broad
compartment_manifest <- bundle$compartment_manifest
sample_inclusion_rule <- bundle$sample_inclusion_rule
treatment_map <- bundle$treatment_map

feature_ids <- safe_char(feature_metadata$feature_id)
initial_score_cols <- paste0("initial__", feature_ids)
missing_initial_score_cols <- initial_score_cols[initial_score_cols %notin% colnames(cell_score_table)]
if (length(missing_initial_score_cols) > 0L) {
  stop("Initial scoring table is missing columns: ", paste(missing_initial_score_cols, collapse = ", "))
}

required_cell_cols <- c(
  "cell_id", "study", "cancer_type", "donor_id", "Sample", "donor_id_pre_post",
  "sample_unit_key", "pre_post", "outcome", "Combined_outcome", "treatment_raw",
  "treatment_class_broad", "cell_types", "cell_types_update", "Cell_type_broad_original",
  "project_broad_compartment", "study_compartment_key", "compartment_role"
)
missing_cell_cols <- required_cell_cols[required_cell_cols %notin% colnames(cell_score_table)]
if (length(missing_cell_cols) > 0L) {
  stop("Initial scoring table is missing metadata columns: ", paste(missing_cell_cols, collapse = ", "))
}

if (!all(cell_score_table$cell_id %in% obj_cell_ids)) {
  stop("Some cell IDs in initial scoring bundle are not present in the Seurat object.")
}

study_filter <- parse_filter_values(study_filter_arg)
study_compartment_filter <- parse_filter_values(study_compartment_filter_arg)

if (length(study_filter) > 0L) {
  keep_units <- safe_char(study_compartment_manifest$study) %in% study_filter
  study_compartment_manifest <- study_compartment_manifest[keep_units, , drop = FALSE]
}
if (length(study_compartment_filter) > 0L) {
  keep_units <- safe_char(study_compartment_manifest$study_compartment_key) %in% study_compartment_filter
  study_compartment_manifest <- study_compartment_manifest[keep_units, , drop = FALSE]
}
if (nrow(study_compartment_manifest) == 0L) {
  stop("No study-compartment units remain after filtering.")
}

unit_keys_to_keep <- safe_char(study_compartment_manifest$study_compartment_key)
cell_score_table <- cell_score_table[safe_char(cell_score_table$study_compartment_key) %in% unit_keys_to_keep, , drop = FALSE]
if (nrow(cell_score_table) == 0L) {
  stop("No cell-level rows remain after filtering the initial scoring table.")
}

compartment_levels <- safe_char(study_compartment_manifest$project_broad_compartment)
compartment_levels <- compartment_levels[!duplicated(compartment_levels)]

correlation_rows <- list()
feature_summary_rows <- list()
refined_gene_set_rows <- list()
refined_score_summary_rows <- list()
manifest_rows <- list()
cell_rows <- list()

retained_gene_sets_by_unit <- list()
scoring_gene_sets_by_unit <- list()

for (idx in seq_len(nrow(study_compartment_manifest))) {
  unit_row <- study_compartment_manifest[idx, , drop = FALSE]
  unit_key <- safe_char(unit_row$study_compartment_key[1])
  unit_study <- safe_char(unit_row$study[1])
  unit_cancer_types <- safe_char(unit_row$cancer_type_summary[1])
  unit_compartment <- safe_char(unit_row$project_broad_compartment[1])
  unit_role <- safe_char(unit_row$compartment_role[1])

  unit_meta <- cell_score_table[safe_char(cell_score_table$study_compartment_key) == unit_key, , drop = FALSE]
  unit_meta <- unit_meta[order(unit_meta$sample_unit_key, unit_meta$cell_id), , drop = FALSE]
  unit_cells <- unit_meta$cell_id
  n_unit_cells <- nrow(unit_meta)
  if (n_unit_cells == 0L) {
    stop("No cells found for study_compartment_key: ", unit_key)
  }

  unit_expr <- expr_mat[, unit_cells, drop = FALSE]
  auc_max_rank <- as.integer(unit_row$auc_max_rank[1])
  if (!is.finite(auc_max_rank) || auc_max_rank <= 0L) {
    auc_max_rank <- ceiling(auc_rank_fraction * nrow(unit_expr))
  }
  auc_max_rank <- min(auc_max_rank, nrow(unit_expr))

  message("Refining AUCell gene sets for ", unit_key, " (cells=", n_unit_cells, ").")
  rankings <- AUCell::AUCell_buildRankings(
    exprMat = unit_expr,
    featureType = "genes",
    plotStats = rank_plot_stats,
    splitByBlocks = rank_split_by_blocks,
    keepZeroesAsNA = rank_keep_zeroes_as_na,
    verbose = FALSE
  )

  original_gene_sets_unit <- lapply(feature_ids, function(fid) {
    genes_now <- safe_char(gene_sets_present[[fid]])
    genes_now[genes_now != "" & genes_now != "NA"]
  })
  names(original_gene_sets_unit) <- feature_ids

  loo_signature_list <- list()
  for (fid in feature_ids) {
    genes_use <- original_gene_sets_unit[[fid]]
    if (length(genes_use) <= 1L) {
      next
    }
    for (gene in genes_use) {
      loo_signature_list[[paste(fid, gene, sep = "__loo__")]] <- setdiff(genes_use, gene)
    }
  }

  if (length(loo_signature_list) > 0L) {
    loo_auc_mat <- score_gene_sets_auc(
      rankings = rankings,
      gene_sets = loo_signature_list,
      auc_max_rank = auc_max_rank,
      auc_norm = auc_norm
    )
    rownames(loo_auc_mat) <- unit_cells
  } else {
    loo_auc_mat <- matrix(numeric(0), nrow = n_unit_cells, ncol = 0L, dimnames = list(unit_cells, character(0)))
  }

  refined_gene_sets_unit <- list()
  scoring_gene_sets_unit <- list()
  feature_result_map <- list()
  low_support_features <- character(0)
  unavailable_features <- character(0)
  fallback_features <- character(0)
  unscored_features <- character(0)

  for (fid in feature_ids) {
    feature_short_label <- feature_metadata$short_label[feature_metadata$feature_id == fid][1]
    feature_pathway <- feature_metadata$pathway[feature_metadata$feature_id == fid][1]
    genes_use <- original_gene_sets_unit[[fid]]

    if (length(genes_use) == 0L) {
      stop("Feature ", fid, " has zero assay-present genes in unit ", unit_key, ".")
    }

    if (length(genes_use) == 1L) {
      gene <- genes_use[1]
      gene_detected_cells <- as.integer(sum(unit_expr[gene, , drop = FALSE] > 0))
      corr_df <- data.frame(
        study = unit_study,
        cancer_type_summary = unit_cancer_types,
        project_broad_compartment = unit_compartment,
        compartment_role = unit_role,
        study_compartment_key = unit_key,
        feature_id = fid,
        short_label = feature_short_label,
        pathway = feature_pathway,
        gene = gene,
        leave_one_out_rho = NA_real_,
        leave_one_out_p_value = NA_real_,
        n_cells_used = n_unit_cells,
        gene_detected_cells = gene_detected_cells,
        gene_detected_fraction = round_or_na(gene_detected_cells / n_unit_cells),
        keep_gene_refined = gene_detected_cells > 0L,
        gene_decision_reason = if (gene_detected_cells > 0L) {
          "single_gene_feature_retained"
        } else {
          "undetected_in_unit"
        },
        stringsAsFactors = FALSE
      )
    } else {
      corr_df <- do.call(
        rbind,
        lapply(genes_use, function(gene) {
          loo_signature_id <- paste(fid, gene, sep = "__loo__")
          gene_vec <- as.numeric(unit_expr[gene, , drop = TRUE])
          loo_vec <- as.numeric(loo_auc_mat[, loo_signature_id])
          cor_stats <- calc_spearman_safe(gene_vec, loo_vec, min_cells = min_cells_for_correlation)
          gene_detected_cells <- as.integer(sum(gene_vec > 0))
          data.frame(
            study = unit_study,
            cancer_type_summary = unit_cancer_types,
            project_broad_compartment = unit_compartment,
            compartment_role = unit_role,
            study_compartment_key = unit_key,
            feature_id = fid,
            short_label = feature_short_label,
            pathway = feature_pathway,
            gene = gene,
            leave_one_out_rho = as.numeric(cor_stats["rho"]),
            leave_one_out_p_value = as.numeric(cor_stats["p_value"]),
            n_cells_used = as.integer(cor_stats["n"]),
            gene_detected_cells = gene_detected_cells,
            gene_detected_fraction = round_or_na(gene_detected_cells / n_unit_cells),
            stringsAsFactors = FALSE
          )
        })
      )

      corr_df$keep_gene_refined <- TRUE
      corr_df$gene_decision_reason <- "retained"
      corr_df$keep_gene_refined[is.na(corr_df$leave_one_out_rho) & corr_df$gene_detected_cells == 0L] <- FALSE
      corr_df$gene_decision_reason[is.na(corr_df$leave_one_out_rho) & corr_df$gene_detected_cells == 0L] <- "undetected_in_unit"
      corr_df$keep_gene_refined[!is.na(corr_df$leave_one_out_rho) & corr_df$leave_one_out_rho < negative_correlation_threshold] <- FALSE
      corr_df$gene_decision_reason[!is.na(corr_df$leave_one_out_rho) & corr_df$leave_one_out_rho < negative_correlation_threshold] <- "negative_leave_one_out_correlation"
      corr_df$gene_decision_reason[is.na(corr_df$leave_one_out_rho) & corr_df$gene_detected_cells > 0L] <- "retained_na_correlation"
    }

    retained_genes <- unique(corr_df$gene[corr_df$keep_gene_refined])
    feature_refine_status <- if (length(retained_genes) == 0L) {
      "unavailable"
    } else if (length(retained_genes) < low_support_gene_threshold) {
      "low_support"
    } else {
      "usable"
    }

    direct_refined_usable <- length(retained_genes) >= min_retained_genes_for_direct_refined_scoring
    if (direct_refined_usable) {
      scoring_genes <- retained_genes
      scoring_gene_set_source <- "refined_retained"
    } else if (allow_original_fallback && length(genes_use) > 0L) {
      scoring_genes <- genes_use
      scoring_gene_set_source <- if (length(retained_genes) == 0L) {
        "original_fallback_zero_retained"
      } else {
        "original_fallback_below_min_retained"
      }
    } else {
      scoring_genes <- character(0)
      scoring_gene_set_source <- "unscored"
    }

    if (feature_refine_status == "low_support") {
      low_support_features <- c(low_support_features, fid)
    }
    if (feature_refine_status == "unavailable") {
      unavailable_features <- c(unavailable_features, fid)
    }
    if (grepl("^original_fallback", scoring_gene_set_source)) {
      fallback_features <- c(fallback_features, fid)
    }
    if (scoring_gene_set_source == "unscored") {
      unscored_features <- c(unscored_features, fid)
    }

    refined_gene_sets_unit[[fid]] <- retained_genes
    scoring_gene_sets_unit[[fid]] <- scoring_genes
    feature_result_map[[fid]] <- list(
      corr_df = corr_df,
      retained_genes = retained_genes,
      scoring_genes = scoring_genes,
      feature_refine_status = feature_refine_status,
      scoring_gene_set_source = scoring_gene_set_source,
      n_original_genes = length(genes_use),
      n_retained_genes = length(retained_genes),
      n_scoring_genes = length(scoring_genes)
    )
  }

  retained_gene_sets_by_unit[[unit_key]] <- refined_gene_sets_unit
  scoring_gene_sets_by_unit[[unit_key]] <- scoring_gene_sets_unit

  refined_auc_mat <- matrix(
    NA_real_,
    nrow = n_unit_cells,
    ncol = length(feature_ids),
    dimnames = list(unit_cells, feature_ids)
  )
  scoreable_gene_sets <- scoring_gene_sets_unit[lengths(scoring_gene_sets_unit) > 0L]
  if (length(scoreable_gene_sets) > 0L) {
    scoreable_auc_mat <- score_gene_sets_auc(
      rankings = rankings,
      gene_sets = scoreable_gene_sets,
      auc_max_rank = auc_max_rank,
      auc_norm = auc_norm
    )
    rownames(scoreable_auc_mat) <- unit_cells
    refined_auc_mat[unit_cells, colnames(scoreable_auc_mat)] <- scoreable_auc_mat
  }

  unit_cell_table <- unit_meta
  for (fid in feature_ids) {
    feature_res <- feature_result_map[[fid]]
    corr_key <- paste(unit_key, fid, sep = "__")
    correlation_rows[[corr_key]] <- feature_res$corr_df

    init_vec <- as.numeric(unit_meta[[paste0("initial__", fid)]])
    ref_vec <- as.numeric(refined_auc_mat[, fid])
    init_ref_cor <- if (sum(is.finite(init_vec) & is.finite(ref_vec)) >= min_cells_for_correlation) {
      suppressWarnings(stats::cor(init_vec, ref_vec, method = "spearman", use = "complete.obs"))
    } else {
      NA_real_
    }

    feature_summary_rows[[corr_key]] <- data.frame(
      study = unit_study,
      cancer_type_summary = unit_cancer_types,
      project_broad_compartment = unit_compartment,
      compartment_role = unit_role,
      study_compartment_key = unit_key,
      feature_id = fid,
      short_label = feature_metadata$short_label[feature_metadata$feature_id == fid][1],
      pathway = feature_metadata$pathway[feature_metadata$feature_id == fid][1],
      n_original_genes = feature_res$n_original_genes,
      n_tested_genes = sum(!is.na(feature_res$corr_df$leave_one_out_rho)),
      n_negative_cor_genes = sum(feature_res$corr_df$leave_one_out_rho < negative_correlation_threshold, na.rm = TRUE),
      n_undetected_genes = sum(feature_res$corr_df$gene_detected_cells == 0L),
      n_retained_genes = feature_res$n_retained_genes,
      n_scoring_genes = feature_res$n_scoring_genes,
      feature_refine_status = feature_res$feature_refine_status,
      feature_unavailable = feature_res$feature_refine_status == "unavailable",
      low_support_after_refine = feature_res$feature_refine_status == "low_support",
      scoring_gene_set_source = feature_res$scoring_gene_set_source,
      initial_refined_score_rho = round_or_na(init_ref_cor),
      stringsAsFactors = FALSE
    )

    scoring_genes_now <- feature_res$scoring_genes
    if (length(scoring_genes_now) > 0L) {
      refined_gene_set_rows[[corr_key]] <- data.frame(
        study = unit_study,
        cancer_type_summary = unit_cancer_types,
        project_broad_compartment = unit_compartment,
        compartment_role = unit_role,
        study_compartment_key = unit_key,
        feature_id = fid,
        short_label = feature_metadata$short_label[feature_metadata$feature_id == fid][1],
        pathway = feature_metadata$pathway[feature_metadata$feature_id == fid][1],
        feature_refine_status = feature_res$feature_refine_status,
        scoring_gene_set_source = feature_res$scoring_gene_set_source,
        gene = scoring_genes_now,
        gene_membership_type = ifelse(scoring_genes_now %in% feature_res$retained_genes, "retained_refined", "original_fallback_only"),
        stringsAsFactors = FALSE
      )
    }

    ref_finite <- ref_vec[is.finite(ref_vec)]
    if (length(ref_finite) == 0L) {
      score_stats <- c(mean = NA_real_, median = NA_real_, sd = NA_real_, q05 = NA_real_, q95 = NA_real_)
    } else {
      score_stats <- c(
        mean = mean(ref_finite),
        median = stats::median(ref_finite),
        sd = stats::sd(ref_finite),
        q05 = stats::quantile(ref_finite, probs = 0.05, names = FALSE),
        q95 = stats::quantile(ref_finite, probs = 0.95, names = FALSE)
      )
    }
    refined_score_summary_rows[[corr_key]] <- data.frame(
      study = unit_study,
      cancer_type_summary = unit_cancer_types,
      project_broad_compartment = unit_compartment,
      compartment_role = unit_role,
      study_compartment_key = unit_key,
      feature_id = fid,
      short_label = feature_metadata$short_label[feature_metadata$feature_id == fid][1],
      pathway = feature_metadata$pathway[feature_metadata$feature_id == fid][1],
      feature_refine_status = feature_res$feature_refine_status,
      scoring_gene_set_source = feature_res$scoring_gene_set_source,
      n_cells = length(ref_vec),
      n_cells_scored = length(ref_finite),
      mean_refined_score = round_or_na(score_stats["mean"]),
      median_refined_score = round_or_na(score_stats["median"]),
      sd_refined_score = round_or_na(score_stats["sd"]),
      q05_refined_score = round_or_na(score_stats["q05"]),
      q95_refined_score = round_or_na(score_stats["q95"]),
      stringsAsFactors = FALSE
    )

    unit_cell_table[[paste0("refined__", fid)]] <- ref_vec
  }

  unit_cell_table$score_method_03 <- "AUCell"
  unit_cell_table$score_assay_03 <- score_assay
  unit_cell_table$score_slot_03 <- score_slot_used
  unit_cell_table$refine_status_03 <- "processed"
  cell_rows[[unit_key]] <- unit_cell_table

  manifest_rows[[unit_key]] <- data.frame(
    study = unit_study,
    cancer_type_summary = unit_cancer_types,
    project_broad_compartment = unit_compartment,
    compartment_role = unit_role,
    study_compartment_key = unit_key,
    n_cells = n_unit_cells,
    n_sample_units = length(unique(unit_meta$sample_unit_key)),
    n_donors = length(unique(unit_meta$donor_id)),
    n_unique_cell_types_update = length(unique(unit_meta$cell_types_update)),
    score_method = "AUCell",
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    auc_max_rank = auc_max_rank,
    negative_correlation_threshold = negative_correlation_threshold,
    min_cells_for_correlation = min_cells_for_correlation,
    low_support_gene_threshold = low_support_gene_threshold,
    min_retained_genes_for_direct_refined_scoring = min_retained_genes_for_direct_refined_scoring,
    allow_original_fallback = allow_original_fallback,
    n_features_usable = sum(vapply(feature_result_map, function(x) x$feature_refine_status == "usable", logical(1))),
    n_features_low_support = sum(vapply(feature_result_map, function(x) x$feature_refine_status == "low_support", logical(1))),
    n_features_unavailable = sum(vapply(feature_result_map, function(x) x$feature_refine_status == "unavailable", logical(1))),
    n_features_fallback_to_original = length(unique(fallback_features)),
    n_features_unscored = length(unique(unscored_features)),
    low_support_features = collapse_or_none(low_support_features),
    unavailable_features = collapse_or_none(unavailable_features),
    fallback_features = collapse_or_none(fallback_features),
    unscored_features = collapse_or_none(unscored_features),
    refine_status = "processed",
    refine_note = paste0(
      "low_support=", collapse_or_none(low_support_features),
      "; unavailable=", collapse_or_none(unavailable_features),
      "; fallback=", collapse_or_none(fallback_features),
      "; unscored=", collapse_or_none(unscored_features)
    ),
    stringsAsFactors = FALSE
  )

  rm(rankings, loo_auc_mat, refined_auc_mat, unit_expr, unit_meta, unit_cell_table)
  gc()
}

study_compartment_refine_manifest <- do.call(rbind, manifest_rows)
study_compartment_refine_manifest <- study_compartment_refine_manifest[order(
  study_compartment_refine_manifest$study,
  match(study_compartment_refine_manifest$project_broad_compartment, compartment_levels)
), , drop = FALSE]
write_tsv(study_compartment_refine_manifest, file.path(outdir, "study_compartment_refine_manifest.tsv"))

correlation_table <- do.call(rbind, correlation_rows)
correlation_table <- correlation_table[order(
  correlation_table$study,
  match(correlation_table$project_broad_compartment, compartment_levels),
  correlation_table$study_compartment_key,
  match(correlation_table$feature_id, feature_ids),
  correlation_table$gene
), , drop = FALSE]
write_tsv(correlation_table, file.path(outdir, "study_compartment_feature_gene_correlations.tsv.gz"), compress = TRUE)

feature_refine_summary <- do.call(rbind, feature_summary_rows)
feature_refine_summary <- feature_refine_summary[order(
  feature_refine_summary$study,
  match(feature_refine_summary$project_broad_compartment, compartment_levels),
  feature_refine_summary$study_compartment_key,
  match(feature_refine_summary$feature_id, feature_ids)
), , drop = FALSE]
write_tsv(feature_refine_summary, file.path(outdir, "study_compartment_feature_gene_set_refine_summary.tsv"))

if (length(refined_gene_set_rows) > 0L) {
  refined_gene_set_table <- do.call(rbind, refined_gene_set_rows)
  refined_gene_set_table <- refined_gene_set_table[order(
    refined_gene_set_table$study,
    match(refined_gene_set_table$project_broad_compartment, compartment_levels),
    refined_gene_set_table$study_compartment_key,
    match(refined_gene_set_table$feature_id, feature_ids),
    refined_gene_set_table$gene
  ), , drop = FALSE]
} else {
  refined_gene_set_table <- data.frame(
    study = character(0),
    cancer_type_summary = character(0),
    project_broad_compartment = character(0),
    compartment_role = character(0),
    study_compartment_key = character(0),
    feature_id = character(0),
    short_label = character(0),
    pathway = character(0),
    feature_refine_status = character(0),
    scoring_gene_set_source = character(0),
    gene = character(0),
    gene_membership_type = character(0),
    stringsAsFactors = FALSE
  )
}
write_tsv(refined_gene_set_table, file.path(outdir, "study_compartment_refined_gene_sets_long.tsv"))

refined_score_summary <- do.call(rbind, refined_score_summary_rows)
refined_score_summary <- refined_score_summary[order(
  refined_score_summary$study,
  match(refined_score_summary$project_broad_compartment, compartment_levels),
  refined_score_summary$study_compartment_key,
  match(refined_score_summary$feature_id, feature_ids)
), , drop = FALSE]
write_tsv(refined_score_summary, file.path(outdir, "study_compartment_feature_refined_score_summary.tsv"))

cell_level_rescored <- do.call(rbind, cell_rows)
cell_level_rescored <- cell_level_rescored[order(
  cell_level_rescored$study,
  match(cell_level_rescored$project_broad_compartment, compartment_levels),
  cell_level_rescored$sample_unit_key,
  cell_level_rescored$cell_id
), , drop = FALSE]
write_tsv(cell_level_rescored, file.path(outdir, "cell_level_lm12_refined_scores.tsv.gz"), compress = TRUE)

refine_bundle <- list(
  parameters = list(
    score_method = "AUCell",
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    auc_rank_fraction = auc_rank_fraction,
    negative_correlation_threshold = negative_correlation_threshold,
    min_cells_for_correlation = min_cells_for_correlation,
    low_support_gene_threshold = low_support_gene_threshold,
    min_retained_genes_for_direct_refined_scoring = min_retained_genes_for_direct_refined_scoring,
    allow_original_fallback = allow_original_fallback,
    study_filter = study_filter,
    study_compartment_filter = study_compartment_filter
  ),
  feature_metadata = feature_metadata,
  gene_sets_present = gene_sets_present,
  study_compartment_refine_manifest = study_compartment_refine_manifest,
  study_compartment_feature_gene_correlations = correlation_table,
  study_compartment_feature_gene_set_refine_summary = feature_refine_summary,
  study_compartment_refined_gene_sets_long = refined_gene_set_table,
  study_compartment_feature_refined_score_summary = refined_score_summary,
  cell_level_refined_score_table = cell_level_rescored,
  retained_gene_sets_by_unit = retained_gene_sets_by_unit,
  scoring_gene_sets_by_unit = scoring_gene_sets_by_unit,
  sample_metadata = sample_metadata,
  celltype_to_broad = celltype_to_broad,
  compartment_manifest = compartment_manifest,
  sample_inclusion_rule = sample_inclusion_rule,
  treatment_map = treatment_map
)
saveRDS(refine_bundle, file.path(outdir, "studywise_refine_bundle.rds"))

summary_lines <- c(
  "03_studywise_within_compartment_lm12_refine_and_rescore completed.",
  paste0("Input Seurat object: ", input_rds),
  paste0("Input bundle: ", input_bundle),
  paste0("Score assay: ", score_assay),
  paste0("Score slot used: ", score_slot_used),
  paste0("Score method: AUCell"),
  paste0("AUC rank fraction: ", auc_rank_fraction),
  paste0("Negative correlation threshold: rho < ", negative_correlation_threshold),
  paste0("Min cells for gene-score correlation: ", min_cells_for_correlation),
  paste0("Low-support gene threshold: retained genes < ", low_support_gene_threshold),
  paste0("Min retained genes for direct refined scoring: ", min_retained_genes_for_direct_refined_scoring),
  paste0("Allow original fallback: ", allow_original_fallback),
  paste0("Studies processed: ", length(unique(study_compartment_refine_manifest$study))),
  paste0("Study x compartment units processed: ", nrow(study_compartment_refine_manifest)),
  paste0("Cells rescored: ", nrow(cell_level_rescored)),
  paste0("Output cell-level table: ", file.path(outdir, "cell_level_lm12_refined_scores.tsv.gz")),
  paste0("Output bundle: ", file.path(outdir, "studywise_refine_bundle.rds"))
)
writeLines(summary_lines, con = file.path(outdir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
}

# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/01_assign_study_compartment_ig_states.R
stage_single_cell_01_assign_study_compartment_ig_states <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("refine_bundle")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# Source-preserved defaults; no study filter is applied.
options(stringsAsFactors = FALSE)
outdir <- out_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
input_bundle <- unname(paths["refine_bundle"])
k_clusters <- 2L
min_usable_features_for_clustering <- 6L
feature_status_for_clustering <- "usable"
kmeans_nstart <- 100L
kmeans_iter_max <- 200L
random_seed <- 20260615L
min_cells_for_clustering <- 50L
study_filter_arg <- ""
study_compartment_filter_arg <- ""
`write_tsv` <- bind_stage_helper(shared_write_tsv_006, environment())
`safe_char` <- bind_stage_helper(shared_safe_char_042, environment())
`round_or_na` <- bind_stage_helper(shared_round_or_na_045, environment())
`collapse_or_none` <- bind_stage_helper(shared_collapse_or_none_043, environment())
`parse_filter_values` <- bind_stage_helper(shared_parse_filter_values_040, environment())
`%notin%` <- bind_stage_helper(shared__notin__055, environment())

required_files <- c(input_bundle)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

message("Loading 03 studywise refine bundle: ", input_bundle)
bundle <- readRDS(input_bundle)

required_components <- c(
  "feature_metadata",
  "study_compartment_refine_manifest",
  "study_compartment_feature_gene_set_refine_summary",
  "cell_level_refined_score_table"
)
missing_components <- required_components[required_components %notin% names(bundle)]
if (length(missing_components) > 0L) {
  stop("03 refine bundle is missing required components: ", paste(missing_components, collapse = ", "))
}

feature_metadata <- bundle$feature_metadata
feature_ids <- safe_char(feature_metadata$feature_id)
study_compartment_refine_manifest <- bundle$study_compartment_refine_manifest
feature_refine_summary <- bundle$study_compartment_feature_gene_set_refine_summary
cell_level_scores <- bundle$cell_level_refined_score_table

immune_axis_features <- setdiff(feature_ids, "LM12_Antimicrobial_Peptides")

required_cell_cols <- c(
  "study_compartment_key", "cell_id", "study", "cancer_type", "sample_unit_key",
  "project_broad_compartment", "compartment_role"
)
missing_cell_cols <- required_cell_cols[required_cell_cols %notin% colnames(cell_level_scores)]
if (length(missing_cell_cols) > 0L) {
  stop("cell_level_refined_score_table is missing required columns: ", paste(missing_cell_cols, collapse = ", "))
}

refined_score_cols <- paste0("refined__", feature_ids)
missing_score_cols <- refined_score_cols[refined_score_cols %notin% colnames(cell_level_scores)]
if (length(missing_score_cols) > 0L) {
  stop("cell_level_refined_score_table is missing refined score columns: ", paste(missing_score_cols, collapse = ", "))
}

study_filter <- parse_filter_values(study_filter_arg)
study_compartment_filter <- parse_filter_values(study_compartment_filter_arg)

if (length(study_filter) > 0L) {
  keep_units <- safe_char(study_compartment_refine_manifest$study) %in% study_filter
  study_compartment_refine_manifest <- study_compartment_refine_manifest[keep_units, , drop = FALSE]
}
if (length(study_compartment_filter) > 0L) {
  keep_units <- safe_char(study_compartment_refine_manifest$study_compartment_key) %in% study_compartment_filter
  study_compartment_refine_manifest <- study_compartment_refine_manifest[keep_units, , drop = FALSE]
}
if (nrow(study_compartment_refine_manifest) == 0L) {
  stop("No study-compartment units remain after filtering.")
}

unit_keys <- safe_char(study_compartment_refine_manifest$study_compartment_key)
feature_refine_summary <- feature_refine_summary[
  safe_char(feature_refine_summary$study_compartment_key) %in% unit_keys,
  ,
  drop = FALSE
]
cell_level_scores <- cell_level_scores[
  safe_char(cell_level_scores$study_compartment_key) %in% unit_keys,
  ,
  drop = FALSE
]

metadata_cols <- intersect(c(
  "study_compartment_key", "cell_id", "study", "cancer_type", "donor_id", "Sample",
  "donor_id_pre_post", "sample_unit_key", "pre_post", "outcome", "Combined_outcome",
  "treatment_raw", "treatment_class_broad", "cell_types", "cell_types_update",
  "Cell_type_broad_original", "project_broad_compartment", "compartment_role",
  "nCount_RNA", "nFeature_RNA"
), colnames(cell_level_scores))

cell_state_rows <- list()
cluster_summary_rows <- list()
cluster_center_rows <- list()
manifest_rows <- list()
sample_fraction_rows <- list()

split_keys <- split(seq_len(nrow(cell_level_scores)), cell_level_scores$study_compartment_key)

for (idx in seq_len(nrow(study_compartment_refine_manifest))) {
  manifest_row <- study_compartment_refine_manifest[idx, , drop = FALSE]
  key <- safe_char(manifest_row$study_compartment_key[1])
  message("Assigning IG2 states for study-compartment: ", key)

  cell_idx <- split_keys[[key]]
  df <- cell_level_scores[cell_idx, , drop = FALSE]
  df <- df[order(df$sample_unit_key, df$cell_id), , drop = FALSE]

  unit_feature_status <- feature_refine_summary[
    safe_char(feature_refine_summary$study_compartment_key) == key,
    ,
    drop = FALSE
  ]
  unit_feature_status <- unit_feature_status[order(match(unit_feature_status$feature_id, feature_ids)), , drop = FALSE]

  usable_features <- safe_char(unit_feature_status$feature_id[
    unit_feature_status$feature_refine_status == feature_status_for_clustering
  ])
  low_support_features <- safe_char(unit_feature_status$feature_id[
    unit_feature_status$feature_refine_status == "low_support"
  ])
  unavailable_features <- safe_char(unit_feature_status$feature_id[
    unit_feature_status$feature_refine_status == "unavailable"
  ])
  fallback_features <- safe_char(unit_feature_status$feature_id[
    grepl("^original_fallback", unit_feature_status$scoring_gene_set_source)
  ])
  unscored_features <- safe_char(unit_feature_status$feature_id[
    unit_feature_status$scoring_gene_set_source == "unscored"
  ])

  cell_out <- df[, metadata_cols, drop = FALSE]
  axis_cols <- paste0("refined__", immune_axis_features)
  axis_mat <- as.matrix(df[, axis_cols, drop = FALSE])
  storage.mode(axis_mat) <- "numeric"
  axis_counts <- rowSums(is.finite(axis_mat))
  immune_axis <- rowMeans(axis_mat, na.rm = TRUE)
  immune_axis[axis_counts == 0L] <- NA_real_

  cell_out$lm12_inflammatory_axis_04 <- round_or_na(immune_axis)
  cell_out$cluster_status_04 <- "not_run"
  cell_out$clusterable_for_04 <- FALSE
  cell_out$n_features_usable_04 <- length(usable_features)
  cell_out$n_features_used_04 <- 0L
  cell_out$ig2_raw_cluster <- NA_character_
  cell_out$ig2_raw_cluster_id <- NA_character_
  cell_out$IG <- NA_character_
  cell_out$IG_Label <- NA_character_
  cell_out$IG_axis_rank <- NA_integer_

  zero_variance_features <- character(0)
  features_used <- character(0)
  cluster_status <- "processed"
  cluster_note <- "Study-wise within-compartment IG2 clustering completed."
  clusterable_cells <- rep(FALSE, nrow(df))
  n_complete_cells <- 0L
  kmeans_total_withinss <- NA_real_
  kmeans_betweenss_fraction <- NA_real_
  cluster_size_min <- NA_integer_
  cluster_size_max <- NA_integer_

  if (nrow(df) < min_cells_for_clustering) {
    cluster_status <- "skipped_few_cells"
    cluster_note <- paste0("Fewer than ", min_cells_for_clustering, " cells.")
  } else if (length(usable_features) < min_usable_features_for_clustering) {
    cluster_status <- "skipped_few_usable_features"
    cluster_note <- paste0("Fewer than ", min_usable_features_for_clustering, " usable refined features.")
  } else {
    score_cols <- paste0("refined__", usable_features)
    score_mat <- as.matrix(df[, score_cols, drop = FALSE])
    storage.mode(score_mat) <- "numeric"
    colnames(score_mat) <- usable_features

    clusterable_cells <- stats::complete.cases(score_mat)
    n_complete_cells <- sum(clusterable_cells)
    if (n_complete_cells < min_cells_for_clustering) {
      cluster_status <- "skipped_few_complete_cells"
      cluster_note <- paste0("Fewer than ", min_cells_for_clustering, " complete-score cells.")
    } else if (n_complete_cells < k_clusters) {
      cluster_status <- "skipped_lt_k_complete_cells"
      cluster_note <- paste0("Fewer than ", k_clusters, " complete-score cells.")
    } else {
      score_mat_complete <- score_mat[clusterable_cells, , drop = FALSE]
      feature_sds <- apply(score_mat_complete, 2, stats::sd)
      features_used <- names(feature_sds)[is.finite(feature_sds) & feature_sds > 0]
      zero_variance_features <- setdiff(usable_features, features_used)

      if (length(features_used) < min_usable_features_for_clustering) {
        cluster_status <- "skipped_few_variable_features"
        cluster_note <- paste0("Fewer than ", min_usable_features_for_clustering, " non-zero-variance usable features.")
      } else {
        score_mat_used <- score_mat_complete[, features_used, drop = FALSE]
        scaled_mat <- scale(score_mat_used)

        set.seed(random_seed + idx)
        kmeans_res <- stats::kmeans(
          x = scaled_mat,
          centers = k_clusters,
          nstart = kmeans_nstart,
          iter.max = kmeans_iter_max
        )

        raw_clusters <- paste0("K", kmeans_res$cluster)
        raw_cluster_levels <- paste0("K", seq_len(k_clusters))
        cluster_axis_means <- vapply(raw_cluster_levels, function(cl) {
          vals <- immune_axis[clusterable_cells][raw_clusters == cl]
          if (all(!is.finite(vals))) {
            return(NA_real_)
          }
          mean(vals, na.rm = TRUE)
        }, numeric(1))
        mean_z_by_cluster <- rowMeans(kmeans_res$centers)
        if (any(!is.finite(cluster_axis_means)) || length(unique(round(cluster_axis_means, 12))) < k_clusters) {
          cluster_axis_means <- mean_z_by_cluster
        }
        rank_by_axis <- rank(-cluster_axis_means, ties.method = "first")
        names(rank_by_axis) <- raw_cluster_levels
        ig_lookup <- ifelse(rank_by_axis == 1L, "IG1", "IG2")
        names(ig_lookup) <- raw_cluster_levels
        ig_label_lookup <- c(
          IG1 = "lm12_inflammatory_remodeling_high",
          IG2 = "lm12_inflammatory_remodeling_low"
        )

        cell_out$cluster_status_04 <- "clustered"
        cell_out$clusterable_for_04 <- clusterable_cells
        cell_out$n_features_used_04 <- length(features_used)
        cell_out$ig2_raw_cluster[clusterable_cells] <- raw_clusters
        cell_out$ig2_raw_cluster_id[clusterable_cells] <- paste(key, raw_clusters, sep = "__")
        cell_out$IG[clusterable_cells] <- ig_lookup[raw_clusters]
        cell_out$IG_Label[clusterable_cells] <- ig_label_lookup[cell_out$IG[clusterable_cells]]
        cell_out$IG_axis_rank[clusterable_cells] <- as.integer(rank_by_axis[raw_clusters])

        kmeans_total_withinss <- kmeans_res$tot.withinss
        kmeans_betweenss_fraction <- if (is.finite(kmeans_res$totss) && kmeans_res$totss > 0) {
          kmeans_res$betweenss / kmeans_res$totss
        } else {
          NA_real_
        }

        size_tab <- table(factor(raw_clusters, levels = raw_cluster_levels))
        cluster_size_min <- min(as.integer(size_tab))
        cluster_size_max <- max(as.integer(size_tab))

        for (cluster_name in raw_cluster_levels) {
          cells_in_cluster <- which(clusterable_cells)[raw_clusters == cluster_name]
          ig_state <- ig_lookup[[cluster_name]]
          ig_label <- ig_label_lookup[[ig_state]]

          cluster_summary_rows[[paste(key, cluster_name, sep = "__")]] <- data.frame(
            study = manifest_row$study,
            cancer_type_summary = manifest_row$cancer_type_summary,
            project_broad_compartment = manifest_row$project_broad_compartment,
            compartment_role = manifest_row$compartment_role,
            study_compartment_key = key,
            ig2_raw_cluster = cluster_name,
            ig2_raw_cluster_id = paste(key, cluster_name, sep = "__"),
            IG = ig_state,
            IG_Label = ig_label,
            IG_axis_rank = as.integer(rank_by_axis[[cluster_name]]),
            n_cells = length(cells_in_cluster),
            fraction_of_clusterable_cells = round_or_na(length(cells_in_cluster) / n_complete_cells),
            n_sample_units = length(unique(df$sample_unit_key[cells_in_cluster])),
            mean_lm12_inflammatory_axis = round_or_na(mean(immune_axis[cells_in_cluster], na.rm = TRUE)),
            median_lm12_inflammatory_axis = round_or_na(stats::median(immune_axis[cells_in_cluster], na.rm = TRUE)),
            mean_used_z_center = round_or_na(mean_z_by_cluster[as.integer(sub("^K", "", cluster_name))]),
            stringsAsFactors = FALSE
          )

          for (fid in feature_ids) {
            raw_col <- paste0("refined__", fid)
            raw_values <- df[[raw_col]][cells_in_cluster]
            raw_values <- raw_values[is.finite(raw_values)]
            feature_status_now <- unit_feature_status$feature_refine_status[match(fid, unit_feature_status$feature_id)]
            scoring_source_now <- unit_feature_status$scoring_gene_set_source[match(fid, unit_feature_status$feature_id)]
            cluster_center_rows[[paste(key, cluster_name, fid, sep = "__")]] <- data.frame(
              study = manifest_row$study,
              cancer_type_summary = manifest_row$cancer_type_summary,
              project_broad_compartment = manifest_row$project_broad_compartment,
              compartment_role = manifest_row$compartment_role,
              study_compartment_key = key,
              ig2_raw_cluster = cluster_name,
              ig2_raw_cluster_id = paste(key, cluster_name, sep = "__"),
              IG = ig_state,
              IG_Label = ig_label,
              IG_axis_rank = as.integer(rank_by_axis[[cluster_name]]),
              feature_id = fid,
              short_label = feature_metadata$short_label[match(fid, feature_metadata$feature_id)],
              pathway = feature_metadata$pathway[match(fid, feature_metadata$feature_id)],
              feature_refine_status = feature_status_now,
              scoring_gene_set_source = scoring_source_now,
              used_for_clustering_04 = fid %in% features_used,
              mean_refined_score = if (length(raw_values) == 0L) NA_real_ else round_or_na(mean(raw_values)),
              mean_scaled_score = if (fid %in% features_used) {
                round_or_na(kmeans_res$centers[as.integer(sub("^K", "", cluster_name)), fid])
              } else {
                NA_real_
              },
              stringsAsFactors = FALSE
            )
          }
        }

        sample_tmp <- cell_out[cell_out$clusterable_for_04, c(
          "study", "cancer_type", "sample_unit_key", "project_broad_compartment",
          "study_compartment_key", "IG", "IG_Label"
        ), drop = FALSE]
        sample_tmp$n_cells <- 1L
        count_df <- stats::aggregate(
          n_cells ~ study + cancer_type + sample_unit_key + project_broad_compartment +
            study_compartment_key + IG + IG_Label,
          data = sample_tmp,
          FUN = sum
        )
        total_df <- stats::aggregate(
          n_cells ~ study + cancer_type + sample_unit_key + project_broad_compartment +
            study_compartment_key,
          data = sample_tmp,
          FUN = sum
        )
        colnames(total_df)[colnames(total_df) == "n_cells"] <- "n_clusterable_cells_in_sample_compartment"
        count_df <- merge(
          count_df,
          total_df,
          by = c("study", "cancer_type", "sample_unit_key", "project_broad_compartment", "study_compartment_key"),
          all.x = TRUE,
          sort = FALSE
        )
        count_df$fraction_of_sample_compartment <- round_or_na(
          count_df$n_cells / count_df$n_clusterable_cells_in_sample_compartment
        )
        sample_fraction_rows[[key]] <- count_df[order(
          count_df$sample_unit_key,
          count_df$IG
        ), , drop = FALSE]
      }
    }
  }

  if (cluster_status != "processed") {
    cell_out$cluster_status_04 <- cluster_status
    cell_out$clusterable_for_04 <- clusterable_cells
    cell_out$n_features_used_04 <- length(features_used)
  }

  manifest_rows[[key]] <- data.frame(
    study = manifest_row$study,
    cancer_type_summary = manifest_row$cancer_type_summary,
    project_broad_compartment = manifest_row$project_broad_compartment,
    compartment_role = manifest_row$compartment_role,
    study_compartment_key = key,
    n_cells = nrow(df),
    n_sample_units = manifest_row$n_sample_units,
    n_donors = manifest_row$n_donors,
    n_unique_cell_types_update = manifest_row$n_unique_cell_types_update,
    k_clusters = k_clusters,
    cluster_status = cluster_status,
    cluster_note = cluster_note,
    n_clusterable_cells = n_complete_cells,
    n_features_usable = length(usable_features),
    n_features_used = length(features_used),
    min_cells_for_clustering = min_cells_for_clustering,
    min_usable_features_for_clustering = min_usable_features_for_clustering,
    feature_status_for_clustering = feature_status_for_clustering,
    kmeans_nstart = kmeans_nstart,
    kmeans_iter_max = kmeans_iter_max,
    random_seed_unit = random_seed + idx,
    kmeans_total_withinss = round_or_na(kmeans_total_withinss),
    kmeans_betweenss_fraction = round_or_na(kmeans_betweenss_fraction),
    cluster_size_min = cluster_size_min,
    cluster_size_max = cluster_size_max,
    low_support_features = collapse_or_none(low_support_features),
    unavailable_features = collapse_or_none(unavailable_features),
    fallback_features = collapse_or_none(fallback_features),
    unscored_features = collapse_or_none(unscored_features),
    zero_variance_features = collapse_or_none(zero_variance_features),
    features_used = collapse_or_none(features_used),
    stringsAsFactors = FALSE
  )

  cell_state_rows[[key]] <- cell_out
}

cell_level_ig_states <- do.call(rbind, cell_state_rows)
ig_state_manifest <- do.call(rbind, manifest_rows)
ig_cluster_summary <- if (length(cluster_summary_rows) > 0L) {
  do.call(rbind, cluster_summary_rows)
} else {
  data.frame()
}
ig_cluster_centers <- if (length(cluster_center_rows) > 0L) {
  do.call(rbind, cluster_center_rows)
} else {
  data.frame()
}
sample_compartment_ig_fractions <- if (length(sample_fraction_rows) > 0L) {
  do.call(rbind, sample_fraction_rows)
} else {
  data.frame()
}

write_tsv(
  cell_level_ig_states,
  file.path(outdir, "cell_level_local_ig_states.tsv.gz"),
  compress = TRUE
)
write_tsv(
  ig_state_manifest,
  file.path(outdir, "study_compartment_ig_state_manifest.tsv")
)
write_tsv(
  ig_cluster_summary,
  file.path(outdir, "study_compartment_ig_cluster_summary.tsv")
)
write_tsv(
  ig_cluster_centers,
  file.path(outdir, "study_compartment_ig_cluster_centers.tsv")
)
write_tsv(
  sample_compartment_ig_fractions,
  file.path(outdir, "sample_compartment_ig_cell_fractions.tsv")
)

bundle_out <- list(
  parameters = list(
    input_bundle = input_bundle,
    outdir = outdir,
    k_clusters = k_clusters,
    kmeans_nstart = kmeans_nstart,
    kmeans_iter_max = kmeans_iter_max,
    random_seed = random_seed,
    min_cells_for_clustering = min_cells_for_clustering,
    min_usable_features_for_clustering = min_usable_features_for_clustering,
    feature_status_for_clustering = feature_status_for_clustering,
    immune_axis_features = immune_axis_features
  ),
  feature_metadata = feature_metadata,
  study_compartment_ig_state_manifest = ig_state_manifest,
  study_compartment_ig_cluster_summary = ig_cluster_summary,
  study_compartment_ig_cluster_centers = ig_cluster_centers,
  sample_compartment_ig_cell_fractions = sample_compartment_ig_fractions,
  cell_level_ig_states = cell_level_ig_states
)
saveRDS(bundle_out, file.path(outdir, "studywise_within_compartment_ig2_bundle.rds"))

run_summary_path <- file.path(outdir, "run_summary.txt")
summary_lines <- c(
  "04_studywise_within_compartment_ig2_state_assignment completed.",
  paste0("Input bundle: ", input_bundle),
  paste0("K clusters: ", k_clusters),
  paste0("Feature status for clustering: ", feature_status_for_clustering),
  paste0("Min cells for clustering: ", min_cells_for_clustering),
  paste0("Min usable features for clustering: ", min_usable_features_for_clustering),
  paste0("Studies processed: ", length(unique(ig_state_manifest$study))),
  paste0("Study x compartment units processed: ", sum(ig_state_manifest$cluster_status == "processed")),
  paste0("Study x compartment units skipped: ", sum(ig_state_manifest$cluster_status != "processed")),
  paste0("Cells assigned IG: ", sum(!is.na(cell_level_ig_states$IG))),
  paste0("Output cell-level table: ", file.path(outdir, "cell_level_local_ig_states.tsv.gz")),
  paste0("Output bundle: ", file.path(outdir, "studywise_within_compartment_ig2_bundle.rds"))
)
writeLines(summary_lines, con = run_summary_path)

message(paste(summary_lines, collapse = "\n"))
}

# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/02_score_acid_assign_ag_states.R
stage_single_cell_02_score_acid_assign_ag_states <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "refine_bundle", "selected_genes")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
input_refine_bundle <- unname(paths["refine_bundle"])
acid_gene_file <- unname(paths["selected_genes"])
score_assay <- "RNA"
rank_plot_stats <- FALSE
rank_split_by_blocks <- TRUE
rank_keep_zeroes_as_na <- FALSE
auc_norm <- TRUE
auc_rank_fraction <- 0.05
k_clusters <- 3L
min_present_acid_genes <- 1L
kmeans_nstart <- 100L
kmeans_iter_max <- 200L
random_seed <- 20260615L
min_cells_for_clustering <- 50L
study_filter_arg <- ""
study_compartment_filter_arg <- ""
`write_tsv` <- bind_stage_helper(shared_write_tsv_006, environment())
`safe_char` <- bind_stage_helper(shared_safe_char_042, environment())
`round_or_na` <- bind_stage_helper(shared_round_or_na_045, environment())
`collapse_or_none` <- bind_stage_helper(shared_collapse_or_none_043, environment())
`parse_filter_values` <- bind_stage_helper(shared_parse_filter_values_040, environment())
`%notin%` <- bind_stage_helper(shared__notin__055, environment())
`extract_scoring_matrix` <- bind_stage_helper(shared_extract_scoring_matrix_035, environment())

score_gene_sets_auc <- function(rankings, gene_sets, auc_max_rank, auc_norm = TRUE) {
  auc_obj <- AUCell::AUCell_calcAUC(
    geneSets = gene_sets,
    rankings = rankings,
    nCores = 1,
    normAUC = auc_norm,
    aucMaxRank = auc_max_rank,
    verbose = FALSE
  )
  auc_mat <- AUCell::getAUC(auc_obj)
  auc_mat <- t(as.matrix(auc_mat))
  auc_mat
}

required_files <- c(input_rds, input_refine_bundle, acid_gene_file)
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

message("Loading 03 studywise refine bundle: ", input_refine_bundle)
refine_bundle <- readRDS(input_refine_bundle)
required_components <- c(
  "study_compartment_refine_manifest",
  "cell_level_refined_score_table",
  "sample_metadata"
)
missing_components <- required_components[required_components %notin% names(refine_bundle)]
if (length(missing_components) > 0L) {
  stop("03 refine bundle is missing required components: ", paste(missing_components, collapse = ", "))
}

study_compartment_manifest <- refine_bundle$study_compartment_refine_manifest
cell_metadata <- refine_bundle$cell_level_refined_score_table
sample_metadata <- refine_bundle$sample_metadata

required_cell_cols <- c(
  "study_compartment_key", "cell_id", "study", "cancer_type", "sample_unit_key",
  "project_broad_compartment", "compartment_role"
)
missing_cell_cols <- required_cell_cols[required_cell_cols %notin% colnames(cell_metadata)]
if (length(missing_cell_cols) > 0L) {
  stop("cell_level_refined_score_table is missing required columns: ", paste(missing_cell_cols, collapse = ", "))
}
if (!all(cell_metadata$cell_id %in% obj_cell_ids)) {
  stop("Some cell IDs in the 03 refine bundle are not present in the Seurat object.")
}

# 隔离加载参考对象，避免覆盖命令行配置和路径保护变量。
selected_gene_env <- new.env(parent = emptyenv())
loaded_selected_names <- load(acid_gene_file, envir = selected_gene_env)
if (!"final_selected_genes" %in% loaded_selected_names) {
  stop("Object `final_selected_genes` not found in ", acid_gene_file)
}
final_selected_genes <- get("final_selected_genes", envir = selected_gene_env, inherits = FALSE)
if (!"Gene" %in% colnames(final_selected_genes)) {
  stop("`final_selected_genes` must contain a Gene column.")
}
acid_genes <- sort(unique(safe_char(final_selected_genes$Gene)))
acid_genes <- acid_genes[acid_genes != "" & acid_genes != "NA"]
if (length(acid_genes) == 0L) {
  stop("No usable acid genes found in final_selected_genes.")
}

study_filter <- parse_filter_values(study_filter_arg)
study_compartment_filter <- parse_filter_values(study_compartment_filter_arg)

if (length(study_filter) > 0L) {
  keep_units <- safe_char(study_compartment_manifest$study) %in% study_filter
  study_compartment_manifest <- study_compartment_manifest[keep_units, , drop = FALSE]
}
if (length(study_compartment_filter) > 0L) {
  keep_units <- safe_char(study_compartment_manifest$study_compartment_key) %in% study_compartment_filter
  study_compartment_manifest <- study_compartment_manifest[keep_units, , drop = FALSE]
}
if (nrow(study_compartment_manifest) == 0L) {
  stop("No study-compartment units remain after filtering.")
}

unit_keys <- safe_char(study_compartment_manifest$study_compartment_key)
cell_metadata <- cell_metadata[
  safe_char(cell_metadata$study_compartment_key) %in% unit_keys,
  ,
  drop = FALSE
]
if (nrow(cell_metadata) == 0L) {
  stop("No cell-level rows remain after filtering the 03 refine bundle.")
}

metadata_cols <- intersect(c(
  "study_compartment_key", "cell_id", "study", "cancer_type", "donor_id", "Sample",
  "donor_id_pre_post", "sample_unit_key", "pre_post", "outcome", "Combined_outcome",
  "treatment_raw", "treatment_class_broad", "cell_types", "cell_types_update",
  "Cell_type_broad_original", "project_broad_compartment", "compartment_role",
  "nCount_RNA", "nFeature_RNA"
), colnames(cell_metadata))

expr_genes <- rownames(expr_mat)
acid_genes_present_overall <- intersect(acid_genes, expr_genes)
if (length(acid_genes_present_overall) < min_present_acid_genes) {
  stop(
    "Fewer than ", min_present_acid_genes,
    " acid genes are present in the selected assay."
  )
}

cell_rows <- list()
coverage_rows <- list()
gene_coverage_rows <- list()
cluster_summary_rows <- list()
manifest_rows <- list()
sample_fraction_rows <- list()

split_keys <- split(seq_len(nrow(cell_metadata)), cell_metadata$study_compartment_key)

for (idx in seq_len(nrow(study_compartment_manifest))) {
  manifest_row <- study_compartment_manifest[idx, , drop = FALSE]
  key <- safe_char(manifest_row$study_compartment_key[1])
  message("Scoring acid and assigning AG3 for study-compartment: ", key)

  cell_idx <- split_keys[[key]]
  df <- cell_metadata[cell_idx, , drop = FALSE]
  df <- df[order(df$sample_unit_key, df$cell_id), , drop = FALSE]
  unit_cells <- df$cell_id
  unit_expr <- expr_mat[, unit_cells, drop = FALSE]

  acid_genes_present <- intersect(acid_genes, rownames(unit_expr))
  acid_genes_missing <- setdiff(acid_genes, acid_genes_present)
  auc_max_rank <- ceiling(auc_rank_fraction * nrow(unit_expr))
  auc_max_rank <- min(auc_max_rank, nrow(unit_expr))

  coverage_rows[[key]] <- data.frame(
    study = manifest_row$study,
    cancer_type_summary = manifest_row$cancer_type_summary,
    project_broad_compartment = manifest_row$project_broad_compartment,
    compartment_role = manifest_row$compartment_role,
    study_compartment_key = key,
    n_cells = nrow(df),
    n_sample_units = manifest_row$n_sample_units,
    n_donors = manifest_row$n_donors,
    n_unique_cell_types_update = manifest_row$n_unique_cell_types_update,
    acid_genes_total = length(acid_genes),
    acid_genes_present = length(acid_genes_present),
    acid_gene_coverage = round_or_na(length(acid_genes_present) / length(acid_genes)),
    acid_genes_present_list = collapse_or_none(acid_genes_present),
    acid_genes_missing_list = collapse_or_none(acid_genes_missing),
    score_method = "AUCell",
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    auc_max_rank = auc_max_rank,
    stringsAsFactors = FALSE
  )

  for (gene in acid_genes) {
    if (gene %in% rownames(unit_expr)) {
      expr_vec <- unit_expr[gene, , drop = TRUE]
      detected_cells <- as.integer(sum(expr_vec > 0))
    } else {
      detected_cells <- NA_integer_
    }
    gene_coverage_rows[[paste(key, gene, sep = "__")]] <- data.frame(
      study = manifest_row$study,
      cancer_type_summary = manifest_row$cancer_type_summary,
      project_broad_compartment = manifest_row$project_broad_compartment,
      compartment_role = manifest_row$compartment_role,
      study_compartment_key = key,
      gene = gene,
      gene_present_in_assay = gene %in% rownames(unit_expr),
      n_cells = nrow(df),
      detected_cells = detected_cells,
      detected_fraction = if (is.na(detected_cells)) NA_real_ else round_or_na(detected_cells / nrow(df)),
      stringsAsFactors = FALSE
    )
  }

  cell_out <- df[, metadata_cols, drop = FALSE]
  cell_out$Acid_Score_05 <- NA_real_
  cell_out$acid_score_status_05 <- "not_run"
  cell_out$cluster_status_05 <- "not_run"
  cell_out$clusterable_for_05 <- FALSE
  cell_out$acid_raw_cluster <- NA_character_
  cell_out$acid_raw_cluster_id <- NA_character_
  cell_out$AG <- NA_character_
  cell_out$AG_Label <- NA_character_
  cell_out$AG_acid_rank <- NA_integer_

  cluster_status <- "processed"
  cluster_note <- "Study-wise within-compartment AG3 clustering completed."
  score_status <- "scored"
  n_scored_cells <- 0L
  n_clusterable_cells <- 0L
  kmeans_total_withinss <- NA_real_
  kmeans_betweenss_fraction <- NA_real_
  cluster_size_min <- NA_integer_
  cluster_size_max <- NA_integer_

  if (length(acid_genes_present) < min_present_acid_genes) {
    score_status <- "skipped_low_gene_coverage"
    cluster_status <- "skipped_low_gene_coverage"
    cluster_note <- paste0("Fewer than ", min_present_acid_genes, " acid genes present.")
  } else {
    rankings <- AUCell::AUCell_buildRankings(
      exprMat = unit_expr,
      featureType = "genes",
      plotStats = rank_plot_stats,
      splitByBlocks = rank_split_by_blocks,
      keepZeroesAsNA = rank_keep_zeroes_as_na,
      verbose = FALSE
    )
    auc_mat <- score_gene_sets_auc(
      rankings = rankings,
      gene_sets = list(Acid_Score = acid_genes_present),
      auc_max_rank = auc_max_rank,
      auc_norm = auc_norm
    )
    acid_score <- as.numeric(auc_mat[unit_cells, "Acid_Score"])
    names(acid_score) <- unit_cells
    cell_out$Acid_Score_05 <- round_or_na(acid_score)
    cell_out$acid_score_status_05 <- score_status
    n_scored_cells <- sum(is.finite(acid_score))

    if (nrow(df) < min_cells_for_clustering) {
      cluster_status <- "skipped_few_cells"
      cluster_note <- paste0("Fewer than ", min_cells_for_clustering, " cells.")
    } else {
      clusterable_cells <- is.finite(acid_score)
      n_clusterable_cells <- sum(clusterable_cells)
      if (n_clusterable_cells < min_cells_for_clustering) {
        cluster_status <- "skipped_few_complete_cells"
        cluster_note <- paste0("Fewer than ", min_cells_for_clustering, " cells with finite acid score.")
      } else if (length(unique(round(acid_score[clusterable_cells], 12))) < k_clusters) {
        cluster_status <- "skipped_few_unique_scores"
        cluster_note <- paste0("Fewer than ", k_clusters, " unique acid scores.")
      } else {
        score_mat <- matrix(acid_score[clusterable_cells], ncol = 1L)
        colnames(score_mat) <- "Acid_Score"

        set.seed(random_seed + idx)
        kmeans_res <- stats::kmeans(
          x = score_mat,
          centers = k_clusters,
          nstart = kmeans_nstart,
          iter.max = kmeans_iter_max
        )

        raw_clusters <- paste0("K", kmeans_res$cluster)
        raw_cluster_levels <- paste0("K", seq_len(k_clusters))
        cluster_acid_means <- vapply(raw_cluster_levels, function(cl) {
          vals <- acid_score[clusterable_cells][raw_clusters == cl]
          if (all(!is.finite(vals))) {
            return(NA_real_)
          }
          mean(vals, na.rm = TRUE)
        }, numeric(1))
        rank_by_acid <- rank(-cluster_acid_means, ties.method = "first")
        names(rank_by_acid) <- raw_cluster_levels
        ag_lookup <- c("AG1", "AG2", "AG3")[rank_by_acid]
        names(ag_lookup) <- raw_cluster_levels
        ag_label_lookup <- c(
          AG1 = "acid_high",
          AG2 = "acid_intermediate",
          AG3 = "acid_low"
        )

        cell_out$cluster_status_05 <- "clustered"
        cell_out$clusterable_for_05 <- clusterable_cells
        cell_out$acid_raw_cluster[clusterable_cells] <- raw_clusters
        cell_out$acid_raw_cluster_id[clusterable_cells] <- paste(key, raw_clusters, sep = "__")
        cell_out$AG[clusterable_cells] <- ag_lookup[raw_clusters]
        cell_out$AG_Label[clusterable_cells] <- ag_label_lookup[cell_out$AG[clusterable_cells]]
        cell_out$AG_acid_rank[clusterable_cells] <- as.integer(rank_by_acid[raw_clusters])

        kmeans_total_withinss <- kmeans_res$tot.withinss
        kmeans_betweenss_fraction <- if (is.finite(kmeans_res$totss) && kmeans_res$totss > 0) {
          kmeans_res$betweenss / kmeans_res$totss
        } else {
          NA_real_
        }

        size_tab <- table(factor(raw_clusters, levels = raw_cluster_levels))
        cluster_size_min <- min(as.integer(size_tab))
        cluster_size_max <- max(as.integer(size_tab))

        for (cluster_name in raw_cluster_levels) {
          cells_in_cluster <- which(clusterable_cells)[raw_clusters == cluster_name]
          ag_state <- ag_lookup[[cluster_name]]
          ag_label <- ag_label_lookup[[ag_state]]
          raw_values <- acid_score[cells_in_cluster]
          cluster_summary_rows[[paste(key, cluster_name, sep = "__")]] <- data.frame(
            study = manifest_row$study,
            cancer_type_summary = manifest_row$cancer_type_summary,
            project_broad_compartment = manifest_row$project_broad_compartment,
            compartment_role = manifest_row$compartment_role,
            study_compartment_key = key,
            acid_raw_cluster = cluster_name,
            acid_raw_cluster_id = paste(key, cluster_name, sep = "__"),
            AG = ag_state,
            AG_Label = ag_label,
            AG_acid_rank = as.integer(rank_by_acid[[cluster_name]]),
            n_cells = length(cells_in_cluster),
            fraction_of_clusterable_cells = round_or_na(length(cells_in_cluster) / n_clusterable_cells),
            n_sample_units = length(unique(df$sample_unit_key[cells_in_cluster])),
            mean_acid_score = round_or_na(mean(raw_values, na.rm = TRUE)),
            median_acid_score = round_or_na(stats::median(raw_values, na.rm = TRUE)),
            min_acid_score = round_or_na(min(raw_values, na.rm = TRUE)),
            max_acid_score = round_or_na(max(raw_values, na.rm = TRUE)),
            stringsAsFactors = FALSE
          )
        }

        sample_tmp <- cell_out[cell_out$clusterable_for_05, c(
          "study", "cancer_type", "sample_unit_key", "project_broad_compartment",
          "study_compartment_key", "AG", "AG_Label"
        ), drop = FALSE]
        sample_tmp$n_cells <- 1L
        count_df <- stats::aggregate(
          n_cells ~ study + cancer_type + sample_unit_key + project_broad_compartment +
            study_compartment_key + AG + AG_Label,
          data = sample_tmp,
          FUN = sum
        )
        total_df <- stats::aggregate(
          n_cells ~ study + cancer_type + sample_unit_key + project_broad_compartment +
            study_compartment_key,
          data = sample_tmp,
          FUN = sum
        )
        colnames(total_df)[colnames(total_df) == "n_cells"] <- "n_clusterable_cells_in_sample_compartment"
        count_df <- merge(
          count_df,
          total_df,
          by = c("study", "cancer_type", "sample_unit_key", "project_broad_compartment", "study_compartment_key"),
          all.x = TRUE,
          sort = FALSE
        )
        count_df$fraction_of_sample_compartment <- round_or_na(
          count_df$n_cells / count_df$n_clusterable_cells_in_sample_compartment
        )
        sample_fraction_rows[[key]] <- count_df[order(
          count_df$sample_unit_key,
          count_df$AG
        ), , drop = FALSE]
      }
    }
  }

  if (cluster_status != "processed") {
    cell_out$cluster_status_05 <- cluster_status
    cell_out$clusterable_for_05 <- FALSE
  }

  manifest_rows[[key]] <- data.frame(
    study = manifest_row$study,
    cancer_type_summary = manifest_row$cancer_type_summary,
    project_broad_compartment = manifest_row$project_broad_compartment,
    compartment_role = manifest_row$compartment_role,
    study_compartment_key = key,
    n_cells = nrow(df),
    n_sample_units = manifest_row$n_sample_units,
    n_donors = manifest_row$n_donors,
    n_unique_cell_types_update = manifest_row$n_unique_cell_types_update,
    score_method = "AUCell",
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    auc_rank_fraction = auc_rank_fraction,
    auc_max_rank = auc_max_rank,
    acid_genes_total = length(acid_genes),
    acid_genes_present = length(acid_genes_present),
    acid_gene_coverage = round_or_na(length(acid_genes_present) / length(acid_genes)),
    score_status = score_status,
    n_scored_cells = n_scored_cells,
    k_clusters = k_clusters,
    cluster_status = cluster_status,
    cluster_note = cluster_note,
    n_clusterable_cells = n_clusterable_cells,
    min_cells_for_clustering = min_cells_for_clustering,
    min_present_acid_genes = min_present_acid_genes,
    kmeans_nstart = kmeans_nstart,
    kmeans_iter_max = kmeans_iter_max,
    random_seed_unit = random_seed + idx,
    kmeans_total_withinss = round_or_na(kmeans_total_withinss),
    kmeans_betweenss_fraction = round_or_na(kmeans_betweenss_fraction),
    cluster_size_min = cluster_size_min,
    cluster_size_max = cluster_size_max,
    acid_genes_present_list = collapse_or_none(acid_genes_present),
    acid_genes_missing_list = collapse_or_none(acid_genes_missing),
    stringsAsFactors = FALSE
  )

  cell_rows[[key]] <- cell_out
}

cell_level_ag_states <- do.call(rbind, cell_rows)
acid_gene_coverage <- do.call(rbind, coverage_rows)
acid_gene_detection <- do.call(rbind, gene_coverage_rows)
ag_state_manifest <- do.call(rbind, manifest_rows)
ag_cluster_summary <- if (length(cluster_summary_rows) > 0L) {
  do.call(rbind, cluster_summary_rows)
} else {
  data.frame()
}
sample_compartment_ag_fractions <- if (length(sample_fraction_rows) > 0L) {
  do.call(rbind, sample_fraction_rows)
} else {
  data.frame()
}

write_tsv(
  cell_level_ag_states,
  file.path(outdir, "cell_level_acid_scores_and_ag_states.tsv.gz"),
  compress = TRUE
)
write_tsv(
  acid_gene_coverage,
  file.path(outdir, "study_compartment_acid_gene_coverage.tsv")
)
write_tsv(
  acid_gene_detection,
  file.path(outdir, "study_compartment_acid_gene_detection_long.tsv")
)
write_tsv(
  ag_state_manifest,
  file.path(outdir, "study_compartment_ag_state_manifest.tsv")
)
write_tsv(
  ag_cluster_summary,
  file.path(outdir, "study_compartment_ag_cluster_summary.tsv")
)
write_tsv(
  sample_compartment_ag_fractions,
  file.path(outdir, "sample_compartment_ag_cell_fractions.tsv")
)

bundle_out <- list(
  parameters = list(
    input_rds = input_rds,
    input_refine_bundle = input_refine_bundle,
    acid_gene_file = acid_gene_file,
    outdir = outdir,
    score_assay = score_assay,
    score_slot_used = score_slot_used,
    score_method = "AUCell",
    auc_norm = auc_norm,
    auc_rank_fraction = auc_rank_fraction,
    k_clusters = k_clusters,
    kmeans_nstart = kmeans_nstart,
    kmeans_iter_max = kmeans_iter_max,
    random_seed = random_seed,
    min_cells_for_clustering = min_cells_for_clustering,
    min_present_acid_genes = min_present_acid_genes
  ),
  acid_genes = acid_genes,
  sample_metadata = sample_metadata,
  study_compartment_acid_gene_coverage = acid_gene_coverage,
  study_compartment_acid_gene_detection_long = acid_gene_detection,
  study_compartment_ag_state_manifest = ag_state_manifest,
  study_compartment_ag_cluster_summary = ag_cluster_summary,
  sample_compartment_ag_cell_fractions = sample_compartment_ag_fractions,
  cell_level_ag_states = cell_level_ag_states
)
saveRDS(bundle_out, file.path(outdir, "studywise_within_compartment_ag3_bundle.rds"))

run_summary_path <- file.path(outdir, "run_summary.txt")
summary_lines <- c(
  "05_studywise_within_compartment_acid_scoring_and_ag3_assignment completed.",
  paste0("Input Seurat object: ", input_rds),
  paste0("Input 03 refine bundle: ", input_refine_bundle),
  paste0("Acid gene file: ", acid_gene_file),
  paste0("Acid genes total: ", length(acid_genes)),
  paste0("Score assay: ", score_assay),
  paste0("Score slot used: ", score_slot_used),
  paste0("Score method: AUCell"),
  paste0("AUC rank fraction: ", auc_rank_fraction),
  paste0("K clusters: ", k_clusters),
  paste0("Min cells for clustering: ", min_cells_for_clustering),
  paste0("Studies processed: ", length(unique(ag_state_manifest$study))),
  paste0("Study x compartment units scored: ", sum(ag_state_manifest$score_status == "scored")),
  paste0("Study x compartment units clustered: ", sum(ag_state_manifest$cluster_status == "processed")),
  paste0("Study x compartment units skipped for AG clustering: ", sum(ag_state_manifest$cluster_status != "processed")),
  paste0("Cells scored acid: ", sum(is.finite(cell_level_ag_states$Acid_Score_05))),
  paste0("Cells assigned AG: ", sum(!is.na(cell_level_ag_states$AG))),
  paste0("Output cell-level table: ", file.path(outdir, "cell_level_acid_scores_and_ag_states.tsv.gz")),
  paste0("Output bundle: ", file.path(outdir, "studywise_within_compartment_ag3_bundle.rds"))
)
writeLines(summary_lines, con = run_summary_path)

message(paste(summary_lines, collapse = "\n"))
}

# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/03_map_sample_c1c4_states.R
stage_single_cell_03_map_sample_c1c4_states <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("cell_metadata", "sample_metadata", "ig_cell_states", "ag_cell_states")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# Source-preserved defaults; no study filter is applied.
options(stringsAsFactors = FALSE)
outdir <- out_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cell_metadata_path <- unname(paths["cell_metadata"])
sample_metadata_path <- unname(paths["sample_metadata"])
ig_cell_path <- unname(paths["ig_cell_states"])
ag_cell_path <- unname(paths["ag_cell_states"])
min_cells_for_standard_confidence <- 50L
low_margin_threshold <- 0.10
weak_top1_threshold <- 0.40
low_scorable_fraction_threshold <- 0.50

ag_levels <- c("AG1", "AG2", "AG3")
ig_levels <- c("IG1", "IG2")
agig_levels <- c("AG1_IG1", "AG1_IG2", "AG2_IG1", "AG2_IG2", "AG3_IG1", "AG3_IG2")
c_levels <- c("C1", "C2", "C3", "C4")

ag_label_map <- c(
  AG1 = "acid_high",
  AG2 = "acid_intermediate",
  AG3 = "acid_low"
)
ig_label_map <- c(
  IG1 = "lm12_inflammatory_remodeling_high",
  IG2 = "lm12_inflammatory_remodeling_low"
)
c_map <- c(
  AG1_IG1 = "C1",
  AG1_IG2 = "C2",
  AG2_IG1 = "C3",
  AG2_IG2 = "C2",
  AG3_IG1 = "C3",
  AG3_IG2 = "C4"
)
c_label_map <- c(
  C1 = "high_acid__high_inflammatory_remodeling",
  C2 = "high_or_intermediate_acid__low_inflammatory_remodeling",
  C3 = "intermediate_or_low_acid__high_inflammatory_remodeling",
  C4 = "low_acid__low_inflammatory_remodeling"
)
`write_tsv` <- bind_stage_helper(shared_write_tsv_006, environment())

read_tsv <- function(path, gz = grepl("\\.gz$", path)) {
  con <- if (gz) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, check.names = FALSE)
}
`safe_char` <- bind_stage_helper(shared_safe_char_042, environment())
`round_or_na` <- bind_stage_helper(shared_round_or_na_045, environment())
`collapse_or_none` <- bind_stage_helper(shared_collapse_or_none_043, environment())
`%notin%` <- bind_stage_helper(shared__notin__055, environment())

required_files <- c(cell_metadata_path, sample_metadata_path, ig_cell_path, ag_cell_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

assign_confidence_label <- function(n_cells, top1_prop, margin, tied_vote) {
  if (!is.finite(n_cells) || n_cells == 0L) {
    return("no_axis_cells")
  }
  if (isTRUE(tied_vote)) {
    return("tied_vote")
  }
  if (n_cells < min_cells_for_standard_confidence) {
    return("low_cell_count")
  }
  low_margin <- is.finite(margin) && margin < low_margin_threshold
  weak_top1 <- is.finite(top1_prop) && top1_prop < weak_top1_threshold
  if (low_margin && weak_top1) {
    return("low_margin_and_weak_top1")
  }
  if (low_margin) {
    return("low_margin")
  }
  if (weak_top1) {
    return("weak_top1")
  }
  "standard_confidence"
}

combined_confidence_label <- function(ag_conf, ig_conf, n_both, scorable_fraction) {
  if (!is.finite(n_both) || n_both == 0L) {
    return("no_joint_ag_ig_cells")
  }
  if (!is.finite(scorable_fraction) || scorable_fraction < low_scorable_fraction_threshold) {
    return("low_scorable_non_other_fraction")
  }
  confs <- c(ag_conf, ig_conf)
  if (any(confs %in% c("no_axis_cells", "tied_vote", "low_cell_count"))) {
    return("low_axis_confidence")
  }
  if (any(confs %in% c("low_margin_and_weak_top1", "low_margin", "weak_top1"))) {
    return("moderate_axis_confidence")
  }
  "standard_confidence"
}

summarise_axis_votes <- function(cell_df, state_col, state_levels, axis_name) {
  cell_df <- cell_df[cell_df[[state_col]] %in% state_levels, , drop = FALSE]
  sample_keys <- sort(unique(safe_char(cell_df$sample_unit_key)))
  if (length(sample_keys) == 0L) {
    return(list(wide = data.frame(), long = data.frame()))
  }

  tab <- table(
    factor(cell_df$sample_unit_key, levels = sample_keys),
    factor(cell_df[[state_col]], levels = state_levels)
  )
  counts <- as.data.frame.matrix(tab)
  colnames(counts) <- paste0("n_", state_levels)
  counts$sample_unit_key <- rownames(counts)
  rownames(counts) <- NULL
  counts <- counts[, c("sample_unit_key", paste0("n_", state_levels)), drop = FALSE]
  n_axis_cells <- rowSums(counts[, paste0("n_", state_levels), drop = FALSE])

  wide_rows <- lapply(seq_len(nrow(counts)), function(i) {
    vals <- as.numeric(counts[i, paste0("n_", state_levels)])
    names(vals) <- state_levels
    n_total <- sum(vals)
    props <- if (n_total > 0L) vals / n_total else rep(NA_real_, length(vals))
    names(props) <- state_levels
    order_idx <- order(-props, match(names(props), state_levels), na.last = TRUE)
    ordered_states <- names(props)[order_idx]
    top1_state <- if (n_total > 0L) ordered_states[1] else NA_character_
    top2_state <- if (n_total > 0L && length(ordered_states) >= 2L) ordered_states[2] else NA_character_
    top1_prop <- if (!is.na(top1_state)) props[top1_state] else NA_real_
    top2_prop <- if (!is.na(top2_state)) props[top2_state] else NA_real_
    margin <- top1_prop - top2_prop
    tied_vote <- if (n_total > 0L) {
      sum(abs(props - max(props, na.rm = TRUE)) < .Machine$double.eps^0.5, na.rm = TRUE) > 1L
    } else {
      FALSE
    }

    count_part <- as.list(vals)
    names(count_part) <- paste0("sample_", axis_name, "_n_", state_levels)
    prop_part <- as.list(round_or_na(props))
    names(prop_part) <- paste0("sample_", axis_name, "_prop_", state_levels)

    data.frame(
      sample_unit_key = counts$sample_unit_key[i],
      sample_axis = axis_name,
      n_axis_cells = n_total,
      dominant_state = top1_state,
      top1_fraction = round_or_na(top1_prop),
      top2_state = top2_state,
      top2_fraction = round_or_na(top2_prop),
      dominance_margin = round_or_na(margin),
      tied_vote = tied_vote,
      confidence_label = assign_confidence_label(n_total, top1_prop, margin, tied_vote),
      as.data.frame(count_part, check.names = FALSE),
      as.data.frame(prop_part, check.names = FALSE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  wide <- do.call(rbind, wide_rows)

  long_rows <- list()
  for (i in seq_len(nrow(counts))) {
    n_total <- n_axis_cells[i]
    for (state in state_levels) {
      n_state <- counts[i, paste0("n_", state)]
      long_rows[[paste(counts$sample_unit_key[i], axis_name, state, sep = "__")]] <- data.frame(
        sample_unit_key = counts$sample_unit_key[i],
        sample_axis = axis_name,
        state = state,
        n_cells = n_state,
        n_axis_cells = n_total,
        fraction = if (n_total == 0L) NA_real_ else round_or_na(n_state / n_total),
        stringsAsFactors = FALSE
      )
    }
  }
  long <- do.call(rbind, long_rows)
  list(wide = wide, long = long)
}

message("Loading cell/sample metadata and AG/IG cell-level states.")
cell_metadata <- read_tsv(cell_metadata_path)
sample_metadata <- read_tsv(sample_metadata_path)
ig_cells <- read_tsv(ig_cell_path)
ag_cells <- read_tsv(ag_cell_path)

required_meta_cols <- c("cell_id", "sample_unit_key", "study", "cancer_type", "project_broad_compartment")
missing_meta_cols <- required_meta_cols[required_meta_cols %notin% colnames(cell_metadata)]
if (length(missing_meta_cols) > 0L) {
  stop("analysis_ready_cell_metadata is missing required columns: ", paste(missing_meta_cols, collapse = ", "))
}
if (!"sample_unit_key" %in% colnames(sample_metadata)) {
  stop("analysis_ready_sample_metadata is missing sample_unit_key.")
}
if (!all(c("cell_id", "sample_unit_key", "AG", "AG_Label", "Acid_Score_05") %in% colnames(ag_cells))) {
  stop("AG cell table is missing required columns.")
}
if (!all(c("cell_id", "sample_unit_key", "IG", "IG_Label", "lm12_inflammatory_axis_04") %in% colnames(ig_cells))) {
  stop("IG cell table is missing required columns.")
}

non_other_meta <- cell_metadata[cell_metadata$project_broad_compartment != "Other", , drop = FALSE]
other_meta <- cell_metadata[cell_metadata$project_broad_compartment == "Other", , drop = FALSE]

if (length(setdiff(ag_cells$cell_id, non_other_meta$cell_id)) > 0L) {
  stop("Some AG cells are not present in the Step 01 non-Other metadata.")
}
if (length(setdiff(ig_cells$cell_id, non_other_meta$cell_id)) > 0L) {
  stop("Some IG cells are not present in the Step 01 non-Other metadata.")
}
if (length(setdiff(ag_cells$cell_id, ig_cells$cell_id)) + length(setdiff(ig_cells$cell_id, ag_cells$cell_id)) > 0L) {
  stop("AG and IG cell tables do not contain the same cell_id set.")
}

ag_keep_cols <- intersect(c(
  "cell_id", "AG", "AG_Label", "AG_acid_rank", "Acid_Score_05",
  "cluster_status_05", "clusterable_for_05"
), colnames(ag_cells))
ig_keep_cols <- intersect(c(
  "cell_id", "IG", "IG_Label", "IG_axis_rank", "lm12_inflammatory_axis_04",
  "cluster_status_04", "clusterable_for_04"
), colnames(ig_cells))
meta_keep_cols <- intersect(c(
  "cell_id", "study", "cancer_type", "donor_id", "Sample", "donor_id_pre_post",
  "sample_unit_key", "pre_post", "outcome", "Combined_outcome", "treatment_raw",
  "treatment_class_broad", "cell_types", "cell_types_update",
  "Cell_type_broad_original", "project_broad_compartment", "nCount_RNA", "nFeature_RNA"
), colnames(non_other_meta))

cell_join <- merge(
  non_other_meta[, meta_keep_cols, drop = FALSE],
  ag_cells[, ag_keep_cols, drop = FALSE],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)
cell_join <- merge(
  cell_join,
  ig_cells[, ig_keep_cols, drop = FALSE],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)
cell_join$cell_has_AG <- cell_join$AG %in% ag_levels
cell_join$cell_has_IG <- cell_join$IG %in% ig_levels
cell_join$cell_has_joint_AGIG <- cell_join$cell_has_AG & cell_join$cell_has_IG
cell_join$cell_AGIG <- ifelse(cell_join$cell_has_joint_AGIG, paste(cell_join$AG, cell_join$IG, sep = "_"), NA_character_)

sample_keys_all <- safe_char(sample_metadata$sample_unit_key)
sample_keys_all <- sample_keys_all[sample_keys_all != "" & sample_keys_all != "NA"]

sample_meta_unique <- sample_metadata[!duplicated(sample_metadata$sample_unit_key), , drop = FALSE]

total_counts <- as.data.frame(table(factor(cell_metadata$sample_unit_key, levels = sample_keys_all)))
colnames(total_counts) <- c("sample_unit_key", "n_cells_total_from_cell_metadata")
non_other_counts <- as.data.frame(table(factor(non_other_meta$sample_unit_key, levels = sample_keys_all)))
colnames(non_other_counts) <- c("sample_unit_key", "n_non_other_cells")
other_counts <- as.data.frame(table(factor(other_meta$sample_unit_key, levels = sample_keys_all)))
colnames(other_counts) <- c("sample_unit_key", "n_other_cells")

ag_summary <- summarise_axis_votes(cell_join[cell_join$cell_has_AG, , drop = FALSE], "AG", ag_levels, "AG")
ig_summary <- summarise_axis_votes(cell_join[cell_join$cell_has_IG, , drop = FALSE], "IG", ig_levels, "IG")

sample_level <- merge(sample_meta_unique, total_counts, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
sample_level <- merge(sample_level, non_other_counts, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
sample_level <- merge(sample_level, other_counts, by = "sample_unit_key", all.x = TRUE, sort = FALSE)

ag_wide <- ag_summary$wide
ig_wide <- ig_summary$wide
if (nrow(ag_wide) > 0L) {
  names(ag_wide) <- ifelse(names(ag_wide) == "sample_unit_key", "sample_unit_key", paste0("AG_axis__", names(ag_wide)))
  sample_level <- merge(sample_level, ag_wide, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
}
if (nrow(ig_wide) > 0L) {
  names(ig_wide) <- ifelse(names(ig_wide) == "sample_unit_key", "sample_unit_key", paste0("IG_axis__", names(ig_wide)))
  sample_level <- merge(sample_level, ig_wide, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
}

for (col in c("n_cells_total_from_cell_metadata", "n_non_other_cells", "n_other_cells")) {
  sample_level[[col]][is.na(sample_level[[col]])] <- 0L
}
sample_level$other_fraction <- ifelse(
  sample_level$n_cells_total_from_cell_metadata > 0,
  round_or_na(sample_level$n_other_cells / sample_level$n_cells_total_from_cell_metadata),
  NA_real_
)
sample_level$non_other_fraction <- ifelse(
  sample_level$n_cells_total_from_cell_metadata > 0,
  round_or_na(sample_level$n_non_other_cells / sample_level$n_cells_total_from_cell_metadata),
  NA_real_
)

joint_counts <- stats::aggregate(
  cell_has_joint_AGIG ~ sample_unit_key,
  data = cell_join,
  FUN = sum
)
colnames(joint_counts)[2] <- "n_joint_ag_ig_cells"
sample_level <- merge(sample_level, joint_counts, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
sample_level$n_joint_ag_ig_cells[is.na(sample_level$n_joint_ag_ig_cells)] <- 0L
sample_level$scorable_non_other_fraction <- ifelse(
  sample_level$n_non_other_cells > 0,
  round_or_na(sample_level$n_joint_ag_ig_cells / sample_level$n_non_other_cells),
  NA_real_
)

compartment_summary <- stats::aggregate(
  project_broad_compartment ~ sample_unit_key,
  data = cell_join[cell_join$cell_has_joint_AGIG, , drop = FALSE],
  FUN = function(x) length(unique(x))
)
colnames(compartment_summary)[2] <- "n_contributing_project_broad_compartments"
compartment_list <- stats::aggregate(
  project_broad_compartment ~ sample_unit_key,
  data = cell_join[cell_join$cell_has_joint_AGIG, , drop = FALSE],
  FUN = collapse_or_none
)
colnames(compartment_list)[2] <- "contributing_project_broad_compartments"
sample_level <- merge(sample_level, compartment_summary, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
sample_level <- merge(sample_level, compartment_list, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
sample_level$n_contributing_project_broad_compartments[is.na(sample_level$n_contributing_project_broad_compartments)] <- 0L
sample_level$contributing_project_broad_compartments[is.na(sample_level$contributing_project_broad_compartments)] <- "none"

sample_level$sample_AG <- safe_char(sample_level$AG_axis__dominant_state)
sample_level$sample_AG[!sample_level$sample_AG %in% ag_levels] <- NA_character_
sample_level$sample_IG <- safe_char(sample_level$IG_axis__dominant_state)
sample_level$sample_IG[!sample_level$sample_IG %in% ig_levels] <- NA_character_
sample_level$sample_AG_Label <- unname(ag_label_map[sample_level$sample_AG])
sample_level$sample_IG_Label <- unname(ig_label_map[sample_level$sample_IG])
sample_level$sample_AGIG <- ifelse(
  !is.na(sample_level$sample_AG) & !is.na(sample_level$sample_IG),
  paste(sample_level$sample_AG, sample_level$sample_IG, sep = "_"),
  NA_character_
)
sample_level$sample_C_subtype <- unname(c_map[sample_level$sample_AGIG])
sample_level$sample_C_label <- unname(c_label_map[sample_level$sample_C_subtype])

sample_level$sample_C_confidence_label <- mapply(
  combined_confidence_label,
  ag_conf = sample_level$AG_axis__confidence_label,
  ig_conf = sample_level$IG_axis__confidence_label,
  n_both = sample_level$n_joint_ag_ig_cells,
  scorable_fraction = sample_level$scorable_non_other_fraction,
  USE.NAMES = FALSE
)
sample_level$low_scorable_non_other_fraction_flag <- is.finite(sample_level$scorable_non_other_fraction) &
  sample_level$scorable_non_other_fraction < low_scorable_fraction_threshold

axis_long <- rbind(
  ag_summary$long,
  ig_summary$long
)

atomic_table <- table(
  factor(cell_join$sample_unit_key[cell_join$cell_AGIG %in% agig_levels], levels = sample_keys_all),
  factor(cell_join$cell_AGIG[cell_join$cell_AGIG %in% agig_levels], levels = agig_levels)
)
atomic_counts <- as.data.frame.matrix(atomic_table)
atomic_counts$sample_unit_key <- rownames(atomic_counts)
rownames(atomic_counts) <- NULL
atomic_counts <- atomic_counts[, c("sample_unit_key", agig_levels), drop = FALSE]
colnames(atomic_counts) <- c("sample_unit_key", paste0("n_cell_", agig_levels))
atomic_counts$n_joint_ag_ig_cells <- rowSums(atomic_counts[, paste0("n_cell_", agig_levels), drop = FALSE])
for (lvl in agig_levels) {
  n_col <- paste0("n_cell_", lvl)
  p_col <- paste0("prop_cell_", lvl)
  atomic_counts[[p_col]] <- ifelse(
    atomic_counts$n_joint_ag_ig_cells > 0,
    round_or_na(atomic_counts[[n_col]] / atomic_counts$n_joint_ag_ig_cells),
    NA_real_
  )
}

count_by_group <- function(df, group_col, state_col, state_levels, out_name) {
  df <- df[!is.na(df[[state_col]]) & df[[state_col]] %in% state_levels, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(data.frame())
  }
  group_values <- safe_char(df[[group_col]])
  group_values[group_values == "" | group_values == "NA"] <- "NA"
  split_df <- split(df, group_values, drop = TRUE)
  rows <- list()
  for (key in names(split_df)) {
    sub <- split_df[[key]]
    total <- nrow(sub)
    for (state in state_levels) {
      rows[[paste(key, state, sep = "__")]] <- data.frame(
        summary_type = out_name,
        grouping = group_col,
        group_value = key,
        state = state,
        n_samples = sum(sub[[state_col]] == state, na.rm = TRUE),
        total_samples = total,
        fraction = round_or_na(sum(sub[[state_col]] == state, na.rm = TRUE) / total),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

sample_level$sample_AG <- factor(sample_level$sample_AG, levels = ag_levels)
sample_level$sample_IG <- factor(sample_level$sample_IG, levels = ig_levels)
sample_level$sample_AGIG <- factor(sample_level$sample_AGIG, levels = agig_levels)
sample_level$sample_C_subtype <- factor(sample_level$sample_C_subtype, levels = c_levels)

sample_count_tables <- rbind(
  count_by_group(sample_level, "study", "sample_AG", ag_levels, "AG_by_study"),
  count_by_group(sample_level, "study", "sample_IG", ig_levels, "IG_by_study"),
  count_by_group(sample_level, "study", "sample_AGIG", agig_levels, "AGIG_by_study"),
  count_by_group(sample_level, "study", "sample_C_subtype", c_levels, "C_by_study"),
  count_by_group(sample_level, "cancer_type", "sample_C_subtype", c_levels, "C_by_cancer_type"),
  count_by_group(sample_level, "pre_post", "sample_C_subtype", c_levels, "C_by_pre_post"),
  count_by_group(sample_level, "Combined_outcome", "sample_C_subtype", c_levels, "C_by_combined_outcome"),
  count_by_group(sample_level, "treatment_class_broad", "sample_C_subtype", c_levels, "C_by_treatment_class")
)

sample_level$sample_AG <- as.character(sample_level$sample_AG)
sample_level$sample_IG <- as.character(sample_level$sample_IG)
sample_level$sample_AGIG <- as.character(sample_level$sample_AGIG)
sample_level$sample_C_subtype <- as.character(sample_level$sample_C_subtype)

minimal_cols <- intersect(c(
  "sample_unit_key", "donor_id_pre_post", "donor_id", "study", "cancer_type",
  "pre_post", "outcome", "Combined_outcome", "treatment_class_broad",
  "n_cells_total", "n_cells_total_from_cell_metadata", "n_non_other_cells",
  "n_other_cells", "other_fraction", "n_joint_ag_ig_cells",
  "scorable_non_other_fraction", "n_contributing_project_broad_compartments",
  "sample_AG", "sample_AG_Label", "sample_IG", "sample_IG_Label",
  "sample_AGIG", "sample_C_subtype", "sample_C_label",
  "sample_C_confidence_label",
  "AG_axis__top1_fraction", "AG_axis__dominance_margin",
  "IG_axis__top1_fraction", "IG_axis__dominance_margin"
), colnames(sample_level))

write_tsv(
  cell_join,
  file.path(outdir, "cell_level_ag_ig_joined_states.tsv.gz"),
  compress = TRUE
)
write_tsv(
  sample_level,
  file.path(outdir, "sample_level_ag_ig_c1c4_labels.tsv")
)
write_tsv(
  sample_level[, minimal_cols, drop = FALSE],
  file.path(outdir, "sample_level_ag_ig_c1c4_labels_minimal.tsv")
)
write_tsv(
  axis_long,
  file.path(outdir, "sample_axis_state_fractions_long.tsv")
)
write_tsv(
  atomic_counts,
  file.path(outdir, "sample_cell_level_agig_atomic_fractions.tsv")
)
write_tsv(
  sample_count_tables,
  file.path(outdir, "sample_level_group_count_summary.tsv")
)

bundle_out <- list(
  parameters = list(
    cell_metadata_path = cell_metadata_path,
    sample_metadata_path = sample_metadata_path,
    ig_cell_path = ig_cell_path,
    ag_cell_path = ag_cell_path,
    outdir = outdir,
    ag_levels = ag_levels,
    ig_levels = ig_levels,
    agig_levels = agig_levels,
    c_levels = c_levels,
    c_map = c_map,
    min_cells_for_standard_confidence = min_cells_for_standard_confidence,
    low_margin_threshold = low_margin_threshold,
    weak_top1_threshold = weak_top1_threshold,
    low_scorable_fraction_threshold = low_scorable_fraction_threshold
  ),
  sample_level_ag_ig_c1c4_labels = sample_level,
  sample_level_group_count_summary = sample_count_tables,
  sample_axis_state_fractions_long = axis_long,
  sample_cell_level_agig_atomic_fractions = atomic_counts,
  cell_level_ag_ig_joined_states = cell_join
)
saveRDS(bundle_out, file.path(outdir, "sample_level_ag_ig_c1c4_mapping_bundle.rds"))

run_summary_path <- file.path(outdir, "run_summary.txt")
n_samples_total <- nrow(sample_level)
n_samples_c <- sum(!is.na(sample_level$sample_C_subtype))
summary_lines <- c(
  "06_sample_level_ag_ig_mapping_and_c1c4_recoding completed.",
  paste0("Input AG cell table: ", ag_cell_path),
  paste0("Input IG cell table: ", ig_cell_path),
  paste0("Samples total: ", n_samples_total),
  paste0("Samples with sample_AG: ", sum(!is.na(sample_level$sample_AG))),
  paste0("Samples with sample_IG: ", sum(!is.na(sample_level$sample_IG))),
  paste0("Samples with sample_C_subtype: ", n_samples_c),
  paste0("Cells with joint AG/IG: ", sum(cell_join$cell_has_joint_AGIG)),
  paste0("C1 samples: ", sum(sample_level$sample_C_subtype == "C1", na.rm = TRUE)),
  paste0("C2 samples: ", sum(sample_level$sample_C_subtype == "C2", na.rm = TRUE)),
  paste0("C3 samples: ", sum(sample_level$sample_C_subtype == "C3", na.rm = TRUE)),
  paste0("C4 samples: ", sum(sample_level$sample_C_subtype == "C4", na.rm = TRUE)),
  paste0("Output sample labels: ", file.path(outdir, "sample_level_ag_ig_c1c4_labels.tsv")),
  paste0("Output bundle: ", file.path(outdir, "sample_level_ag_ig_c1c4_mapping_bundle.rds"))
)
writeLines(summary_lines, con = run_summary_path)

message(paste(summary_lines, collapse = "\n"))
}
