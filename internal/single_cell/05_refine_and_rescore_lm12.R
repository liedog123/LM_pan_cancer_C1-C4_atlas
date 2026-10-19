#!/usr/bin/env Rscript
# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Task: single_cell_states; step: s05
# Inputs (DCF fields): seurat_object,initial_bundle
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_states --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 05_refine_and_rescore_lm12.R --config PATH [--check-config]\n       Rscript 05_refine_and_rescore_lm12.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("seurat_object", "initial_bundle", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "initial_bundle")]
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

collapse_or_none <- function(x, n_max = 20L) {
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

  # AUCell 只依赖单细胞内的表达排序，这里优先沿用 counts，避免依赖不稳定的 data 层。
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
