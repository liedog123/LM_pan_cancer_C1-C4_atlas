#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 02_score_acid_assign_ag_states.R --config PATH [--check-config]\n       Rscript 02_score_acid_assign_ag_states.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("seurat_object", "refine_bundle", "selected_genes", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "refine_bundle", "selected_genes")]
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

collapse_or_none <- function(x) {
  x <- unique(safe_char(x))
  x <- x[x != "" & x != "NA"]
  if (length(x) == 0L) {
    return("none")
  }
  paste(sort(x), collapse = ",")
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

  # AUCell 只依赖单细胞内的表达排序，这里与 02/03 保持一致，优先使用 counts。
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
