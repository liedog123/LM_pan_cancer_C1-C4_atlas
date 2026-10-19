#!/usr/bin/env Rscript
# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Task: single_cell_states; step: s01
# Inputs (DCF fields): refine_bundle
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_states --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 01_assign_study_compartment_ig_states.R --config PATH [--check-config]\n       Rscript 01_assign_study_compartment_ig_states.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("refine_bundle", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("refine_bundle")]
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
