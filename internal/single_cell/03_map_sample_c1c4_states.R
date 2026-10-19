#!/usr/bin/env Rscript
# Purpose: Score fixed gene sets and construct cellular and sample-level state labels.
# Task: single_cell_states; step: s03
# Inputs (DCF fields): cell_metadata,sample_metadata,ig_cell_states,ag_cell_states
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_states --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 03_map_sample_c1c4_states.R --config PATH [--check-config]\n       Rscript 03_map_sample_c1c4_states.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("cell_metadata", "sample_metadata", "ig_cell_states", "ag_cell_states", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("cell_metadata", "sample_metadata", "ig_cell_states", "ag_cell_states")]
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

write_tsv <- function(df, path, compress = FALSE) {
  con <- if (compress) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

read_tsv <- function(path, gz = grepl("\\.gz$", path)) {
  con <- if (gz) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, check.names = FALSE)
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

`%notin%` <- function(x, table) !(x %in% table)

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
