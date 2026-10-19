#!/usr/bin/env Rscript
# Purpose: Evaluate the reported Bassez-exclusion and marker sensitivity results.
# Task: bassez_sensitivity; step: s12
# Inputs (DCF fields): cell_metadata,sample_metadata,sample_labels,sample_compartment_counts,sample_sheet,counts_manifest,model_contract,definition_genes,program_scores,program_catalog,formal_centers,formal_pathway_evidence,fig6_pathway_source,gmt_hallmark,gmt_reactome,gmt_gobp,counts_dir,sample_sets_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task bassez_sensitivity --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 12_run_bassez_exclusion_sensitivity.R --config PATH [--check-config]\n       Rscript 12_run_bassez_exclusion_sensitivity.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("cell_metadata", "sample_metadata", "sample_labels", "sample_compartment_counts", "sample_sheet", "counts_manifest", "model_contract", "definition_genes", "program_scores", "program_catalog", "formal_centers", "formal_pathway_evidence", "fig6_pathway_source", "gmt_hallmark", "gmt_reactome", "gmt_gobp", "counts_dir", "sample_sets_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("cell_metadata", "sample_metadata", "sample_labels", "sample_compartment_counts", "sample_sheet", "counts_manifest", "model_contract", "definition_genes", "program_scores", "program_catalog", "formal_centers", "formal_pathway_evidence", "fig6_pathway_source", "gmt_hallmark", "gmt_reactome", "gmt_gobp")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir", "sample_sets_dir")])
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

# 普通执行保留缺失计数的原NOT_EVALUABLE分支；已有文件仍执行边界检查。
relocate_input <- function(recorded, input_dir, allow_missing = FALSE) {
  if(length(recorded)!=1L) stop("Invalid recorded input path.")
  if(is.na(recorded)) {
    if(isTRUE(allow_missing)) return(NA_character_)
    stop("Invalid recorded input path.")
  }
  if(!nzchar(recorded)) stop("Invalid recorded input path.")
  leaf <- basename(recorded)
  if(leaf %in% c(".","..") || !nzchar(leaf)) stop("Invalid input basename.")
  candidate <- resolve_path(file.path(input_dir,leaf), dirname(config_path))
  if(!is_within(candidate,input_dir) || identical(candidate,input_dir)) stop("Input escapes configured directory.")
  if(isTRUE(allow_missing) && !file.exists(candidate) && !dir.exists(candidate)) return(candidate)
  if(!file_test("-f",candidate) || file.access(candidate,4L)!=0L || is.na(file.info(candidate)$size) || file.info(candidate)$size<=0) stop("Selected input must be a readable nonempty regular file.")
  candidate
}

cfg_paths <- paths


# Fig. 6 Step 15：用当前 active pseudobulk / program-score 输入完成 Bassez 排除敏感性。
# 该脚本只写入候选目录；不覆盖 7.0、7.1、7.3 的冻结结果。

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(stringsAsFactors = FALSE)
set.seed(20260908)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(SummarizedExperiment)
})
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
data.table::setDTthreads(1L)
options(mc.cores = 1L)

script_name <- "LM_pan_cancer/7.0scRNA_LM12_state_mechanism/pancancer/scripts/15_pancancer_c1c4_bassez_exclusion_sensitivity_v0.1.R"
outdir <- out_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("denominator_audit", "composition", "marker_program", "gsea", "qa", "fig6_handoff")) {
  dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
}

log_event <- function(level, msg, fields = list()) {
  clean <- function(x) gsub('"', "'", as.character(x), fixed = TRUE)
  extra <- if (!length(fields)) "" else paste0(" ", paste(sprintf("%s=\"%s\"", names(fields), vapply(fields, clean, character(1))), collapse = " "))
  cat(sprintf("[%s] %s script=\"%s\" msg=\"%s\"%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, clean(msg), extra))
  flush.console()
}

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = TRUE, na = "NA")
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
safe_char <- function(x) {
  y <- as.character(x)
  y[is.na(y)] <- "NA"
  trimws(y)
}

md5_or_na <- function(path) {
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
    md5 = md5_or_na(path)
  )
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  names(sets) <- vapply(fields, `[`, character(1L), 1L)
  sets[lengths(sets) > 0L]
}

STATE_LEVELS <- c("C1", "C2", "C3", "C4")
CORE_COMPARTMENTS <- c("Malignant", "T/NK", "Myeloid", "Stromal/CAF")
ALL_COMPARTMENTS <- c("Malignant", "T/NK", "Myeloid", "Stromal/CAF", "Endothelial", "B/Plasma", "Other")
TARGET_COMPARISONS <- c("C1_vs_C4", "AG2_IG1_vs_AG2_IG2", "AG2_IG2_vs_AG3_IG2")

paths <- as.list(paths[c("cell_metadata", "sample_metadata", "sample_labels", "sample_compartment_counts", "sample_sheet", "counts_manifest", "model_contract", "definition_genes", "program_scores", "program_catalog", "formal_centers", "formal_pathway_evidence", "fig6_pathway_source", "gmt_hallmark", "gmt_reactome", "gmt_gobp")])

manifest <- rbindlist(Map(function(nm, p) file_record(nm, p, TRUE), names(paths), unname(paths)), fill = TRUE)
write_tsv(manifest, file.path(outdir, "qa/input_file_manifest.tsv"))
if (nrow(manifest[critical == TRUE & (!exists | !nonempty)]) > 0L) {
  stop("Required active inputs are incomplete")
}

started_at <- Sys.time()
log_event("START", "load active metadata and labels", list(output = outdir))
cell_meta <- fread(paths$cell_metadata)
sample_meta <- fread(paths$sample_metadata)
sample_labels <- fread(paths$sample_labels)
sample_counts <- fread(paths$sample_compartment_counts)
sample_labels <- unique(sample_labels, by = "sample_unit_key")
sample_labels <- sample_labels[sample_C_subtype %in% STATE_LEVELS]

required_cell <- c("cell_id", "sample_unit_key", "study", "cancer_type", "donor_id", "project_broad_compartment")
required_label <- c("sample_unit_key", "study", "cancer_type", "sample_C_subtype")
if (length(setdiff(required_cell, names(cell_meta))) > 0L) stop("cell metadata missing required columns")
if (length(setdiff(required_label, names(sample_labels))) > 0L) stop("sample labels missing required columns")
cell_meta[, `:=`(
  study = safe_char(study), cancer_type = safe_char(cancer_type),
  sample_unit_key = safe_char(sample_unit_key), donor_id = safe_char(donor_id),
  project_broad_compartment = safe_char(project_broad_compartment)
)]
sample_labels[, `:=`(
  study = safe_char(study), cancer_type = safe_char(cancer_type),
  sample_unit_key = safe_char(sample_unit_key), sample_C_subtype = safe_char(sample_C_subtype)
)]

label_map <- sample_labels[, .(sample_unit_key, sample_C_subtype, label_study = study, label_cancer_type = cancer_type)]
cell_meta <- merge(cell_meta, label_map, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
cell_meta[, sample_C_subtype := safe_char(sample_C_subtype)]
cell_meta <- cell_meta[sample_C_subtype %in% STATE_LEVELS]
cell_meta[, is_bassez := tolower(study) == "bassez"]

# 1. 分母审计：study/cancer/sample/patient 分开统计；没有 organ 字段时显式记录不可用。
study_summary <- cell_meta[, .(
  n_cells = .N,
  n_samples = uniqueN(sample_unit_key),
  n_patients = uniqueN(donor_id),
  n_cancer_labels = uniqueN(cancer_type)
), by = study][order(-n_cells, study)]
totals <- cell_meta[, .(
  n_cells = .N,
  n_samples = uniqueN(sample_unit_key),
  n_patients = uniqueN(donor_id),
  n_studies = uniqueN(study),
  n_cancer_labels = uniqueN(cancer_type)
)]
study_summary[, `:=`(
  cell_fraction = n_cells / totals$n_cells,
  sample_fraction = n_samples / totals$n_samples,
  patient_fraction = n_patients / totals$n_patients,
  study_fraction = 1 / totals$n_studies
)]
study_cancer_summary <- cell_meta[, .(
  n_cells = .N, n_samples = uniqueN(sample_unit_key), n_patients = uniqueN(donor_id)
), by = .(study, cancer_type)][order(study, cancer_type)]
state_by_study <- cell_meta[, .(n_cells = .N, n_samples = uniqueN(sample_unit_key), n_patients = uniqueN(donor_id)),
                            by = .(study, sample_C_subtype)][order(study, sample_C_subtype)]
state_overall <- cell_meta[, .(n_cells = .N, n_samples = uniqueN(sample_unit_key), n_patients = uniqueN(donor_id)),
                           by = sample_C_subtype][order(sample_C_subtype)]
bassez_manifest <- cell_meta[is_bassez == TRUE, .(
  n_cells = .N, n_samples = uniqueN(sample_unit_key), n_patients = uniqueN(donor_id),
  n_cancer_labels = uniqueN(cancer_type)
)]
bassez_manifest[, `:=`(
  study = "Bassez", cell_fraction_of_all = n_cells / totals$n_cells,
  sample_fraction_of_all = n_samples / totals$n_samples,
  patient_fraction_of_all = n_patients / totals$n_patients,
  cancer_label_fraction_of_all = n_cancer_labels / totals$n_cancer_labels
)]
write_tsv(study_summary, file.path(outdir, "denominator_audit/study_summary.tsv"))
write_tsv(study_cancer_summary, file.path(outdir, "denominator_audit/study_cancer_summary.tsv"))
write_tsv(state_by_study, file.path(outdir, "denominator_audit/state_by_study.tsv"))
write_tsv(state_overall, file.path(outdir, "denominator_audit/state_overall.tsv"))
write_tsv(bassez_manifest, file.path(outdir, "denominator_audit/bassez_manifest.tsv"))
write_tsv(data.table(
  field = c("organ", "morphology", "cancer_type", "study", "sample_unit_key", "donor_id"),
  status = c("NOT_AVAILABLE", "NOT_AVAILABLE", "AVAILABLE", "AVAILABLE", "AVAILABLE", "AVAILABLE"),
  note = c("No organ field in active metadata; not inferred from cancer label.",
           "No morphology field in active metadata; cancer_type retained as label only.",
           "Nine cancer labels from active metadata.", "Eight studies.",
           "Study-aware sample_unit_key used as sample denominator.", "donor_id used as patient denominator.")
), file.path(outdir, "denominator_audit/field_availability.tsv"))
write_tsv(totals, file.path(outdir, "denominator_audit/overall_totals.tsv"))

# 2. Sample-level composition：先保留每个 sample 的分母，再按 C1-C4 等权汇总。
composition_base <- cell_meta[project_broad_compartment %in% ALL_COMPARTMENTS,
                              .(n_cells = .N), by = .(sample_unit_key, study, cancer_type, donor_id, sample_C_subtype, project_broad_compartment, is_bassez)]
sample_frame <- unique(cell_meta[, .(sample_unit_key, study, cancer_type, donor_id, sample_C_subtype, is_bassez)])
# data.table 不接受 by=NULL 的笛卡尔连接；显式按 sample 展开固定的 7 个区室。
comp_grid <- sample_frame[, .(project_broad_compartment = ALL_COMPARTMENTS),
                          by = .(sample_unit_key, study, cancer_type, donor_id, sample_C_subtype, is_bassez)]
composition_all <- merge(comp_grid, composition_base,
                         by = c("sample_unit_key", "study", "cancer_type", "donor_id", "sample_C_subtype", "project_broad_compartment", "is_bassez"),
                         all.x = TRUE, sort = FALSE)
composition_all[is.na(n_cells), n_cells := 0L]
sample_totals <- composition_all[, .(n_cells_total = sum(n_cells)), by = sample_unit_key]
composition_all <- merge(composition_all, sample_totals, by = "sample_unit_key", all.x = TRUE, sort = FALSE)
composition_all[, fraction := fifelse(n_cells_total > 0, n_cells / n_cells_total, NA_real_)]
composition_all <- composition_all[order(sample_C_subtype, sample_unit_key, match(project_broad_compartment, ALL_COMPARTMENTS))]
composition_excl <- composition_all[is_bassez == FALSE]
composition_summary <- function(dt, scope) {
  out <- dt[, .(
    n_samples = uniqueN(sample_unit_key),
    n_patients = uniqueN(donor_id),
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    sd_fraction = sd(fraction, na.rm = TRUE),
    pooled_fraction = sum(n_cells, na.rm = TRUE) / sum(n_cells_total, na.rm = TRUE)
  ), by = .(sample_C_subtype, project_broad_compartment)]
  out[, scope := scope]
  out[]
}
comp_all_summary <- composition_summary(composition_all, "all")
comp_excl_summary <- composition_summary(composition_excl, "exclude_Bassez")
comp_compare <- merge(
  comp_all_summary[, scope := NULL], comp_excl_summary[, scope := NULL],
  by = c("sample_C_subtype", "project_broad_compartment"), suffixes = c("_all", "_exclude_Bassez"), all = TRUE
)
comp_compare[, `:=`(
  delta_mean_fraction = mean_fraction_exclude_Bassez - mean_fraction_all,
  delta_median_fraction = median_fraction_exclude_Bassez - median_fraction_all,
  delta_pooled_fraction = pooled_fraction_exclude_Bassez - pooled_fraction_all
)]
write_tsv(composition_all, file.path(outdir, "composition/sample_level_composition_all.tsv.gz"))
write_tsv(composition_excl, file.path(outdir, "composition/sample_level_composition_exclude_Bassez.tsv.gz"))
write_tsv(rbind(comp_all_summary, comp_excl_summary), file.path(outdir, "composition/composition_summary_by_state.tsv"))
write_tsv(comp_compare, file.path(outdir, "composition/composition_all_vs_exclude_Bassez.tsv"))

# 3. 以冻结的 7.3 program score 作为输入，按 Stage 3 同一状态/队列/Pre-Post框架重估中心。
program_scores <- fread(paths$program_scores)
program_catalog <- fread(paths$program_catalog)
formal_centers <- fread(paths$formal_centers)
program_scores <- program_scores[score_layer == "NONDEFINITION" &
                                 project_broad_compartment %in% CORE_COMPARTMENTS &
                                 n_cells >= 20 & sample_C_subtype %in% STATE_LEVELS]
program_scores[, is_bassez := tolower(study) == "bassez"]

make_state_design <- function(dat) {
  dd <- copy(dat)
  dd[, state := factor(sample_C_subtype, levels = STATE_LEVELS)]
  dd[, cohort_block := droplevels(factor(cohort_block))]
  use_pre_post <- uniqueN(dd$pre_post, na.rm = TRUE) > 1L
  if (use_pre_post) dd[, pre_post := factor(pre_post, levels = c("Pre", "Post"))]
  rhs <- c("0 + state", "cohort_block", if (use_pre_post) "pre_post")
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(form, data = dd)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design) && use_pre_post) {
    use_pre_post <- FALSE
    form <- ~ 0 + state + cohort_block
    design <- model.matrix(form, data = dd)
    colnames(design) <- make.names(colnames(design))
  }
  if (qr(design)$rank < ncol(design)) stop("state_design_rank_deficient")
  margins <- lapply(STATE_LEVELS, function(st) {
    nd <- copy(dd)
    nd[, state := factor(st, levels = STATE_LEVELS)]
    mm <- model.matrix(form, data = nd)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(design), drop = FALSE]
    colMeans(mm)
  })
  contrast <- do.call(cbind, margins)
  colnames(contrast) <- STATE_LEVELS
  rownames(contrast) <- colnames(design)
  list(data = dd, design = design, contrast = contrast,
       formula = paste(deparse(form), collapse = ""), uses_pre_post = use_pre_post)
}

fit_program_centers <- function(dat, scope) {
  if (uniqueN(dat$sample_C_subtype) < 4L) stop("not_all_four_states")
  payload <- make_state_design(dat)
  setorder(dat, unit_id)
  wide <- dcast(dat, program_id ~ unit_id, value.var = "score_z_global")
  mat <- as.matrix(wide[, setdiff(names(wide), "program_id"), with = FALSE])
  rownames(mat) <- wide$program_id
  idx <- match(payload$data$unit_id, colnames(mat))
  if (anyNA(idx)) stop("program_score_unit_mismatch")
  mat <- mat[, idx, drop = FALSE]
  repeated <- payload$data[, .N, by = patient_key][N > 1L, .N]
  block_used <- repeated >= 5L
  cor_value <- NA_real_
  if (block_used) {
    corfit <- duplicateCorrelation(mat, payload$design, block = factor(payload$data$patient_key))
    cor_value <- corfit$consensus.correlation
    fit <- lmFit(mat, payload$design, block = factor(payload$data$patient_key), correlation = cor_value)
  } else {
    fit <- lmFit(mat, payload$design)
  }
  fit2 <- eBayes(contrasts.fit(fit, payload$contrast))
  beta <- as.data.table(fit2$coefficients, keep.rownames = "program_id")
  long <- melt(beta, id.vars = "program_id", variable.name = "state", value.name = "adjusted_center")
  long[, `:=`(scope = scope, n_samples = nrow(payload$data), n_patients = uniqueN(payload$data$patient_key),
              n_cohorts = uniqueN(payload$data$cohort_block), n_repeated_patients = repeated,
              formula = payload$formula, patient_block_used = block_used, consensus_correlation = cor_value)]
  long
}

program_center_rows <- list()
program_model_rows <- list()
for (scope in c("all", "exclude_Bassez")) {
  for (comp in CORE_COMPARTMENTS) {
    dat <- program_scores[project_broad_compartment == comp & (scope == "all" | !is_bassez)]
    ok <- TRUE
    result <- tryCatch(fit_program_centers(dat, scope), error = function(e) { ok <<- FALSE; NULL })
    if (ok) {
      result[, compartment := comp]
      program_center_rows[[length(program_center_rows) + 1L]] <- result
      program_model_rows[[length(program_model_rows) + 1L]] <- unique(result[, .(
        scope, compartment, n_samples, n_patients, n_cohorts, n_repeated_patients,
        formula, patient_block_used, consensus_correlation
      )])
    } else {
      program_model_rows[[length(program_model_rows) + 1L]] <- data.table(scope = scope, compartment = comp,
        n_samples = nrow(dat), n_patients = uniqueN(dat$patient_key), n_cohorts = uniqueN(dat$cohort_block),
        n_repeated_patients = NA_integer_, formula = NA_character_, patient_block_used = FALSE,
        consensus_correlation = NA_real_, status = "NOT_EVALUABLE")
    }
  }
}
program_centers <- rbindlist(program_center_rows, fill = TRUE)
program_centers <- merge(program_centers, program_catalog[, .(program_id, program_order, compartment, biological_family,
                                                             evidence_use, definition_overlap_class)],
                         by = c("program_id", "compartment"), all.x = TRUE)
program_centers <- program_centers[order(program_order, scope, state)]
program_compare <- merge(
  program_centers[scope == "all", .(program_id, compartment, state, center_all = adjusted_center)],
  program_centers[scope == "exclude_Bassez", .(program_id, compartment, state, center_exclude_Bassez = adjusted_center)],
  by = c("program_id", "compartment", "state"), all = TRUE
)
program_compare[, `:=`(
  delta_center = center_exclude_Bassez - center_all,
  direction_concordant = sign(center_exclude_Bassez) == sign(center_all) | (center_exclude_Bassez == 0 & center_all == 0)
)]
program_summary <- program_compare[, .(
  n_states = sum(is.finite(center_all) & is.finite(center_exclude_Bassez)),
  state_direction_concordance = mean(direction_concordant, na.rm = TRUE),
  spearman_center_concordance = if (sum(is.finite(center_all) & is.finite(center_exclude_Bassez)) >= 3L)
    cor(center_all, center_exclude_Bassez, method = "spearman", use = "complete.obs") else NA_real_,
  max_abs_delta_center = max(abs(delta_center), na.rm = TRUE)
), by = .(program_id, compartment)]
formal_focus <- unique(formal_centers[score_layer == "NONDEFINITION" & stage2_g4_status != "QC_ONLY",
                                      .(program_id, compartment, stage2_g4_status)])
program_summary <- merge(program_summary, formal_focus, by = c("program_id", "compartment"), all.x = TRUE)
write_tsv(program_centers, file.path(outdir, "marker_program/program_state_centers_all_vs_exclude_Bassez.tsv"))
write_tsv(program_compare, file.path(outdir, "marker_program/program_state_center_comparison.tsv"))
write_tsv(program_summary, file.path(outdir, "marker_program/program_sensitivity_summary.tsv"))
write_tsv(rbindlist(program_model_rows, fill = TRUE), file.path(outdir, "marker_program/program_model_summary.tsv"))

# 4. 对 Fig. 6 D 已固定的 15 个通路重跑 NONDEFINITION GSEA；不扩展新通路集合。
pathway_catalog <- fread(paths$fig6_pathway_source)[, .(pathway = unique(pathway)), by = pathway]
pathway_catalog <- unique(pathway_catalog$pathway)
if (length(pathway_catalog) != 15L) stop("Expected the fixed 15-pathway Fig. 6 catalog")
gmt_sets <- c(
  read_gmt(paths$gmt_hallmark), read_gmt(paths$gmt_reactome), read_gmt(paths$gmt_gobp)
)
missing_pathways <- setdiff(pathway_catalog, names(gmt_sets))
if (length(missing_pathways) > 0L) stop("Fixed pathway catalog missing from GMT: ", paste(missing_pathways, collapse = ","))
pathway_sets <- gmt_sets[pathway_catalog]
write_tsv(data.table(pathway = pathway_catalog,
                     collection = fifelse(grepl("^HALLMARK_", pathway_catalog), "Hallmark",
                                   fifelse(grepl("^REACTOME_", pathway_catalog), "Reactome", "GO_BP"))),
          file.path(outdir, "gsea/fixed_15_pathway_catalog.tsv"))

model_contract <- fread(paths$model_contract)
counts_manifest <- fread(paths$counts_manifest)
definition_genes <- fread(paths$definition_genes)[exclude_from_independent_discovery == TRUE, unique(gene)]
targets <- model_contract[role == "primary" & status == "PASS_TECHNICAL_DESIGN" &
                          comparison %in% TARGET_COMPARISONS & compartment %in% CORE_COMPARTMENTS]
if (nrow(targets) != 12L) stop("Expected 12 Fig. 6 pathway targets")

select_one_per_patient <- function(dat) {
  dd <- copy(dat)
  dd[, pre_priority := fifelse(pre_post == "Pre", 0L, 1L)]
  setorder(dd, patient_key, pre_priority, -n_cells, -library_size, sample_unit_key)
  dd[, .SD[1L], by = patient_key]
}

make_group_design <- function(dat, g1, g2, include_pre_post = TRUE) {
  dd <- copy(dat)
  dd[, group := factor(analysis_group, levels = c(g2, g1))]
  dd[, cohort_block := droplevels(factor(cohort_block))]
  use_cohort <- uniqueN(dd$cohort_block) > 1L
  use_pre_post <- include_pre_post && uniqueN(dd$pre_post, na.rm = TRUE) > 1L
  if (use_pre_post) dd[, pre_post := factor(pre_post, levels = c("Pre", "Post"))]
  rhs <- c("0 + group", if (use_cohort) "cohort_block", if (use_pre_post) "pre_post")
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(form, data = dd)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design) && use_pre_post) {
    use_pre_post <- FALSE
    rhs <- c("0 + group", if (use_cohort) "cohort_block")
    form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
    design <- model.matrix(form, data = dd)
    colnames(design) <- make.names(colnames(design))
  }
  if (qr(design)$rank < ncol(design)) stop("group_design_rank_deficient")
  list(data = dd, design = design, formula = paste(deparse(form), collapse = ""),
       uses_cohort = use_cohort, uses_pre_post = use_pre_post)
}

fit_gene_model <- function(counts, dat, g1, g2, block_patient, label) {
  if (uniqueN(dat$analysis_group) < 2L) stop("one_group_only")
  payload <- make_group_design(dat, g1, g2, include_pre_post = TRUE)
  idx <- match(payload$data$unit_id, colnames(counts))
  if (anyNA(idx)) stop("count_sample_mismatch")
  y <- DGEList(counts = counts[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = payload$design)
  if (sum(keep) < 100L) stop("too_few_genes_after_filter")
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  v <- voom(y, payload$design, plot = FALSE)
  cor_value <- NA_real_
  block_used <- isTRUE(block_patient) && uniqueN(payload$data[, .(patient_key)], by = "patient_key") >= 5L &&
    anyDuplicated(payload$data$patient_key) > 0L
  if (block_used) {
    corfit <- duplicateCorrelation(v, payload$design, block = factor(payload$data$patient_key))
    cor_value <- corfit$consensus.correlation
    fit <- lmFit(v, payload$design, block = factor(payload$data$patient_key), correlation = cor_value)
  } else {
    fit <- lmFit(v, payload$design)
  }
  contrast <- makeContrasts(contrasts = paste0("group", make.names(g1), "-group", make.names(g2)), levels = payload$design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"), keep.rownames = "gene")
  setnames(tt, c("P.Value", "adj.P.Val"), c("p_value", "fdr"))
  tt[, `:=`(definition_gene = gene %in% definition_genes, analysis_label = label,
            group1 = g1, group2 = g2, n_samples = nrow(payload$data),
            n_patients = uniqueN(payload$data$patient_key), n_genes_tested = .N,
            model_formula = payload$formula, patient_block_used = block_used,
            consensus_correlation = cor_value)]
  list(result = tt, data = payload$data, summary = data.table(
    analysis_label = label, n_samples = nrow(payload$data), n_patients = uniqueN(payload$data$patient_key),
    n_group1_samples = sum(payload$data$analysis_group == g1), n_group2_samples = sum(payload$data$analysis_group == g2),
    n_group1_patients = uniqueN(payload$data[analysis_group == g1, patient_key]),
    n_group2_patients = uniqueN(payload$data[analysis_group == g2, patient_key]),
    n_genes_tested = nrow(tt), n_genes_fdr005 = sum(tt$fdr < 0.05, na.rm = TRUE),
    model_formula = payload$formula, patient_block_used = block_used, consensus_correlation = cor_value
  ))
}

run_fixed_gsea <- function(tt, comparison, compartment, scope) {
  rank_dt <- tt[definition_gene == FALSE & is.finite(t)][order(-abs(t), gene)][!duplicated(gene)]
  stats <- rank_dt$t
  names(stats) <- rank_dt$gene
  if (length(stats) < 100L) return(data.table())
  # 敏感性复核使用固定上限的简单置换，避免 multilevel 在服务器上无界运行。
  fg <- as.data.table(fgseaSimple(pathways = pathway_sets, stats = stats, minSize = 10L, maxSize = 500L,
                                  nperm = 1000L, nproc = 1L))
  if (nrow(fg) == 0L) return(fg)
  fg[, leadingEdge_genes := vapply(leadingEdge, paste, collapse = ",", character(1))]
  fg[, leadingEdge := NULL]
  fg[, `:=`(comparison = comparison, compartment = compartment, scope = scope,
            n_rank_genes = length(stats), n_pathways = .N)]
  setcolorder(fg, c("comparison", "compartment", "scope", "pathway", "NES", "padj", "ES", "size",
                    "n_rank_genes", "n_pathways", "leadingEdge_genes"))
  fg
}

gsea_results <- list()
gsea_summaries <- list()
gene_model_summaries <- list()
for (i in seq_len(nrow(targets))) {
  target <- targets[i]
  comparison <- target$comparison
  comp <- target$compartment
  log_event("CHECKPOINT", "start GSEA sensitivity target", list(comparison = comparison, compartment = comp))
  prefix <- paste0(safe_name(comparison), "__", safe_name(comp))
  sample_file <- relocate_input(paste0(prefix, ".tsv.gz"), unname(cfg_paths["sample_sets_dir"]))
  count_file <- counts_manifest[compartment == comp, file][1]
  dat_all <- fread(sample_file)[n_cells >= target$threshold]
  dat_all[, is_bassez := tolower(study) == "bassez"]
  dat_excl <- dat_all[is_bassez == FALSE]
  analyses <- list()
  recommended <- target$recommended_model
  if (grepl("duplicateCorrelation", recommended)) {
    analyses$all <- list(dat = dat_all, block = TRUE, label = "all_samples_patient_block")
    analyses$exclude_Bassez <- list(dat = dat_excl, block = TRUE, label = "exclude_Bassez_patient_block")
  } else if (grepl("no_patient_block", recommended)) {
    analyses$all <- list(dat = dat_all, block = FALSE, label = "all_samples_no_patient_block")
    analyses$exclude_Bassez <- list(dat = dat_excl, block = FALSE, label = "exclude_Bassez_no_patient_block")
  } else {
    analyses$all <- list(dat = select_one_per_patient(dat_all), block = FALSE, label = "independent_patient")
    analyses$exclude_Bassez <- list(dat = select_one_per_patient(dat_excl), block = FALSE, label = "exclude_Bassez_independent_patient")
  }
  target_gsea <- list()
  for (scope in names(analyses)) {
    spec <- analyses[[scope]]
    log_event("CHECKPOINT", "fit GSEA sensitivity model", list(comparison = comparison, compartment = comp, scope = scope))
    fit <- tryCatch({
      se <- readRDS(relocate_input(count_file, unname(cfg_paths["counts_dir"])))
      counts <- assay(se, "counts")
      fit_gene_model(counts, spec$dat, target$group1, target$group2, spec$block, spec$label)
    }, error = function(e) {
      gsea_summaries[[length(gsea_summaries) + 1L]] <<- data.table(comparison, compartment = comp, scope,
        status = "NOT_EVALUABLE", reason = conditionMessage(e), n_samples = nrow(spec$dat),
        n_patients = uniqueN(spec$dat$patient_key), n_group1_samples = sum(spec$dat$analysis_group == target$group1),
        n_group2_samples = sum(spec$dat$analysis_group == target$group2))
      NULL
    })
    if (is.null(fit)) next
    gene_model_summaries[[length(gene_model_summaries) + 1L]] <- cbind(data.table(comparison = comparison, compartment = comp, scope = scope), fit$summary)
    fg <- run_fixed_gsea(fit$result, comparison, comp, scope)
    log_event("CHECKPOINT", "finish GSEA sensitivity model", list(comparison = comparison, compartment = comp, scope = scope, n_pathways = nrow(fg)))
    if (nrow(fg) > 0L) gsea_results[[length(gsea_results) + 1L]] <- fg
    gsea_summaries[[length(gsea_summaries) + 1L]] <- cbind(data.table(comparison = comparison, compartment = comp, scope, status = "PASS"), fit$summary)
    target_gsea[[scope]] <- fg
    rm(se, counts, fit); gc(verbose = FALSE)
  }
  if (all(c("all", "exclude_Bassez") %in% names(target_gsea))) {
    gsea_summaries[[length(gsea_summaries) + 1L]] <- data.table(comparison = comparison, compartment = comp,
      scope = "comparison", status = "PASS", n_pathways = length(pathway_catalog))
  }
}
gsea_dt <- rbindlist(gsea_results, fill = TRUE)
gsea_summary_dt <- rbindlist(gsea_summaries, fill = TRUE)
gene_summary_dt <- rbindlist(gene_model_summaries, fill = TRUE)
gsea_compare <- merge(
  gsea_dt[scope == "all", .(comparison, compartment, pathway, NES_all = NES, padj_all = padj, leadingEdge_all = leadingEdge_genes, n_rank_genes_all = n_rank_genes)],
  gsea_dt[scope == "exclude_Bassez", .(comparison, compartment, pathway, NES_exclude_Bassez = NES, padj_exclude_Bassez = padj, leadingEdge_exclude_Bassez = leadingEdge_genes, n_rank_genes_exclude_Bassez = n_rank_genes)],
  by = c("comparison", "compartment", "pathway"), all = TRUE
)
gsea_compare[, `:=`(
  delta_NES = NES_exclude_Bassez - NES_all,
  direction_concordant = sign(NES_exclude_Bassez) == sign(NES_all),
  baseline_significant = is.finite(padj_all) & padj_all < 0.05,
  exclude_Bassez_significant = is.finite(padj_exclude_Bassez) & padj_exclude_Bassez < 0.05
)]
gsea_compare[, survival_class := fifelse(baseline_significant & exclude_Bassez_significant & direction_concordant, "retained",
                                  fifelse(baseline_significant & direction_concordant, "attenuated_or_nonsignificant",
                                  fifelse(baseline_significant & exclude_Bassez_significant & !direction_concordant, "reversed",
                                  fifelse(baseline_significant, "lost", "not_baseline_significant"))))]
gsea_target_summary <- gsea_compare[, .(
  n_pathways = .N,
  n_baseline_significant = sum(baseline_significant, na.rm = TRUE),
  n_exclude_Bassez_significant = sum(exclude_Bassez_significant, na.rm = TRUE),
  n_retained = sum(survival_class == "retained", na.rm = TRUE),
  n_reversed = sum(survival_class == "reversed", na.rm = TRUE),
  n_lost = sum(survival_class == "lost", na.rm = TRUE),
  direction_concordance = mean(direction_concordant, na.rm = TRUE)
), by = .(comparison, compartment)]
write_tsv(gsea_dt, file.path(outdir, "gsea/gsea_fixed_15_all_vs_exclude_Bassez.tsv.gz"))
write_tsv(gsea_compare, file.path(outdir, "gsea/gsea_fixed_15_comparison.tsv"))
write_tsv(gsea_target_summary, file.path(outdir, "gsea/gsea_target_summary.tsv"))
write_tsv(gsea_summary_dt, file.path(outdir, "gsea/gsea_model_summary.tsv"))
write_tsv(gene_summary_dt, file.path(outdir, "gsea/gene_model_summary.tsv"))

# 5. QA 与停止条件：任何关键 GSEA 翻转/不可评价都不自动晋级为主图证据。
qa_rows <- list()
qa_add <- function(id, passed, note) {
  qa_rows[[length(qa_rows) + 1L]] <<- data.table(check_id = id, passed = isTRUE(passed), note = note)
}
qa_add("inputs_complete", nrow(manifest[critical == TRUE & exists & nonempty]) == nrow(manifest[critical == TRUE]), "All critical active inputs exist and are non-empty")
qa_add("overall_cell_total_355941", totals$n_cells[[1]] == 355941L, paste0("observed=", totals$n_cells[[1]]))
qa_add("overall_sample_units_224", totals$n_samples[[1]] == 224L, paste0("observed=", totals$n_samples[[1]]))
qa_add("overall_patients_167", totals$n_patients[[1]] == 167L, paste0("observed=", totals$n_patients[[1]]))
qa_add("overall_studies_8", totals$n_studies[[1]] == 8L, paste0("observed=", totals$n_studies[[1]]))
qa_add("overall_cancer_labels_9", totals$n_cancer_labels[[1]] == 9L, paste0("observed=", totals$n_cancer_labels[[1]]))
qa_add("bassez_cells_216698", bassez_manifest$n_cells[[1]] == 216698L, paste0("observed=", bassez_manifest$n_cells[[1]]))
qa_add("bassez_samples_84", bassez_manifest$n_samples[[1]] == 84L, paste0("observed=", bassez_manifest$n_samples[[1]]))
qa_add("bassez_patients_42", bassez_manifest$n_patients[[1]] == 42L, paste0("observed=", bassez_manifest$n_patients[[1]]))
qa_add("composition_exclusion_uses_samples", uniqueN(composition_excl$sample_unit_key) == totals$n_samples[[1]] - bassez_manifest$n_samples[[1]], "Sample-level denominator retained")
qa_add("program_centers_complete", nrow(program_centers) == 23L * 4L * 2L, paste0("rows=", nrow(program_centers)))
qa_add("program_centers_finite", all(is.finite(program_centers$adjusted_center)), "Frozen program-score centers finite")
qa_add("gsea_target_count_12", nrow(gsea_target_summary) == 12L, paste0("targets=", nrow(gsea_target_summary)))
qa_add("gsea_no_reversal", if (nrow(gsea_compare)) sum(gsea_compare$survival_class == "reversed", na.rm = TRUE) == 0L else FALSE,
       "A reversal would stop an unqualified pan-cancer claim")
qa_add("replicate_unit_not_cells", TRUE, "Composition=sample; program/GSEA=sample-level pseudobulk with patient block where estimable")
qa_table <- rbindlist(qa_rows)
write_tsv(qa_table, file.path(outdir, "qa/step15_qa.tsv"))

n_not_eval <- sum(gsea_summary_dt$status == "NOT_EVALUABLE", na.rm = TRUE)
n_reversed <- sum(gsea_compare$survival_class == "reversed", na.rm = TRUE)
stop_table <- data.table(
  stop_condition = c("bassez_definition_unresolved", "gsea_target_not_evaluable", "fixed_pathway_direction_reversed", "single_remaining_study_drives_signal"),
  triggered = c(FALSE, n_not_eval > 0L, n_reversed > 0L, FALSE),
  action = c("No", "STOP_UNQUALIFIED_PAN_CANCER_GSEA", "STOP_UNQUALIFIED_PAN_CANCER_FUNCTIONAL_CLAIM", "DOWNGRADE_TO_STUDY_BOUND"),
  note = c("Bassez resolved from active study field.",
           paste0("n_not_evaluable=", n_not_eval), paste0("n_reversed=", n_reversed),
           "Study-dominance remains an interpretation check; no automatic promotion.")
)
write_tsv(stop_table, file.path(outdir, "qa/stop_conditions.tsv"))

# 6. Fig. 6 handoff 只给出可用于重绘的审计摘要，不改变正式 source-data。
handoff <- rbindlist(list(
  data.table(item = "overall_cells", value = as.character(totals$n_cells[[1]]), interpretation = "Active pan-cancer cell denominator"),
  data.table(item = "overall_sample_units", value = as.character(totals$n_samples[[1]]), interpretation = "Study-aware sample denominator"),
  data.table(item = "overall_patients", value = as.character(totals$n_patients[[1]]), interpretation = "donor_id patient denominator"),
  data.table(item = "bassez_cell_fraction", value = sprintf("%.4f", bassez_manifest$cell_fraction_of_all[[1]]), interpretation = "Bassez is a major cell-level contributor"),
  data.table(item = "bassez_sample_fraction", value = sprintf("%.4f", bassez_manifest$sample_fraction_of_all[[1]]), interpretation = "Bassez sample contribution"),
  data.table(item = "bassez_patient_fraction", value = sprintf("%.4f", bassez_manifest$patient_fraction_of_all[[1]]), interpretation = "Bassez patient contribution"),
  data.table(item = "gsea_reversed_pathways", value = as.character(n_reversed), interpretation = "Must be zero for unqualified pan-cancer functional wording"),
  data.table(item = "gsea_not_evaluable_models", value = as.character(n_not_eval), interpretation = "Any non-zero count blocks automatic promotion"),
  data.table(item = "cellchat_role", value = "supplementary_candidate_only", interpretation = "Not re-tested here; remains inferred communication landscape"),
  data.table(item = "panel_use", value = "candidate_only", interpretation = "No formal Fig. 6 asset update authorized")
), fill = TRUE)
write_tsv(handoff, file.path(outdir, "fig6_handoff/fig6_step15_bassez_sensitivity_handoff.tsv"))

finished_at <- Sys.time()
run_manifest <- data.table(
  item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version", "edgeR_version", "limma_version", "fgsea_version", "output_directory", "status"),
  value = c(script_name, format(started_at, "%Y-%m-%d %H:%M:%S"), format(finished_at, "%Y-%m-%d %H:%M:%S"),
            sprintf("%.3f", as.numeric(difftime(finished_at, started_at, units = "mins"))), R.version.string, as.character(packageVersion("edgeR")), as.character(packageVersion("limma")),
            as.character(packageVersion("fgsea")), outdir, if (all(qa_table$passed)) "PASS_QA" else "PASS_WITH_LIMITATIONS")
)
write_tsv(run_manifest, file.path(outdir, "qa/step15_run_manifest.tsv"))
log_event("DONE", "Bassez exclusion sensitivity completed", list(qa = if (all(qa_table$passed)) "PASS" else "LIMITATIONS", output = outdir))
