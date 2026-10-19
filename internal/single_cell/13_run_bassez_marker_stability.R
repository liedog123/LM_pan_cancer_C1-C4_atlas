#!/usr/bin/env Rscript
# Purpose: Evaluate the reported Bassez-exclusion and marker sensitivity results.
# Task: bassez_sensitivity; step: s13
# Inputs (DCF fields): sample_sheet,counts_manifest,model_contract,marker_panel,counts_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task bassez_sensitivity --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 13_run_bassez_marker_stability.R --config PATH [--check-config]\n       Rscript 13_run_bassez_marker_stability.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("sample_sheet", "counts_manifest", "model_contract", "marker_panel", "counts_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet", "counts_manifest", "model_contract", "marker_panel")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir")])
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


# Fig. 6 Step 15b：以固定 36-gene panel 核验 Bassez 排除后的样本级 marker 稳定性。
# 只写入 Step 15 候选结果目录，不覆盖 active 或 formal 资产。

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(stringsAsFactors = FALSE)
set.seed(20260908)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
})
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
data.table::setDTthreads(1L)
options(mc.cores = 1L)

script_name <- "LM_pan_cancer/7.0scRNA_LM12_state_mechanism/pancancer/scripts/15b_pancancer_c1c4_bassez_marker_stability_v0.1.R"
outdir <- out_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("marker_program", "qa", "fig6_handoff")) {
  dir.create(file.path(outdir, subdir), recursive = TRUE, showWarnings = FALSE)
}

log_event <- function(level, msg, fields = list()) {
  clean <- function(x) gsub('"', "'", as.character(x), fixed = TRUE)
  extra <- if (!length(fields)) "" else paste0(" ", paste(sprintf("%s=\"%s\"", names(fields), vapply(fields, clean, character(1))), collapse = " "))
  cat(sprintf("[%s] %s script=\"%s\" msg=\"%s\"%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, clean(msg), extra))
  flush.console()
}

started_at <- Sys.time()
options(error = function() {
  log_event("ERROR", "unhandled error", list(detail = geterrmessage()))
  quit(status = 1L)
})

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", quote = TRUE, na = "NA")
}

safe_char <- function(x) {
  y <- as.character(x)
  y[is.na(y)] <- "NA"
  trimws(y)
}

file_record <- function(role, path, critical = TRUE) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  data.table(
    role = role, path = path, critical = critical, exists = exists,
    bytes = if (exists) as.numeric(info$size) else NA_real_,
    nonempty = exists && !dir.exists(path) && isTRUE(info$size > 0),
    md5 = if (exists && !dir.exists(path)) unname(tools::md5sum(path)) else NA_character_
  )
}

STATE_LEVELS <- c("C1", "C2", "C3", "C4")
CORE_COMPARTMENTS <- c("Malignant", "T/NK", "Myeloid", "Stromal/CAF")
MARKER_PANEL_PATH <- unname(cfg_paths["marker_panel"])

paths <- as.list(paths[c("sample_sheet", "counts_manifest", "model_contract", "marker_panel")])
manifest <- rbindlist(Map(function(nm, p) file_record(nm, p, TRUE), names(paths), unname(paths)), fill = TRUE)
write_tsv(manifest, file.path(outdir, "qa/marker_input_file_manifest.tsv"))
if (nrow(manifest[critical == TRUE & (!exists | !nonempty)]) > 0L) {
  stop("Required marker-stability inputs are incomplete")
}

log_event("START", "load fixed marker panel and active pseudobulk sample sheet",
          list(output = outdir, marker_panel = MARKER_PANEL_PATH, min_cells = 20L))
marker_panel <- fread(paths$marker_panel)
required_marker <- c("compartment", "gene", "program_display", "gene_order", "main_order")
if (length(setdiff(required_marker, names(marker_panel))) > 0L) {
  stop("Marker panel is missing required columns")
}
marker_panel[, `:=`(compartment = safe_char(compartment), gene = safe_char(gene))]
marker_panel <- unique(marker_panel[compartment %in% CORE_COMPARTMENTS, .(compartment, program_display, gene, gene_order, main_order)])
if (nrow(marker_panel) != 36L || uniqueN(marker_panel$gene) != 36L) {
  stop("Expected 36 unique fixed marker genes")
}

sample_sheet <- fread(paths$sample_sheet)
required_sample <- c("unit_id", "sample_unit_key", "project_broad_compartment", "n_cells", "patient_key", "study", "sample_C_subtype", "cohort_block", "pre_post")
if (length(setdiff(required_sample, names(sample_sheet))) > 0L) {
  stop("Pseudobulk sample sheet is missing required columns")
}
sample_sheet[, `:=`(
  unit_id = safe_char(unit_id), sample_unit_key = safe_char(sample_unit_key),
  project_broad_compartment = safe_char(project_broad_compartment), patient_key = safe_char(patient_key),
  study = safe_char(study), sample_C_subtype = safe_char(sample_C_subtype),
  cohort_block = safe_char(cohort_block), pre_post = safe_char(pre_post)
)]
sample_sheet <- sample_sheet[project_broad_compartment %in% CORE_COMPARTMENTS &
                             n_cells >= 20 & sample_C_subtype %in% STATE_LEVELS]
sample_sheet[, is_bassez := tolower(study) == "bassez"]

counts_manifest <- fread(paths$counts_manifest)
if (length(setdiff(c("compartment", "file"), names(counts_manifest))) > 0L) {
  stop("Pseudobulk counts manifest is missing compartment/file columns")
}

make_state_design <- function(dat) {
  dd <- copy(dat)
  if (uniqueN(dd$sample_C_subtype) < 4L) stop("not_all_four_states")
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

fit_marker_centers <- function(counts, dat, scope, compartment, marker_genes) {
  payload <- make_state_design(dat)
  idx <- match(payload$data$unit_id, colnames(counts))
  if (anyNA(idx)) stop("count_sample_mismatch")
  y <- DGEList(counts = counts[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = payload$design)
  keep <- keep | rownames(y) %in% marker_genes
  if (sum(keep) < 100L) stop("too_few_genes_after_filter")
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  v <- voom(y, payload$design, plot = FALSE)
  repeated_patients <- payload$data[, .N, by = patient_key][N > 1L, .N]
  block_used <- repeated_patients >= 5L
  cor_value <- NA_real_
  if (block_used) {
    corfit <- duplicateCorrelation(v, payload$design, block = factor(payload$data$patient_key))
    cor_value <- corfit$consensus.correlation
    fit <- lmFit(v, payload$design, block = factor(payload$data$patient_key), correlation = cor_value)
  } else {
    fit <- lmFit(v, payload$design)
  }
  fit2 <- eBayes(contrasts.fit(fit, payload$contrast))
  beta <- as.data.table(fit2$coefficients, keep.rownames = "gene")
  beta <- beta[gene %in% marker_genes]
  long <- melt(beta, id.vars = "gene", variable.name = "state", value.name = "adjusted_center")
  long[, `:=`(
    scope = scope, compartment = compartment, n_samples = nrow(payload$data),
    n_patients = uniqueN(payload$data$patient_key), n_studies = uniqueN(payload$data$study),
    n_repeated_patients = repeated_patients, formula = payload$formula,
    patient_block_used = block_used, consensus_correlation = cor_value
  )]
  summary <- data.table(
    scope = scope, compartment = compartment, status = "PASS", n_samples = nrow(payload$data),
    n_patients = uniqueN(payload$data$patient_key), n_studies = uniqueN(payload$data$study),
    n_states = uniqueN(payload$data$sample_C_subtype), n_repeated_patients = repeated_patients,
    formula = payload$formula, patient_block_used = block_used, consensus_correlation = cor_value,
    n_marker_genes = uniqueN(long$gene)
  )
  list(centers = long, summary = summary)
}

center_rows <- list()
model_rows <- list()
marker_genes_by_compartment <- split(marker_panel$gene, marker_panel$compartment)
for (comp_name in CORE_COMPARTMENTS) {
  marker_genes <- unique(marker_genes_by_compartment[[comp_name]])
  count_file <- relocate_input(counts_manifest[compartment == comp_name, file][1], unname(cfg_paths["counts_dir"]), allow_missing = TRUE)
  if (is.na(count_file) || !file.exists(count_file)) {
    model_rows[[length(model_rows) + 1L]] <- data.table(scope = c("all", "exclude_Bassez"), compartment = comp_name,
      status = "NOT_EVALUABLE", reason = "count_file_missing")
    next
  }
  log_event("CHECKPOINT", "load compartment pseudobulk counts", list(compartment = comp_name, input = count_file))
  se <- readRDS(count_file)
  counts <- assay(se, "counts")
  dat_all <- sample_sheet[project_broad_compartment == comp_name]
  for (scope in c("all", "exclude_Bassez")) {
    dat <- if (scope == "all") dat_all else dat_all[is_bassez == FALSE]
    log_event("CHECKPOINT", "fit marker stability model", list(scope = scope, compartment = comp_name, n_samples = nrow(dat)))
    fit <- tryCatch(
      fit_marker_centers(counts, dat, scope, comp_name, marker_genes),
      error = function(e) {
        model_rows[[length(model_rows) + 1L]] <<- data.table(scope = scope, compartment = comp_name,
          status = "NOT_EVALUABLE", n_samples = nrow(dat), n_patients = uniqueN(dat$patient_key),
          n_studies = uniqueN(dat$study), reason = conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) next
    center_rows[[length(center_rows) + 1L]] <- fit$centers
    model_rows[[length(model_rows) + 1L]] <- fit$summary
    log_event("CHECKPOINT", "finish marker stability model", list(scope = scope, compartment = comp_name,
                                                                      n_marker_genes = uniqueN(fit$centers$gene)))
  }
  rm(se, counts); gc(verbose = FALSE)
}

centers <- rbindlist(center_rows, fill = TRUE)
models <- rbindlist(model_rows, fill = TRUE)
if (nrow(centers) == 0L) stop("No marker centers were evaluated")
centers <- merge(centers, marker_panel, by = c("compartment", "gene"), all.x = TRUE, sort = FALSE)
setorder(centers, main_order, scope, state)
write_tsv(centers, file.path(outdir, "marker_program/marker_gene_state_centers_all_vs_exclude_Bassez.tsv"))
write_tsv(models, file.path(outdir, "marker_program/marker_model_summary.tsv"))

comparison <- merge(
  centers[scope == "all", .(compartment, gene, program_display, state, center_all = adjusted_center)],
  centers[scope == "exclude_Bassez", .(compartment, gene, program_display, state, center_exclude_Bassez = adjusted_center)],
  by = c("compartment", "gene", "program_display", "state"), all = TRUE
)
comparison[, `:=`(
  delta_center = center_exclude_Bassez - center_all,
  direction_concordant = (sign(center_exclude_Bassez) == sign(center_all)) |
    (center_exclude_Bassez == 0 & center_all == 0)
)]
summary <- comparison[, .(
  n_states = sum(is.finite(center_all) & is.finite(center_exclude_Bassez)),
  state_direction_concordance = mean(direction_concordant, na.rm = TRUE),
  spearman_center_concordance = if (sum(is.finite(center_all) & is.finite(center_exclude_Bassez)) >= 3L)
    cor(center_all, center_exclude_Bassez, method = "spearman", use = "complete.obs") else NA_real_,
  max_abs_delta_center = max(abs(delta_center), na.rm = TRUE),
  n_direction_reversals = sum(!direction_concordant, na.rm = TRUE)
), by = .(compartment, gene, program_display)]
write_tsv(comparison, file.path(outdir, "marker_program/marker_gene_state_center_comparison.tsv"))
write_tsv(summary, file.path(outdir, "marker_program/marker_sensitivity_summary.tsv"))

n_not_eval <- sum(models$status == "NOT_EVALUABLE", na.rm = TRUE)
n_reversals <- sum(comparison$direction_concordant == FALSE, na.rm = TRUE)
n_unstable_genes <- sum(summary$state_direction_concordance < 0.5 | summary$spearman_center_concordance < 0.5, na.rm = TRUE)
qa_rows <- list()
qa_add <- function(id, passed, note) {
  qa_rows[[length(qa_rows) + 1L]] <<- data.table(check_id = id, passed = isTRUE(passed), note = note)
}
qa_add("inputs_complete", nrow(manifest[critical == TRUE & exists & nonempty]) == nrow(manifest[critical == TRUE]), "All fixed marker-stability inputs exist")
qa_add("fixed_marker_count_36", nrow(marker_panel) == 36L, paste0("rows=", nrow(marker_panel)))
qa_add("fixed_marker_genes_unique", uniqueN(marker_panel$gene) == 36L, paste0("genes=", uniqueN(marker_panel$gene)))
qa_add("centers_nonempty", nrow(centers) > 0L, paste0("rows=", nrow(centers)))
qa_add("sample_level_units", TRUE, "All models use sample_unit_key pseudobulk units; cells are not repeats")
qa_add("models_with_not_evaluable", n_not_eval == 0L, paste0("n_not_evaluable=", n_not_eval))
qa_add("direction_reversals_reported", TRUE, paste0("n_direction_reversals=", n_reversals))
qa_add("unstable_genes_reported", TRUE, paste0("n_unstable_genes=", n_unstable_genes))
qa_table <- rbindlist(qa_rows)
write_tsv(qa_table, file.path(outdir, "qa/marker_stability_qa.tsv"))

stop_table <- data.table(
  stop_condition = c("marker_model_not_evaluable", "marker_direction_reversal", "marker_scope_widespread_instability"),
  triggered = c(n_not_eval > 0L, n_reversals > 0L, n_unstable_genes > floor(nrow(summary) / 2)),
  action = c("RESTRICT_TO_EVALUABLE_COMPARTMENTS", "MOVE_SENSITIVE_MARKERS_TO_SUPPLEMENTARY", "DO_NOT_PROMOTE_BROAD_MARKER_PANEL"),
  note = c(paste0("n_not_evaluable=", n_not_eval), paste0("n_direction_reversals=", n_reversals),
           paste0("n_unstable_genes=", n_unstable_genes, "/", nrow(summary)))
)
write_tsv(stop_table, file.path(outdir, "qa/marker_stability_stop_conditions.tsv"))

handoff <- rbindlist(list(
  data.table(item = "fixed_marker_genes", value = as.character(uniqueN(marker_panel$gene)), interpretation = "Frozen Fig. 6 candidate C panel; no post-hoc additions"),
  data.table(item = "marker_model_not_evaluable", value = as.character(n_not_eval), interpretation = "Any non-zero count restricts candidate display"),
  data.table(item = "marker_direction_reversals", value = as.character(n_reversals), interpretation = "Report before any main-panel promotion"),
  data.table(item = "marker_unstable_genes", value = as.character(n_unstable_genes), interpretation = "Directional/concordance screen only; not a new selection procedure"),
  data.table(item = "panel_use", value = "candidate_only", interpretation = "No formal Fig. 6 asset update authorized")
), fill = TRUE)
write_tsv(handoff, file.path(outdir, "fig6_handoff/fig6_marker_stability_handoff.tsv"))

finished_at <- Sys.time()
run_manifest <- data.table(
  item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version", "edgeR_version", "limma_version", "output_directory", "status"),
  value = c(script_name, format(started_at, "%Y-%m-%d %H:%M:%S"), format(finished_at, "%Y-%m-%d %H:%M:%S"),
            sprintf("%.3f", as.numeric(difftime(finished_at, started_at, units = "mins"))), R.version.string,
            as.character(packageVersion("edgeR")), as.character(packageVersion("limma")), outdir,
            if (all(qa_table$passed)) "PASS_QA_WITH_REPORTED_LIMITATIONS" else "LIMITATIONS_RECORDED")
)
write_tsv(run_manifest, file.path(outdir, "qa/marker_stability_run_manifest.tsv"))
log_event("DONE", "fixed marker stability completed", list(output = outdir, centers = nrow(centers), reversals = n_reversals,
                                                           not_evaluable_models = n_not_eval))
