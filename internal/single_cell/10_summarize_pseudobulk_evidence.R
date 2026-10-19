#!/usr/bin/env Rscript
# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Task: single_cell_pseudobulk; step: s10
# Inputs (DCF fields): primary_gate,stability_gate,candidates,sensitivity,cohort_meta,loso,target_stability_summary,legacy_deg,primary_gene_results_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 10_summarize_pseudobulk_evidence.R --config PATH [--check-config]\n       Rscript 10_summarize_pseudobulk_evidence.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("primary_gate", "stability_gate", "candidates", "sensitivity", "cohort_meta", "loso", "target_stability_summary", "legacy_deg", "primary_gene_results_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("primary_gate", "stability_gate", "candidates", "sensitivity", "cohort_meta", "loso", "target_stability_summary", "legacy_deg")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("primary_gene_results_dir")])
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

# 普通执行仅在原读取位置解析选中产物；元信息检查不枚举输入目录。
relocate_input <- function(recorded, input_dir) {
  if(length(recorded)!=1L || is.na(recorded) || !nzchar(recorded)) stop("Invalid recorded input path.")
  leaf <- basename(recorded)
  if(leaf %in% c(".","..") || !nzchar(leaf)) stop("Invalid input basename.")
  candidate <- resolve_path(file.path(input_dir,leaf), dirname(config_path))
  if(!is_within(candidate,input_dir) || identical(candidate,input_dir)) stop("Input escapes configured directory.")
  if(!file_test("-f",candidate) || file.access(candidate,4L)!=0L || is.na(file.info(candidate)$size) || file.info(candidate)$size<=0) stop("Selected input must be a readable nonempty regular file.")
  candidate
}



suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02e_summarize_pseudobulk_evidence.R"

log_event <- function(level, msg, ...) {
  extras <- list(...)
  extra_text <- if (length(extras) == 0L) "" else paste(
    vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)),
    collapse = " "
  )
  cat(sprintf("[%s] %s script=%s msg=%s%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg),
              if (nzchar(extra_text)) " " else "", extra_text))
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  primary_gene_results_dir <- unname(paths["primary_gene_results_dir"])
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (fread(unname(paths["primary_gate"]))$status[[1]] != "PASS") stop("G2B gate is not PASS")
  if (fread(unname(paths["stability_gate"]))$status[[1]] != "PASS") stop("G2C gate is not PASS")
  fwrite(data.table(gate = "G2D_gene_evidence_synthesis", status = "RUNNING",
                    note = "RUNNING prevents stale evidence status after interruption."),
         file.path(outdir, "g2d_evidence_gate.tsv"), sep = "\t", na = "NA")
  log_event("START", "summarize pseudobulk stability and legacy survival evidence")

  candidates <- fread(unname(paths["candidates"]))
  sensitivity <- fread(unname(paths["sensitivity"]))
  cohort_meta <- fread(unname(paths["cohort_meta"]))
  loso <- fread(unname(paths["loso"]))

  sens_sum <- sensitivity[, .(
    n_sensitivity = sum(!is.na(sensitivity_logFC)),
    sensitivity_sign_fraction = mean(direction_concordant, na.rm = TRUE),
    sensitivity_fdr010_fraction = mean(retained_fdr_010, na.rm = TRUE),
    sensitivity_median_logFC = median(sensitivity_logFC, na.rm = TRUE)
  ), by = .(comparison, compartment, gene)]

  evidence <- merge(candidates, sens_sum, by = c("comparison", "compartment", "gene"), all.x = TRUE)
  evidence <- merge(evidence, cohort_meta, by = c("comparison", "compartment", "gene"), all.x = TRUE)
  evidence <- merge(evidence, loso, by = c("comparison", "compartment", "gene"), all.x = TRUE)
  evidence[, random_direction_concordant := sign(main_logFC) == sign(beta_random)]
  evidence[, evidence_tier := fcase(
    main_fdr < 0.10 & sensitivity_sign_fraction >= 0.75 & loso_sign_fraction >= 0.80 &
      n_cohorts >= 2 & cohort_sign_fraction >= 0.67 & random_direction_concordant & p_random < 0.10,
    "A_cross_cohort_stable",
    main_fdr < 0.10 & sensitivity_sign_fraction >= 0.75 & loso_sign_fraction >= 0.80,
    "B_model_stable_context_dependent",
    main_fdr < 0.10,
    "C_overall_only",
    default = "D_leading_nonsignificant"
  )]
  evidence[, evidence_tier := factor(evidence_tier, levels = c(
    "A_cross_cohort_stable", "B_model_stable_context_dependent",
    "C_overall_only", "D_leading_nonsignificant"
  ))]
  evidence[, abs_main_logFC_sort := abs(main_logFC)]
  setorder(evidence, comparison, compartment, evidence_tier, main_fdr, -abs_main_logFC_sort)
  evidence[, abs_main_logFC_sort := NULL]
  evidence[, evidence_tier := as.character(evidence_tier)]
  fwrite(evidence, file.path(outdir, "integrated_gene_evidence.tsv.gz"),
         sep = "\t", na = "NA", compress = "gzip")

  target_summary <- evidence[, .(
    n_candidates = .N,
    n_main_fdr010_in_candidates = sum(main_fdr < 0.10, na.rm = TRUE),
    n_tier_a = sum(evidence_tier == "A_cross_cohort_stable"),
    n_tier_b = sum(evidence_tier == "B_model_stable_context_dependent"),
    n_tier_c = sum(evidence_tier == "C_overall_only"),
    n_tier_d = sum(evidence_tier == "D_leading_nonsignificant"),
    median_sensitivity_sign_fraction_main_sig = median(sensitivity_sign_fraction[main_fdr < 0.10], na.rm = TRUE),
    median_loso_sign_fraction_main_sig = median(loso_sign_fraction[main_fdr < 0.10], na.rm = TRUE),
    median_cohort_sign_fraction_main_sig = median(cohort_sign_fraction[main_fdr < 0.10], na.rm = TRUE)
  ), by = .(comparison, compartment)]
  total_signal <- fread(unname(paths["target_stability_summary"]))[,
    .(comparison, compartment, n_main_fdr010_total_nondefinition = n_main_fdr_010_nondefinition)]
  target_summary <- merge(target_summary, total_signal, by = c("comparison", "compartment"), all.x = TRUE)
  target_summary[, target_evidence_status := fcase(
    n_tier_a >= 10L, "CROSS_COHORT_SUPPORTED",
    n_tier_a > 0L | n_tier_b >= 10L, "MODEL_STABLE_CONTEXT_DEPENDENT",
    n_main_fdr010_total_nondefinition > 0L, "OVERALL_ONLY_NOT_STABLE",
    default = "NO_MAIN_FDR_SIGNAL"
  )]
  fwrite(target_summary, file.path(outdir, "target_gene_evidence_summary.tsv"), sep = "\t", na = "NA")

  # 旧 pooled-cell 证据只有 C1_vs_C4 与当前主端点完全对应；原子比较不强行映射。
  legacy_file <- unname(paths["legacy_deg"])
  legacy <- fread(legacy_file)[comparison == "C1_vs_C4" &
                                project_broad_compartment %in% c("Malignant", "T/NK", "Myeloid", "Stromal/CAF") &
                                p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25]
  legacy_rows <- list()
  for (comp in unique(legacy$project_broad_compartment)) {
    prefix <- paste0("C1_vs_C4__", safe_name(comp))
    new <- fread(relocate_input(paste0(prefix, "__main.tsv.gz"), primary_gene_results_dir))[,
      .(gene, new_logFC = logFC, new_p_value = p_value, new_fdr = fdr, definition_gene)]
    old <- legacy[project_broad_compartment == comp, .(
      comparison, compartment = project_broad_compartment, gene,
      legacy_logFC = avg_log2FC, legacy_fdr = p_val_adj
    )]
    audit <- merge(old, new, by = "gene", all.x = TRUE)
    audit[, direction_concordant := sign(legacy_logFC) == sign(new_logFC)]
    audit[, survival_class := fcase(
      is.na(new_logFC), "lost_not_tested",
      direction_concordant & new_fdr < 0.10, "retained",
      direction_concordant, "attenuated",
      !direction_concordant & new_fdr < 0.10, "reversed",
      default = "lost"
    )]
    legacy_rows[[length(legacy_rows) + 1L]] <- audit
  }
  legacy_audit <- rbindlist(legacy_rows, fill = TRUE)
  fwrite(legacy_audit, file.path(outdir, "legacy_c1c4_pooled_deg_survival.tsv.gz"),
         sep = "\t", na = "NA", compress = "gzip")
  legacy_summary <- legacy_audit[, .(
    n_legacy_significant = .N,
    n_retained = sum(survival_class == "retained"),
    n_attenuated = sum(survival_class == "attenuated"),
    n_reversed = sum(survival_class == "reversed"),
    n_lost = sum(survival_class %in% c("lost", "lost_not_tested")),
    retained_fraction = mean(survival_class == "retained"),
    direction_concordant_fraction = mean(direction_concordant, na.rm = TRUE),
    n_definition_genes = sum(definition_gene == TRUE, na.rm = TRUE)
  ), by = .(comparison, compartment)]
  fwrite(legacy_summary, file.path(outdir, "legacy_c1c4_survival_summary.tsv"), sep = "\t", na = "NA")

  out_gate <- data.table(
    gate = "G2D_gene_evidence_synthesis", status = "PASS",
    n_targets = nrow(target_summary), n_targets_cross_cohort_supported = sum(target_summary$target_evidence_status == "CROSS_COHORT_SUPPORTED"),
    n_targets_model_stable_context_dependent = sum(target_summary$target_evidence_status == "MODEL_STABLE_CONTEXT_DEPENDENT"),
    note = "Gene-level synthesis completed; pathway interpretation and final G2 decision remain pending."
  )
  fwrite(out_gate, file.path(outdir, "g2d_evidence_gate.tsv"), sep = "\t", na = "NA")
  fwrite(data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string, outdir)
  ), file.path(outdir, "02e_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pseudobulk gene evidence synthesis completed", status = "PASS", output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
