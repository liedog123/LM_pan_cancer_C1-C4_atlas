#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

usage <- function() {
  cat("Usage: Rscript 01_c1c4_from_corrected_predictions.R --config PATH [--check-config]\n",
      "       Rscript 01_c1c4_from_corrected_predictions.R --help\n\n",
      "Continue Stage 5 C1/C4 statistics from 20 corrected prediction matrices.\n",
      "DCF fields: canonical_labels, model_input, prediction_dir, drug_annotation, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks path metadata only; it reads no research data and creates no outputs.\n",
      "An ordinary configured invocation runs the full analysis. Outputs must be independent and empty.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (!(length(args) %in% c(2L, 3L)) || args[1L] != "--config" ||
    (length(args) == 3L && args[3L] != "--check-config")) {
  usage()
  stop("Expected --config PATH, optionally followed by --check-config.", call. = FALSE)
}

# Resolve relative paths from the configuration directory and resolve existing ancestors for symlinks.
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
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
config_path <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
cfg <- read.dcf(config_path)
required_fields <- c("canonical_labels", "model_input", "prediction_dir",
                     "drug_annotation", "output_dir")
if (nrow(cfg) != 1L || !setequal(colnames(cfg), required_fields)) {
  stop("Configuration must contain one DCF record with exactly: ",
       paste(required_fields, collapse = ", "), call. = FALSE)
}
cfg <- setNames(as.character(cfg[1L, required_fields]), required_fields)
if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME", cfg, fixed = TRUE))) {
  stop("Replace all example placeholders with explicit paths.", call. = FALSE)
}
paths <- vapply(cfg, resolve_path, character(1L), base_dir = dirname(config_path))
canonical_path <- unname(paths["canonical_labels"])
model_input_path <- unname(paths["model_input"])
prediction_dir <- unname(paths["prediction_dir"])
annotation_path <- unname(paths["drug_annotation"])
stage5_root <- unname(paths["output_dir"])

# Reject input/output overlap, broad root paths and nonempty output directories; no overwrite option.
protected_paths <- c(dirname(canonical_path), dirname(model_input_path),
                     prediction_dir, dirname(annotation_path), config_path)
if (stage5_root == "/" || stage5_root == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, stage5_root) || is_within(stage5_root, p)
    }, logical(1L)))) {
  stop("output_dir must be an independent directory outside all inputs and the configuration file.",
       call. = FALSE)
}
if (file.exists(stage5_root) && !dir.exists(stage5_root)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(stage5_root) &&
    length(list.files(stage5_root, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("output_dir is nonempty; choose a new independent directory.", call. = FALSE)
}
input_files <- c(canonical_path, model_input_path, annotation_path)
if (any(!file.exists(input_files)) || any(dir.exists(input_files)) ||
    any(file.access(input_files, mode = 4L) != 0L) ||
    any(is.na(file.info(input_files)$size) | file.info(input_files)$size <= 0)) {
  stop("All three input files must be readable, nonempty regular files.", call. = FALSE)
}
if (!dir.exists(prediction_dir) || file.access(prediction_dir, mode = 4L) != 0L) {
  stop("prediction_dir must be a readable directory.", call. = FALSE)
}
prediction_paths <- list.files(prediction_dir,
  pattern = "^TCGA-[A-Z0-9]+_corrected_DrugPredictions\\.csv$", full.names = TRUE)
if (length(prediction_paths) != 20L || any(dir.exists(prediction_paths)) ||
    any(file.access(prediction_paths, mode = 4L) != 0L) ||
    any(is.na(file.info(prediction_paths)$size) | file.info(prediction_paths)$size <= 0)) {
  stop("Expected exactly 20 readable, nonempty corrected prediction CSV files.", call. = FALSE)
}
# Individual prediction files may be symlinks to locations outside their directory.
prediction_real_paths <- vapply(prediction_paths, normalizePath, character(1L),
                               winslash = "/", mustWork = TRUE)
if (any(vapply(prediction_real_paths, function(p) {
  is_within(p, stage5_root) || is_within(stage5_root, p)
}, logical(1L)))) {
  stop("output_dir overlaps a resolved prediction file.", call. = FALSE)
}
existing_parent <- stage5_root
while (!dir.exists(existing_parent)) existing_parent <- dirname(existing_parent)
if (file.access(existing_parent, mode = 2L) != 0L) {
  stop("output_dir has no writable existing parent.", call. = FALSE)
}
if (length(args) == 3L) {
  cat("Configuration paths passed metadata checks. No analysis was run; no outputs were created.\n")
  quit(status = 0L)
}

suppressPackageStartupMessages({
  library(data.table)
})

version_tag <- "v0.1_20260731"
readiness_dir <- file.path(stage5_root, "results", "00_readiness", version_tag)
analysis_dir <- file.path(stage5_root, "results", "01_c1c4_drug_response", version_tag)
dir.create(readiness_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

min_group_n <- 30L
fdr_cut <- 0.05
delta_cut <- 0.20

log_msg <- function(level, step, ...) {
  cat(sprintf(
    "%s\t%s\t%s\t%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level, step, paste(..., collapse = " ")
  ))
  flush.console()
}

assert_nonempty <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  size <- file.info(path)$size
  if (is.na(size) || size <= 0) stop("Empty input: ", path)
  invisible(size)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

normalize_sample15 <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 1L, 15L)
}

sample_type_priority <- function(sample15) {
  priority <- c("01", "03", "06", "02", "05", "07", "08", "09")
  ans <- match(substr(sample15, 14L, 15L), priority)
  ans[is.na(ans)] <- length(priority) + 1L
  ans
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x)
  n2 <- length(y)
  if (n1 == 0L || n2 == 0L) return(NA_real_)
  ranks <- rank(c(x, y), ties.method = "average")
  u1 <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n2) - 1
}

test_drug <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x)
  n4 <- length(y)
  base <- list(
    status = "ok",
    reason = "",
    n_C1 = n1,
    n_C4 = n4,
    median_C1 = if (n1) median(x) else NA_real_,
    median_C4 = if (n4) median(y) else NA_real_,
    median_diff_C1_minus_C4 = NA_real_,
    hl_shift_C1_minus_C4 = NA_real_,
    hl_ci_low = NA_real_,
    hl_ci_high = NA_real_,
    cliffs_delta_C1_minus_C4 = NA_real_,
    p_value = NA_real_
  )
  if (n1 < min_group_n || n4 < min_group_n) {
    base$status <- "not_tested"
    base$reason <- sprintf("group_size_below_threshold:C1=%d;C4=%d", n1, n4)
    return(base)
  }
  if (length(unique(c(x, y))) < 2L) {
    base$status <- "not_tested"
    base$reason <- "constant_prediction"
    return(base)
  }
  wt <- tryCatch(
    suppressWarnings(wilcox.test(
      x, y,
      alternative = "two.sided",
      exact = FALSE,
      conf.int = TRUE,
      conf.level = 0.95
    )),
    error = function(e) e
  )
  if (inherits(wt, "error")) {
    base$status <- "failed"
    base$reason <- conditionMessage(wt)
    return(base)
  }
  base$median_diff_C1_minus_C4 <- base$median_C1 - base$median_C4
  base$hl_shift_C1_minus_C4 <- if (length(wt$estimate)) unname(wt$estimate)[1L] else NA_real_
  base$hl_ci_low <- if (!is.null(wt$conf.int)) unname(wt$conf.int)[1L] else NA_real_
  base$hl_ci_high <- if (!is.null(wt$conf.int)) unname(wt$conf.int)[2L] else NA_real_
  base$cliffs_delta_C1_minus_C4 <- cliffs_delta(x, y)
  base$p_value <- unname(wt$p.value)
  base
}

tryCatch({
  log_msg("START", "stage5", "C1 vs C4 formal GDSC predicted-response analysis")

  prediction_files <- sort(list.files(
    prediction_dir,
    pattern = "^TCGA-[A-Z0-9]+_corrected_DrugPredictions\\.csv$",
    full.names = TRUE
  ))
  input_paths <- c(canonical_path, model_input_path, annotation_path, prediction_files)
  invisible(lapply(input_paths, assert_nonempty))
  if (length(prediction_files) != 20L) {
    stop("Expected 20 corrected prediction files, observed ", length(prediction_files))
  }

  input_inventory <- data.table(
    object_role = c(
      "canonical_labels",
      "legacy_stageA_model_input",
      "drug_annotation",
      rep("corrected_prediction_matrix", length(prediction_files))
    ),
    path = input_paths,
    size_bytes = file.info(input_paths)$size,
    md5 = unname(tools::md5sum(input_paths))
  )
  write_tsv(input_inventory, file.path(readiness_dir, "stage5_input_inventory.tsv"))

  canonical <- fread(canonical_path)
  required_label_cols <- c("Sample", "Cancer", "C_subtype")
  missing_label_cols <- setdiff(required_label_cols, names(canonical))
  if (length(missing_label_cols)) {
    stop("Canonical labels missing columns: ", paste(missing_label_cols, collapse = ","))
  }
  canonical[, Sample15 := normalize_sample15(Sample)]
  canonical[, Patient_ID := substr(Sample15, 1L, 12L)]
  canonical <- canonical[C_subtype %chin% c("C1", "C4")]
  if (any(!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}$", canonical$Sample15))) {
    stop("Invalid canonical Sample15 format")
  }

  patient_label_audit <- canonical[, .(
    n_label_rows = .N,
    n_cancers = uniqueN(Cancer),
    n_subtypes = uniqueN(C_subtype),
    cancers = paste(sort(unique(Cancer)), collapse = ";"),
    subtypes = paste(sort(unique(C_subtype)), collapse = ";")
  ), by = Patient_ID]
  patient_label_audit[, discordant := n_cancers > 1L | n_subtypes > 1L]
  discordant_patients <- patient_label_audit[discordant == TRUE, Patient_ID]
  write_tsv(patient_label_audit, file.path(readiness_dir, "canonical_c1c4_patient_label_audit.tsv"))

  annotation <- fread(annotation_path)
  if (!"drug" %chin% names(annotation)) stop("Drug annotation missing drug column")
  annotation[, drug := as.character(drug)]
  annotation <- annotation[!duplicated(drug)]

  file_audit_rows <- list()
  coverage_rows <- list()
  analysis_rows <- list()
  matched_patient_rows <- list()

  for (prediction_file in prediction_files) {
    project <- sub("_corrected_DrugPredictions\\.csv$", "", basename(prediction_file))
    cancer <- sub("^TCGA-", "", project)
    log_msg("INFO", "cancer", cancer)

    pred <- fread(prediction_file, check.names = FALSE)
    if (ncol(pred) < 2L) stop("Prediction matrix has fewer than two columns: ", prediction_file)
    sample_col <- names(pred)[1L]
    setnames(pred, sample_col, "source_sample")
    drug_cols <- names(pred)[-1L]
    if (anyDuplicated(drug_cols)) stop("Duplicated drug columns in ", prediction_file)

    pred[, Sample15 := normalize_sample15(source_sample)]
    pred[, Patient_ID := substr(Sample15, 1L, 12L)]
    n_source_rows <- nrow(pred)
    n_duplicate_sample_rows <- sum(duplicated(pred$Sample15))
    pred <- pred[!duplicated(Sample15)]

    matched <- merge(
      pred,
      canonical[
        !Patient_ID %chin% discordant_patients,
        .(Sample15, Patient_ID_label = Patient_ID, Cancer, C_subtype)
      ],
      by = "Sample15",
      all = FALSE,
      sort = FALSE
    )
    matched <- matched[Cancer == cancer]
    matched[, Patient_ID := substr(Sample15, 1L, 12L)]
    matched[, sample_priority := sample_type_priority(Sample15)]
    setorder(matched, Patient_ID, sample_priority, Sample15)
    matched <- matched[, .SD[1L], by = Patient_ID]

    if (anyDuplicated(matched$Patient_ID)) stop("Patient deduplication failed for ", cancer)
    if (any(!matched$C_subtype %chin% c("C1", "C4"))) stop("Non-C1/C4 labels entered ", cancer)

    coverage <- matched[, .N, by = C_subtype]
    n_c1 <- coverage[C_subtype == "C1", sum(N)]
    n_c4 <- coverage[C_subtype == "C4", sum(N)]
    formal_cancer <- n_c1 >= min_group_n && n_c4 >= min_group_n
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      Cancer = cancer,
      n_C1 = n_c1,
      n_C4 = n_c4,
      n_total = n_c1 + n_c4,
      min_group_n = min_group_n,
      formal_eligible = formal_cancer
    )
    matched_patient_rows[[length(matched_patient_rows) + 1L]] <- matched[, .(
      Cancer = cancer,
      Patient_ID,
      Sample15,
      C_subtype
    )]

    finite_fraction <- vapply(
      pred[, ..drug_cols],
      function(z) mean(is.finite(suppressWarnings(as.numeric(z)))),
      numeric(1L)
    )
    file_audit_rows[[length(file_audit_rows) + 1L]] <- data.table(
      Cancer = cancer,
      source_file = prediction_file,
      source_rows = n_source_rows,
      duplicate_sample15_rows = n_duplicate_sample_rows,
      n_drugs = length(drug_cols),
      matched_unique_patients = nrow(matched),
      n_C1 = n_c1,
      n_C4 = n_c4,
      formal_eligible = formal_cancer,
      min_drug_finite_fraction = min(finite_fraction),
      max_drug_finite_fraction = max(finite_fraction)
    )

    for (drug in drug_cols) {
      values <- suppressWarnings(as.numeric(matched[[drug]]))
      group <- matched$C_subtype
      stat <- test_drug(values[group == "C1"], values[group == "C4"])
      analysis_rows[[length(analysis_rows) + 1L]] <- data.table(
        Cancer = cancer,
        drug = drug,
        status = stat$status,
        reason = stat$reason,
        n_C1 = stat$n_C1,
        n_C4 = stat$n_C4,
        median_C1 = stat$median_C1,
        median_C4 = stat$median_C4,
        median_diff_C1_minus_C4 = stat$median_diff_C1_minus_C4,
        hl_shift_C1_minus_C4 = stat$hl_shift_C1_minus_C4,
        hl_ci_low = stat$hl_ci_low,
        hl_ci_high = stat$hl_ci_high,
        cliffs_delta_C1_minus_C4 = stat$cliffs_delta_C1_minus_C4,
        p_value = stat$p_value
      )
    }
  }

  file_audit <- rbindlist(file_audit_rows, fill = TRUE)
  coverage <- rbindlist(coverage_rows, fill = TRUE)
  matched_patients <- rbindlist(matched_patient_rows, fill = TRUE)
  results <- rbindlist(analysis_rows, fill = TRUE)
  setorder(file_audit, Cancer)
  setorder(coverage, Cancer)
  setorder(matched_patients, Cancer, Patient_ID)
  setorder(results, Cancer, drug)

  # 癌种内所有可检验药物构成一个BH统计家族。
  results[, FDR := NA_real_]
  results[status == "ok", FDR := p.adjust(p_value, method = "BH"), by = Cancer]
  results[, more_sensitive_group := fifelse(
    cliffs_delta_C1_minus_C4 < 0, "C1",
    fifelse(cliffs_delta_C1_minus_C4 > 0, "C4", "tie")
  )]
  results[, significant := status == "ok" &
    !is.na(FDR) & FDR < fdr_cut &
    !is.na(cliffs_delta_C1_minus_C4) &
    abs(cliffs_delta_C1_minus_C4) >= delta_cut]

  results <- merge(results, annotation, by = "drug", all.x = TRUE, sort = FALSE)
  setorder(results, Cancer, drug)

  cross <- results[, .(
    n_available_cancers = uniqueN(Cancer),
    n_tested_cancers = uniqueN(Cancer[status == "ok"]),
    n_sig_cancers = uniqueN(Cancer[significant == TRUE]),
    n_C1_sensitive = uniqueN(Cancer[significant == TRUE & more_sensitive_group == "C1"]),
    n_C4_sensitive = uniqueN(Cancer[significant == TRUE & more_sensitive_group == "C4"]),
    sig_cancers_C1_sensitive = paste(sort(unique(Cancer[
      significant == TRUE & more_sensitive_group == "C1"
    ])), collapse = ";"),
    sig_cancers_C4_sensitive = paste(sort(unique(Cancer[
      significant == TRUE & more_sensitive_group == "C4"
    ])), collapse = ";"),
    min_FDR = if (any(status == "ok")) min(FDR[status == "ok"], na.rm = TRUE) else NA_real_,
    max_abs_delta_sig = if (any(significant == TRUE)) {
      max(abs(cliffs_delta_C1_minus_C4[significant == TRUE]), na.rm = TRUE)
    } else {
      NA_real_
    },
    median_abs_delta_sig = if (any(significant == TRUE)) {
      median(abs(cliffs_delta_C1_minus_C4[significant == TRUE]), na.rm = TRUE)
    } else {
      NA_real_
    }
  ), by = drug]
  cross[, context_dependent := n_C1_sensitive > 0L & n_C4_sensitive > 0L]
  cross[, recurrence_tier := fifelse(
    n_sig_cancers >= 5L, "HIGH_RECURRENCE_GE5",
    fifelse(
      n_sig_cancers >= 3L, "CROSS_CANCER_GE3",
      fifelse(n_sig_cancers >= 1L, "SINGLE_OR_LIMITED", "NO_SINGLE_CANCER_SIGNAL")
    )
  )]
  cross <- merge(cross, annotation, by = "drug", all.x = TRUE, sort = FALSE)
  setorder(cross, -n_sig_cancers, -max_abs_delta_sig, drug)

  parameters <- data.table(
    parameter = c(
      "formal_comparison",
      "min_group_n",
      "single_cancer_test",
      "within_cancer_fdr_method",
      "within_cancer_fdr_cut",
      "cliffs_delta_cut",
      "cross_cancer_rule",
      "prediction_direction",
      "plotting"
    ),
    value = c(
      "C1_vs_C4_only",
      as.character(min_group_n),
      "two_sided_Wilcoxon_rank_sum",
      "BH_across_all_ok_drugs_within_each_cancer",
      as.character(fdr_cut),
      as.character(delta_cut),
      "count_single_cancer_significant_results_only_no_second_BH_no_meta",
      "lower_predicted_value_means_more_sensitive",
      "not_run_in_this_step"
    )
  )

  qa <- rbindlist(list(
    data.table(
      check = "corrected_prediction_file_count_is_20",
      observed = as.character(nrow(file_audit)),
      expected = "20",
      pass = nrow(file_audit) == 20L
    ),
    data.table(
      check = "all_inputs_nonempty",
      observed = as.character(all(input_inventory$size_bytes > 0)),
      expected = "TRUE",
      pass = all(input_inventory$size_bytes > 0)
    ),
    data.table(
      check = "analysis_contains_only_C1_C4",
      observed = paste(sort(unique(matched_patients$C_subtype)), collapse = ";"),
      expected = "C1;C4",
      pass = identical(sort(unique(matched_patients$C_subtype)), c("C1", "C4"))
    ),
    data.table(
      check = "one_patient_per_cancer",
      observed = as.character(anyDuplicated(matched_patients[, paste(Cancer, Patient_ID)])),
      expected = "0",
      pass = anyDuplicated(matched_patients[, paste(Cancer, Patient_ID)]) == 0L
    ),
    data.table(
      check = "all_20_cancers_have_complete_drug_rows",
      observed = as.character(all(file_audit$n_drugs == file_audit$n_drugs[1L])),
      expected = "TRUE",
      pass = all(file_audit$n_drugs == file_audit$n_drugs[1L])
    ),
    data.table(
      check = "full_results_not_significant_only",
      observed = as.character(any(results$significant == FALSE)),
      expected = "TRUE",
      pass = any(results$significant == FALSE)
    ),
    data.table(
      check = "fdr_present_for_all_ok_rows",
      observed = as.character(all(!is.na(results$FDR[results$status == "ok"]))),
      expected = "TRUE",
      pass = all(!is.na(results$FDR[results$status == "ok"]))
    ),
    data.table(
      check = "cross_summary_covers_all_unique_drugs",
      observed = as.character(nrow(cross)),
      expected = as.character(uniqueN(results$drug)),
      pass = nrow(cross) == uniqueN(results$drug)
    ),
    data.table(
      check = "no_noneligible_cancer_ok_tests",
      observed = as.character(any(
        results$status == "ok" &
          results$Cancer %chin% coverage[formal_eligible == FALSE, Cancer]
      )),
      expected = "FALSE",
      pass = !any(
        results$status == "ok" &
          results$Cancer %chin% coverage[formal_eligible == FALSE, Cancer]
      )
    ),
    data.table(
      check = "no_figures_generated",
      observed = "TRUE",
      expected = "TRUE",
      pass = TRUE
    )
  ))

  write_tsv(file_audit, file.path(readiness_dir, "corrected_prediction_file_audit.tsv"))
  write_tsv(coverage, file.path(readiness_dir, "c1c4_drug_prediction_coverage_by_cancer.tsv"))
  write_tsv(matched_patients, file.path(readiness_dir, "c1c4_matched_patient_manifest.tsv"))
  write_tsv(parameters, file.path(readiness_dir, "stage5_frozen_parameters.tsv"))
  write_tsv(results, file.path(analysis_dir, "drug_c1_vs_c4_by_cancer_complete.tsv"))
  write_tsv(results[significant == TRUE], file.path(analysis_dir, "drug_c1_vs_c4_by_cancer_significant.tsv"))
  write_tsv(cross, file.path(analysis_dir, "drug_cross_cancer_significant_counts.tsv"))
  write_tsv(cross[n_sig_cancers >= 3L], file.path(analysis_dir, "drug_cross_cancer_candidates_ge3.tsv"))
  write_tsv(cross[n_sig_cancers >= 5L], file.path(analysis_dir, "drug_high_recurrence_candidates_ge5.tsv"))
  write_tsv(qa, file.path(analysis_dir, "stage5_drug_analysis_qa.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(analysis_dir, "stage5_sessionInfo.txt"))

  if (!all(qa$pass)) {
    stop("Stage 5 QA failed: ", paste(qa[pass == FALSE, check], collapse = ";"))
  }

  summary <- data.table(
    comparison = "C1_vs_C4",
    n_prediction_files = nrow(file_audit),
    n_formal_cancers = sum(coverage$formal_eligible),
    n_matched_patients = uniqueN(matched_patients$Patient_ID),
    n_unique_drugs = uniqueN(results$drug),
    n_complete_test_rows = sum(results$status == "ok"),
    n_single_cancer_significant_rows = sum(results$significant),
    n_drugs_with_any_signal = sum(cross$n_sig_cancers >= 1L),
    n_drugs_cross_cancer_ge3 = sum(cross$n_sig_cancers >= 3L),
    n_drugs_high_recurrence_ge5 = sum(cross$n_sig_cancers >= 5L),
    n_context_dependent_drugs = sum(cross$context_dependent),
    qa_pass = all(qa$pass),
    interpretation = "predicted_candidate_drug_vulnerability_only"
  )
  write_tsv(summary, file.path(analysis_dir, "stage5_drug_analysis_summary.tsv"))

  log_msg(
    "DONE", "stage5",
    sprintf(
      "formal_cancers=%d drugs=%d significant_rows=%d ge3=%d ge5=%d",
      summary$n_formal_cancers,
      summary$n_unique_drugs,
      summary$n_single_cancer_significant_rows,
      summary$n_drugs_cross_cancer_ge3,
      summary$n_drugs_high_recurrence_ge5
    )
  )
}, error = function(e) {
  log_msg("ERROR", "stage5", conditionMessage(e))
  quit(status = 1L)
})
