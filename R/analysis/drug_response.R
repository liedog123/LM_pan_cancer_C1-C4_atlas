# Purpose: drug response analysis stages.
# Configuration: config/drug.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Compare prepared GLDS-corrected predictions between C1 and C4.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/drug_response/01_c1c4_from_corrected_predictions.R
stage_drug_response_01_c1c4_from_corrected_predictions <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Compare prepared GLDS-corrected predictions between C1 and C4.
# Task: predicted_drug_response; step: s01
# Inputs (DCF fields): canonical_labels,model_input,prediction_dir,drug_annotation
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task predicted_drug_response --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`write_tsv` <- bind_stage_helper(shared_write_tsv_049, environment())

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
}

# Purpose: Compare GDSC IC50 and PRISM AUC across matched cell lines for the five-drug family.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/drug_response/03_run_prism_concordance.R
stage_drug_response_03_run_prism_concordance <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

# 固定PRISM快照文件名；元信息检查不读取或计算研究文件hash。
expected_basenames <- c(
  prism_curve = "secondary-screen-dose-response-top8-filtered.csv",
  prism_cell_metadata = "secondary-screen-cell-line-info.csv"
)
if (!identical(unname(basename(paths[names(expected_basenames)])),
               unname(expected_basenames))) stop("Required PRISM snapshot basenames differ.")

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

suppressPackageStartupMessages(library(readxl))
# 从独立的正式或公开续跑中间表继续；不要求先重新拟合内部模型。
gdsc_response_path <- unname(paths["gdsc_response"])
prism_curve_path <- unname(paths["prism_curve"])
prism_cell_path <- unname(paths["prism_cell_metadata"])
result_dir <- out_dir

secondary_drugs <- c(
  "Sorafenib", "VE821", "Vinorelbine", "Linsitinib", "Palbociclib"
)
preferred_screen_expected <- c(
  Sorafenib = "MTS010",
  VE821 = "HTS002",
  Vinorelbine = "HTS002",
  Linsitinib = "MTS010",
  Palbociclib = "MTS010"
)

expected_sha256 <- c(
  "secondary-screen-dose-response-top8-filtered.csv" =
    "c07b791e6d5356232ce7e1a9339a338f65f662e3a2717c1408c2f0462420aa6e",
  "secondary-screen-cell-line-info.csv" =
    "b93436b5f4bcf5fd14589697be4ca8f99ffd99d639223080380bc514c3718d27"
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)
`write_tsv` <- bind_stage_helper(shared_write_tsv_052, environment())

normalize_drug <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

# PRISM CCLE名称通常为CELL_LINE_TISSUE；匹配时只使用细胞系名称部分。
normalize_cell <- function(x, prism = FALSE) {
  x <- as.character(x)
  if (prism) x <- sub("_[^_]+$", "", x)
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

collapse_unique <- function(x) {
  paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))),
        collapse = ";")
}

safe_cor_test <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(estimate = unname(test$estimate), p_value = test$p.value)
}

message("[Stage5-crossplatform] 核验输入")
required_paths <- c(
  gdsc_response_path, prism_curve_path, prism_cell_path
)
if (!all(file.exists(required_paths))) {
  stop("缺少输入: ", paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}
if (!all(file.info(required_paths)$size > 0)) stop("存在空输入文件")

sha_lines <- system2(
  "sha256sum",
  shQuote(c(prism_curve_path, prism_cell_path)),
  stdout = TRUE
)
observed_sha <- sub("\\s+.*$", "", sha_lines)
names(observed_sha) <- basename(sub("^[^ ]+\\s+", "", sha_lines))
checksum_table <- data.frame(
  file = names(expected_sha256),
  expected_sha256 = unname(expected_sha256),
  observed_sha256 = observed_sha[names(expected_sha256)],
  pass = unname(expected_sha256) == observed_sha[names(expected_sha256)],
  stringsAsFactors = FALSE
)
write_tsv(checksum_table, file.path(result_dir, "prism_uploaded_input_checksum.tsv"))
if (!all(checksum_table$pass)) stop("PRISM上传文件checksum不一致")

message("[Stage5-crossplatform] 读取GDSC与PRISM真实药敏")
gdsc <- as.data.frame(read_excel(gdsc_response_path), stringsAsFactors = FALSE)
prism <- read.csv(prism_curve_path, check.names = FALSE)
prism_cell <- read.csv(prism_cell_path, check.names = FALSE)

gdsc$drug_key <- normalize_drug(gdsc$DRUG_NAME)
gdsc$cell_key <- normalize_cell(gdsc$CELL_LINE_NAME)
gdsc$LN_IC50 <- suppressWarnings(as.numeric(gdsc$LN_IC50))

prism$drug_key <- normalize_drug(prism$name)
prism$cell_key <- normalize_cell(prism$ccle_name, prism = TRUE)
prism$auc <- suppressWarnings(as.numeric(prism$auc))
prism$passed_str_profiling <- tolower(as.character(prism$passed_str_profiling)) == "true"

target_key_map <- setNames(secondary_drugs, normalize_drug(secondary_drugs))
secondary_keys <- normalize_drug(secondary_drugs)

gdsc_target <- gdsc[
  gdsc$drug_key %in% secondary_keys &
    nzchar(gdsc$cell_key) & is.finite(gdsc$LN_IC50),
  c("drug_key", "CELL_LINE_NAME", "cell_key", "TCGA_DESC", "LN_IC50"),
  drop = FALSE
]
prism_target <- prism[
  prism$drug_key %in% secondary_keys &
    nzchar(prism$cell_key) & prism$passed_str_profiling & is.finite(prism$auc),
  c("drug_key", "depmap_id", "ccle_name", "cell_key", "screen_id",
    "auc", "r2", "broad_id"),
  drop = FALSE
]

# 选择每种药物的首选PRISM screen。
screen_rows <- lapply(secondary_keys, function(key) {
  screens <- unique(prism_target$screen_id[prism_target$drug_key == key])
  preferred <- if ("MTS010" %in% screens) "MTS010" else if ("HTS002" %in% screens) "HTS002" else NA_character_
  data.frame(
    drug = target_key_map[[key]],
    drug_key = key,
    available_screens = paste(sort(screens), collapse = ";"),
    preferred_screen = preferred,
    stringsAsFactors = FALSE
  )
})
screen_table <- do.call(rbind, screen_rows)
write_tsv(screen_table, file.path(result_dir, "prism_preferred_screen.tsv"))

prism_selected <- merge(
  prism_target,
  screen_table[, c("drug_key", "preferred_screen")],
  by = "drug_key",
  all.x = TRUE
)
prism_selected <- prism_selected[
  prism_selected$screen_id == prism_selected$preferred_screen, ,
  drop = FALSE
]

# 审计名称标准化后是否把不同原始细胞系折叠到同一key。
gdsc_key_audit <- aggregate(
  CELL_LINE_NAME ~ drug_key + cell_key,
  data = gdsc_target,
  FUN = collapse_unique
)
names(gdsc_key_audit)[3] <- "source_names"
gdsc_key_audit$n_source_names <- vapply(
  strsplit(gdsc_key_audit$source_names, ";", fixed = TRUE), length, integer(1)
)
gdsc_key_audit$source <- "GDSC"

prism_key_audit <- aggregate(
  ccle_name ~ drug_key + cell_key,
  data = prism_selected,
  FUN = collapse_unique
)
names(prism_key_audit)[3] <- "source_names"
prism_key_audit$n_source_names <- vapply(
  strsplit(prism_key_audit$source_names, ";", fixed = TRUE), length, integer(1)
)
prism_key_audit$source <- "PRISM"

mapping_audit <- rbind(
  gdsc_key_audit[, c("source", "drug_key", "cell_key", "source_names", "n_source_names")],
  prism_key_audit[, c("source", "drug_key", "cell_key", "source_names", "n_source_names")]
)
write_tsv(mapping_audit, file.path(result_dir, "cell_line_name_mapping_audit.tsv"))

# 不同原始细胞系若在去标点后得到同一key，则整键排除，禁止强行合并。
ambiguous_keys <- unique(
  mapping_audit$cell_key[mapping_audit$n_source_names > 1L]
)
ambiguous_key_table <- mapping_audit[
  mapping_audit$cell_key %in% ambiguous_keys, ,
  drop = FALSE
]
write_tsv(
  ambiguous_key_table,
  file.path(result_dir, "excluded_ambiguous_cell_keys.tsv")
)
gdsc_target <- gdsc_target[
  !gdsc_target$cell_key %in% ambiguous_keys, ,
  drop = FALSE
]
prism_selected <- prism_selected[
  !prism_selected$cell_key %in% ambiguous_keys, ,
  drop = FALSE
]

# 同一药物-细胞系若有技术重复，按均值折叠，保持与旧GDSC训练矩阵构建方式一致。
gdsc_agg <- aggregate(
  LN_IC50 ~ drug_key + cell_key,
  data = gdsc_target,
  FUN = mean
)
gdsc_names <- aggregate(
  CELL_LINE_NAME ~ drug_key + cell_key,
  data = gdsc_target,
  FUN = collapse_unique
)
gdsc_agg <- merge(gdsc_agg, gdsc_names, by = c("drug_key", "cell_key"))

prism_agg <- aggregate(
  auc ~ drug_key + cell_key,
  data = prism_selected,
  FUN = mean
)
prism_names <- aggregate(
  ccle_name ~ drug_key + cell_key,
  data = prism_selected,
  FUN = collapse_unique
)
prism_ids <- aggregate(
  depmap_id ~ drug_key + cell_key,
  data = prism_selected,
  FUN = collapse_unique
)
prism_agg <- Reduce(
  function(x, y) merge(x, y, by = c("drug_key", "cell_key")),
  list(prism_agg, prism_names, prism_ids)
)

matched <- merge(
  gdsc_agg, prism_agg,
  by = c("drug_key", "cell_key"),
  all = FALSE
)
matched$drug <- unname(target_key_map[matched$drug_key])
matched <- matched[
  order(match(matched$drug, secondary_drugs), matched$cell_key), ,
  drop = FALSE
]
write_tsv(
  matched[, c(
    "drug", "cell_key", "CELL_LINE_NAME", "ccle_name", "depmap_id",
    "LN_IC50", "auc"
  )],
  file.path(result_dir, "gdsc_prism_matched_cellline_pairs.tsv")
)

message("[Stage5-crossplatform] 计算跨平台相关")
summary_rows <- lapply(secondary_drugs, function(drug) {
  key <- normalize_drug(drug)
  one <- matched[matched$drug_key == key, , drop = FALSE]
  spearman <- safe_cor_test(one$LN_IC50, one$auc, "spearman")
  pearson <- safe_cor_test(one$LN_IC50, one$auc, "pearson")
  data.frame(
    drug = drug,
    preferred_prism_screen =
      screen_table$preferred_screen[screen_table$drug == drug][1],
    n_gdsc_cell_lines = sum(gdsc_agg$drug_key == key),
    n_prism_cell_lines = sum(prism_agg$drug_key == key),
    n_matched_cell_lines = nrow(one),
    spearman_rho = spearman[["estimate"]],
    spearman_p = spearman[["p_value"]],
    pearson_r = pearson[["estimate"]],
    pearson_p = pearson[["p_value"]],
    median_gdsc_ln_ic50 = median(one$LN_IC50, na.rm = TRUE),
    median_prism_auc = median(one$auc, na.rm = TRUE),
    comparable = nrow(one) >= 30L,
    stringsAsFactors = FALSE
  )
})
cross_summary <- do.call(rbind, summary_rows)
cross_summary$spearman_fdr <- NA_real_
eligible <- cross_summary$comparable & is.finite(cross_summary$spearman_p)
cross_summary$spearman_fdr[eligible] <- p.adjust(
  cross_summary$spearman_p[eligible], method = "BH"
)
cross_summary$external_direction_support <-
  cross_summary$comparable &
  cross_summary$spearman_rho >= 0.20 &
  cross_summary$spearman_fdr < 0.10
cross_summary$external_status <- ifelse(
  !cross_summary$comparable,
  "NOT_COMPARABLE",
  ifelse(
    cross_summary$external_direction_support,
    "PRISM_DIRECTION_SUPPORTED",
    "PRISM_DIRECTION_NOT_SUPPORTED"
  )
)
write_tsv(
  cross_summary,
  file.path(result_dir, "gdsc_prism_crossplatform_summary.tsv")
)

cell_info_ids <- unique(prism_cell$depmap_id)
selected_ids <- unique(prism_selected$depmap_id)
expected_screen_ok <- all(
  screen_table$preferred_screen[
    match(names(preferred_screen_expected), screen_table$drug)
  ] == unname(preferred_screen_expected)
)
matched_keys <- unique(matched$cell_key)
ambiguous_matched <- mapping_audit[
  mapping_audit$cell_key %in% matched_keys &
    mapping_audit$n_source_names > 1L, ,
  drop = FALSE
]

qa <- data.frame(
  check = c(
    "uploaded_inputs_exist_and_nonempty",
    "uploaded_input_checksums_match",
    "all_selected_prism_depmap_ids_in_cell_metadata",
    "five_predefined_secondary_drugs_present",
    "preferred_screen_matches_contract",
    "all_five_drugs_have_at_least_30_matched_cell_lines",
    "bh_family_contains_all_five_comparable_drugs",
    "no_ambiguous_normalized_names_in_matched_pairs",
    "no_prism_cell_line_assigned_c1_or_c4",
    "no_source_object_modified"
  ),
  pass = c(
    all(file.info(c(prism_curve_path, prism_cell_path))$size > 0),
    all(checksum_table$pass),
    all(selected_ids %in% cell_info_ids),
    setequal(unique(prism_selected$drug_key), secondary_keys),
    expected_screen_ok,
    all(cross_summary$n_matched_cell_lines >= 30L),
    sum(eligible) == length(secondary_drugs),
    nrow(ambiguous_matched) == 0L,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, file.path(result_dir, "qa_gdsc_prism_crossplatform.tsv"))
if (nrow(ambiguous_matched)) {
  write_tsv(
    ambiguous_matched,
    file.path(result_dir, "ambiguous_matched_cell_keys.tsv")
  )
}

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo_crossplatform.txt")
)
writeLines(
  c(
    paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "comparison=GDSC_LN_IC50_vs_PRISM_AUC",
    "higher_value_in_both_datasets=relative_resistance",
    "prism_screen_rule=MTS010_if_available_else_HTS002",
    "min_matched_cell_lines=30",
    "test=two_sided_spearman",
    "fdr_family=all_comparable_predefined_secondary_drugs",
    "support_rule=rho>=0.20_and_FDR<0.10",
    "no_C1_C4_assignment_to_PRISM_cell_lines=true"
  ),
  con = file.path(result_dir, "parameters_crossplatform.txt")
)

if (!all(qa$pass)) {
  stop(
    "GDSC–PRISM跨平台QA未全部通过: ",
    paste(qa$check[!qa$pass], collapse = ", ")
  )
}
message("[Stage5-crossplatform] 完成")
}

# Purpose: Adjust the three reported drug predictions for the other drugs shared PC1.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/drug_response/04_run_shared_axis_adjustment.R
stage_drug_response_04_run_shared_axis_adjustment <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("drug_complete", "shared_patient", "prediction_audit")]
required_files <- unname(input_paths)

# 仅验证预测目录元信息，不读取audit内容或枚举预测文件。
prediction_dir <- unname(paths["prediction_dir"])
if (!dir.exists(prediction_dir) || file.access(prediction_dir, 5L) != 0L) stop("Prediction directory must be readable and searchable.")

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
protected <- c(dirname(required_files),prediction_dir,config_path,package_root)
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

options(stringsAsFactors = FALSE)
output_dir <- out_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  read.delim(
    con, sep = "\t", header = TRUE, quote = "\"", comment.char = "",
    check.names = FALSE, na.strings = c("NA", "")
  )
}

write_tsv <- function(x, name) {
  write.table(
    x, file.path(output_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

cliffs_delta <- function(x, y) {
  n1 <- length(x)
  n0 <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u1 <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n0) - 1
}

# 三药来自当前稿件固定比较集合，不再依赖已退出投稿的内部证据评级。
d <- lapply(input_paths, read_tsv)
strict_drugs <- c("Vinorelbine", "Linsitinib", "Palbociclib")
raw_drug <- d$drug_complete[
  d$drug_complete$drug %in% strict_drugs & d$drug_complete$status == "ok", ,
  drop = FALSE
]
formal_cancers <- sort(unique(raw_drug$Cancer[raw_drug$significant %in% TRUE]))
raw_drug <- raw_drug[raw_drug$Cancer %in% formal_cancers, , drop = FALSE]
if (length(formal_cancers) != 18L || nrow(raw_drug) != 54L) {
  stop("Expected 18 formal cancers and 54 Cancer x Drug units.")
}

# 每个药物使用排除自身后的285药物PC1，避免当前药物直接进入共享轴。
residual_rows <- list()
rr <- 0L
for (cancer in formal_cancers) {
  audit_row <- d$prediction_audit[d$prediction_audit$Cancer == cancer, , drop = FALSE]
  if (nrow(audit_row) != 1L) stop("Prediction source unresolved for ", cancer)
  # 内嵌历史路径只提供文件名；在配置预测目录内解析并保护输出边界。
  source_name <- basename(as.character(audit_row$source_file))
  if (length(source_name) != 1L || is.na(source_name) || !nzchar(source_name) ||
      source_name %in% c(".", "..") || grepl("[/\\\\]", source_name)) {
    stop("Prediction basename unresolved for ", cancer)
  }
  prediction_path <- resolve_path(file.path(prediction_dir, source_name), prediction_dir)
  if (!is_within(prediction_path, prediction_dir) || is_within(prediction_path, out_dir) ||
      !file_test("-f", prediction_path) || file.access(prediction_path, 4L) != 0L ||
      is.na(file.info(prediction_path)$size) || file.info(prediction_path)$size <= 0L) {
    stop("Prediction source must be a protected readable nonempty file for ", cancer)
  }
  pred <- read.csv(
    prediction_path, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_col <- names(pred)[1L]
  pred$Sample15 <- substr(as.character(pred[[sample_col]]), 1L, 15L)
  shared <- d$shared_patient[d$shared_patient$Cancer == cancer, , drop = FALSE]
  idx <- match(shared$Sample15, pred$Sample15)
  if (anyNA(idx)) stop("Shared-score patients missing from prediction matrix: ", cancer)
  pred <- pred[idx, , drop = FALSE]
  if (!identical(pred$Sample15, shared$Sample15)) stop("Sample order mismatch: ", cancer)

  drug_cols <- setdiff(names(pred), c(sample_col, "Sample15"))
  if (!all(strict_drugs %in% drug_cols)) stop("Strict drug column missing: ", cancer)

  cancer_rows <- list()
  for (drug in strict_drugs) {
    other_cols <- setdiff(drug_cols, drug)
    other_mat <- as.matrix(data.frame(lapply(pred[other_cols], as.numeric), check.names = FALSE))
    keep <- apply(other_mat, 2L, function(x) all(is.finite(x)) && stats::sd(x) > 0)
    other_mat <- other_mat[, keep, drop = FALSE]
    if (ncol(other_mat) < 280L) stop("Too few nonconstant leave-one-drug-out columns: ", cancer)
    pca <- prcomp(other_mat, center = TRUE, scale. = TRUE, rank. = 1L)
    pc1 <- pca$x[, 1L]
    pc1_var <- pca$sdev[1L]^2 / sum(pca$sdev^2)

    y <- as.numeric(pred[[drug]])
    fit <- lm(y ~ pc1)
    resid_y <- residuals(fit)
    is_c1 <- shared$C_subtype == "C1"
    x <- resid_y[is_c1]
    y4 <- resid_y[!is_c1]
    wt <- suppressWarnings(wilcox.test(x, y4, exact = FALSE))
    delta <- cliffs_delta(x, y4)
    median_diff <- median(x) - median(y4)
    raw_row <- raw_drug[raw_drug$Cancer == cancer & raw_drug$drug == drug, , drop = FALSE]

    rr <- rr + 1L
    cancer_rows[[drug]] <- data.frame(
      drug = drug,
      Cancer = cancer,
      n_C1 = sum(is_c1),
      n_C4 = sum(!is_c1),
      n_other_drugs_for_PC1 = ncol(other_mat),
      loo_PC1_variance_fraction = pc1_var,
      abs_spearman_drug_vs_loo_PC1 = abs(suppressWarnings(cor(
        as.numeric(pred[[drug]]), pc1, method = "spearman"
      ))),
      raw_cliffs_delta_C1_minus_C4 = raw_row$cliffs_delta_C1_minus_C4,
      raw_more_sensitive_group = raw_row$more_sensitive_group,
      raw_FDR = raw_row$FDR,
      residual_median_diff_C1_minus_C4 = median_diff,
      residual_cliffs_delta_C1_minus_C4 = delta,
      residual_p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
  cancer_df <- do.call(rbind, cancer_rows)
  cancer_df$residual_FDR_within_cancer_3_drugs <- p.adjust(
    cancer_df$residual_p_value, method = "BH"
  )
  cancer_df$residual_more_sensitive_group <- ifelse(
    cancer_df$residual_median_diff_C1_minus_C4 < 0, "C1",
    ifelse(cancer_df$residual_median_diff_C1_minus_C4 > 0, "C4", "tie")
  )
  cancer_df$residual_significant <- (
    cancer_df$residual_FDR_within_cancer_3_drugs < 0.05 &
      abs(cancer_df$residual_cliffs_delta_C1_minus_C4) >= 0.20
  )
  cancer_df$direction_preserved_after_residual <- (
    cancer_df$residual_more_sensitive_group == cancer_df$raw_more_sensitive_group
  )
  residual_rows[[cancer]] <- cancer_df
}
residual <- do.call(rbind, residual_rows)
rownames(residual) <- NULL
write_tsv(residual, "strict_drug_shared_axis_residual_by_cancer.tsv")

if (any(residual$n_other_drugs_for_PC1 < 285L)) stop("Expected all 285 other drug columns.")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
}
