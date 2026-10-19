# Purpose: genomics analysis stages.
# Configuration: config/genomics.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Compare patient-level somatic mutation frequencies between C1 and C4.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/genomics/01_run_mutation_c1c4.R
stage_genomics_01_run_mutation_c1c4 <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
callable_dir <- unname(paths["callable_dir"])
if(!dir.exists(callable_dir) || file.access(callable_dir,5L)!=0L) stop("Callable input must be a readable searchable directory.")
callable_names <- c("mutation_callable_decision_v0.1.tsv", "mutation_callable_recovery_qa.tsv",
 "mutation_callable_patient_audit.tsv", "mutation_callable_formal_eligibility_by_cancer.tsv")
required_files <- c(file.path(callable_dir, callable_names), unname(paths["mutation_events"]))
required_files <- vapply(required_files,resolve_path,character(1L),base_dir=dirname(config_path))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
protected <- c(callable_dir,dirname(required_files),config_path,package_root)
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

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))
xena_path <- unname(paths["mutation_events"])
`log_msg` <- bind_stage_helper(shared_log_msg_016, environment())
`assert_nonempty` <- bind_stage_helper(shared_assert_nonempty_021, environment())
`write_tsv` <- bind_stage_helper(shared_write_tsv_049, environment())
`save_session_info` <- bind_stage_helper(shared_save_session_info_019, environment())

nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Translation_Start_Site",
  "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"
)

fisher_gene <- function(a, n1, c, n4) {
  ft <- fisher.test(
    matrix(c(a, n1 - a, c, n4 - c), nrow = 2L, byrow = TRUE),
    alternative = "two.sided"
  )
  c(
    p_value = unname(ft$p.value),
    odds_ratio_C1_vs_C4 = if (length(ft$estimate)) unname(ft$estimate) else NA_real_,
    odds_ratio_ci_low = unname(ft$conf.int[1L]),
    odds_ratio_ci_high = unname(ft$conf.int[2L])
  )
}

tryCatch({
  log_msg("START", "mutation_c1c4", "Formal within-cancer C1 vs C4 mutation analysis")
  required <- c(
    file.path(callable_dir, "mutation_callable_decision_v0.1.tsv"),
    file.path(callable_dir, "mutation_callable_recovery_qa.tsv"),
    file.path(callable_dir, "mutation_callable_patient_audit.tsv"),
    file.path(callable_dir, "mutation_callable_formal_eligibility_by_cancer.tsv"),
    xena_path
  )
  invisible(lapply(required, assert_nonempty))

  decision <- fread(required[1L])
  recovery_qa <- fread(required[2L])
  if (!all(recovery_qa$pass)) stop("Callable recovery QA did not pass")
  if (!isTRUE(decision$accepted_for_formal_denominator[1L])) {
    stop("Callable universe was not accepted for formal mutation analysis")
  }

  callable_audit <- fread(required[3L])
  callable <- callable_audit[
    mc3_callable_exact_sample15 == TRUE,
    .(Sample15 = selected_canonical_sample15, Patient_ID, Cancer, C_subtype)
  ]
  if (anyDuplicated(callable$Patient_ID)) stop("Duplicate callable patients")
  if (any(!callable$C_subtype %chin% c("C1", "C4"))) stop("Non-C1/C4 subtype")

  write_tsv(callable, file.path(out_dir, "mutation_callable_patient_manifest.tsv"))

  log_msg("CHECKPOINT", "event_read", "Reading frozen Xena PASS mutation events")
  events <- fread(xena_path, select = c("sample", "gene", "effect"))
  events[, Sample15 := substr(as.character(sample), 1L, 15L)]
  events <- events[
    effect %chin% nonsilent_classes &
      !is.na(gene) & gene != "" &
      Sample15 %chin% callable$Sample15,
    .(Sample15, Gene = as.character(gene), Variant_Classification = effect)
  ]
  event_class_counts <- events[, .(
    event_rows = .N,
    unique_sample_gene_events = uniqueN(paste(Sample15, Gene, sep = "|")),
    unique_samples = uniqueN(Sample15),
    unique_genes = uniqueN(Gene)
  ), by = Variant_Classification][order(-event_rows)]
  write_tsv(
    event_class_counts,
    file.path(out_dir, "mutation_nonsilent_event_class_counts.tsv")
  )

  sparse_binary <- unique(events[, .(Sample15, Gene)])
  sparse_binary[, mutated := 1L]
  sparse_binary <- merge(
    sparse_binary,
    callable,
    by = "Sample15",
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(sparse_binary$Patient_ID))) {
    stop("Mutation event outside accepted callable manifest")
  }
  setorder(sparse_binary, Cancer, Sample15, Gene)
  write_tsv(
    sparse_binary[, .(
      Sample15, Patient_ID, Cancer, C_subtype, Gene, mutated
    )],
    file.path(out_dir, "mutation_nonsilent_patient_gene_sparse_binary.tsv")
  )
  rm(events)
  gc()

  cancer_coverage <- callable[, .N, by = .(Cancer, C_subtype)]
  cancer_coverage <- dcast(
    cancer_coverage,
    Cancer ~ C_subtype,
    value.var = "N",
    fill = 0L
  )
  if (!"C1" %chin% names(cancer_coverage)) cancer_coverage[, C1 := 0L]
  if (!"C4" %chin% names(cancer_coverage)) cancer_coverage[, C4 := 0L]
  cancer_coverage[, `:=`(
    min_group_n = 30L,
    formal_eligible = C1 >= 30L & C4 >= 30L
  )]
  setorder(cancer_coverage, Cancer)
  write_tsv(
    cancer_coverage,
    file.path(out_dir, "mutation_formal_coverage_by_cancer.tsv")
  )

  eligible_cancers <- cancer_coverage[formal_eligible == TRUE, Cancer]
  result_list <- vector("list", length(eligible_cancers))

  for (i in seq_along(eligible_cancers)) {
    cancer <- eligible_cancers[i]
    group_dat <- callable[Cancer == cancer]
    n1 <- group_dat[C_subtype == "C1", .N]
    n4 <- group_dat[C_subtype == "C4", .N]

    event_dat <- sparse_binary[Cancer == cancer]
    counts <- event_dat[, .(
      mutated_patients = uniqueN(Sample15)
    ), by = .(Gene, C_subtype)]
    counts_wide <- dcast(
      counts,
      Gene ~ C_subtype,
      value.var = "mutated_patients",
      fill = 0L
    )
    if (!"C1" %chin% names(counts_wide)) counts_wide[, C1 := 0L]
    if (!"C4" %chin% names(counts_wide)) counts_wide[, C4 := 0L]
    counts_wide <- counts_wide[C1 + C4 >= 5L]

    if (nrow(counts_wide) == 0L) {
      log_msg("CHECKPOINT", "cancer", paste(cancer, "no gene passed event threshold"))
      next
    }

    test_values <- t(vapply(
      seq_len(nrow(counts_wide)),
      function(j) fisher_gene(counts_wide$C1[j], n1, counts_wide$C4[j], n4),
      numeric(4L)
    ))
    results <- cbind(
      data.table(
        Cancer = cancer,
        Gene = counts_wide$Gene,
        n_C1 = n1,
        n_C4 = n4,
        mutated_C1 = counts_wide$C1,
        mutated_C4 = counts_wide$C4
      ),
      as.data.table(test_values)
    )
    results[, `:=`(
      frequency_C1 = mutated_C1 / n_C1,
      frequency_C4 = mutated_C4 / n_C4
    )]
    results[, frequency_difference_C1_minus_C4 := frequency_C1 - frequency_C4]
    results[, abs_frequency_difference := abs(frequency_difference_C1_minus_C4)]
    results[, FDR_within_cancer := p.adjust(p_value, method = "BH")]
    results[, single_cancer_significant :=
      FDR_within_cancer < 0.05 &
        abs(frequency_difference_C1_minus_C4) >= 0.10
    ]
    results[, direction := fifelse(
      frequency_difference_C1_minus_C4 > 0,
      "C1_higher",
      fifelse(frequency_difference_C1_minus_C4 < 0, "C4_higher", "tie")
    )]
    setorder(results, FDR_within_cancer, -abs_frequency_difference)
    result_list[[i]] <- results
    log_msg(
      "CHECKPOINT", "cancer",
      paste(
        cancer,
        "n_C1:", n1,
        "n_C4:", n4,
        "tested genes:", nrow(results),
        "significant:", sum(results$single_cancer_significant)
      )
    )
  }

  mutation_results <- rbindlist(result_list, fill = TRUE)
  if (nrow(mutation_results) == 0L) stop("No mutation gene tests were generated")
  write_tsv(
    mutation_results,
    file.path(out_dir, "mutation_gene_c1_vs_c4_by_cancer.tsv")
  )

  cross_cancer <- mutation_results[, .(
    n_tested_cancers = uniqueN(Cancer),
    n_significant_cancers = sum(single_cancer_significant),
    n_significant_C1_higher = sum(
      single_cancer_significant & direction == "C1_higher"
    ),
    n_significant_C4_higher = sum(
      single_cancer_significant & direction == "C4_higher"
    ),
    significant_cancers_C1_higher = paste(
      sort(Cancer[single_cancer_significant & direction == "C1_higher"]),
      collapse = ";"
    ),
    significant_cancers_C4_higher = paste(
      sort(Cancer[single_cancer_significant & direction == "C4_higher"]),
      collapse = ";"
    ),
    max_abs_frequency_difference = max(
      abs(frequency_difference_C1_minus_C4),
      na.rm = TRUE
    ),
    min_within_cancer_FDR = min(FDR_within_cancer, na.rm = TRUE)
  ), by = Gene]
  cross_cancer[, context_dependent := n_significant_C1_higher > 0L &
    n_significant_C4_higher > 0L
  ]
  cross_cancer[, dominant_significant_count := pmax(
    n_significant_C1_higher,
    n_significant_C4_higher
  )]
  setorder(
    cross_cancer,
    -n_significant_cancers,
    -dominant_significant_count,
    min_within_cancer_FDR,
    Gene
  )
  write_tsv(
    cross_cancer,
    file.path(out_dir, "mutation_gene_cross_cancer_significant_counts.tsv")
  )
  write_tsv(
    cross_cancer[n_significant_cancers > 0L],
    file.path(out_dir, "mutation_single_cancer_significant_gene_summary.tsv")
  )

  significant_units <- mutation_results[single_cancer_significant == TRUE]
  write_tsv(
    significant_units,
    file.path(out_dir, "mutation_significant_gene_cancer_units.tsv")
  )

  endpoint_summary <- data.table(
    endpoint = "nonsilent_gene_mutation",
    n_callable_patients = nrow(callable),
    n_formal_cancers = length(eligible_cancers),
    n_gene_cancer_tests = nrow(mutation_results),
    n_significant_gene_cancer_units = nrow(significant_units),
    n_genes_significant_in_at_least_one_cancer =
      cross_cancer[n_significant_cancers > 0L, .N],
    n_genes_recurrent_same_direction_ge2_cancers = cross_cancer[
      pmax(n_significant_C1_higher, n_significant_C4_higher) >= 2L,
      .N
    ],
    n_context_dependent_genes = cross_cancer[context_dependent == TRUE, .N]
  )
  write_tsv(
    endpoint_summary,
    file.path(out_dir, "mutation_endpoint_summary.tsv")
  )

  qa <- data.table(
    check = c(
      "callable_recovery_qa_pass",
      "only_exact_callable_samples_used",
      "callable_patient_unique",
      "only_c1_c4",
      "formal_cancers_meet_group_n_30",
      "all_tested_genes_meet_total_mutated_n_5",
      "bh_is_within_cancer",
      "effect_threshold_applied_to_significance",
      "no_cross_cancer_bh",
      "no_figure_output_generated"
    ),
    pass = c(
      all(recovery_qa$pass),
      all(sparse_binary$Sample15 %chin% callable$Sample15),
      !anyDuplicated(callable$Patient_ID),
      all(callable$C_subtype %chin% c("C1", "C4")),
      cancer_coverage[formal_eligible == TRUE, all(C1 >= 30L & C4 >= 30L)],
      mutation_results[, all(mutated_C1 + mutated_C4 >= 5L)],
      mutation_results[, all(
        abs(
          FDR_within_cancer -
            ave(p_value, Cancer, FUN = function(z) p.adjust(z, method = "BH"))
        ) < 1e-12
      )],
      mutation_results[single_cancer_significant == TRUE, all(
        FDR_within_cancer < 0.05 &
          abs(frequency_difference_C1_minus_C4) >= 0.10
      )],
      !any(c("cross_cancer_FDR", "meta_FDR") %chin% names(cross_cancer)),
      !dir.exists(file.path(out_dir, "figures/01_mutation"))
    )
  )
  write_tsv(qa, file.path(out_dir, "mutation_analysis_qa.tsv"))
  if (!all(qa$pass)) {
    stop(
      "Mutation analysis QA failed: ",
      paste(qa[pass == FALSE, check], collapse = ", ")
    )
  }

  save_session_info(file.path(out_dir, "mutation_analysis_sessionInfo.txt"))
  log_msg(
    "DONE", "mutation_c1c4",
    paste(
      "Callable patients:", nrow(callable),
      "formal cancers:", length(eligible_cancers),
      "gene-cancer tests:", nrow(mutation_results),
      "significant units:", nrow(significant_units)
    )
  )
}, error = function(e) {
  log_msg("ERROR", "mutation_c1c4", conditionMessage(e))
  quit(status = 1L)
})
}

# Purpose: Compare mutation burden and MSI measures using the adopted eligibility rules.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/genomics/02_run_tmb_msi_c1c4.R
stage_genomics_02_run_tmb_msi_c1c4 <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
readiness_dir <- unname(paths["readiness_dir"])
if(!dir.exists(readiness_dir) || file.access(readiness_dir,5L)!=0L) stop("Readiness input must be a readable searchable directory.")
required_names <- c("tmb_mantis_c1c4_patient_data.tsv", "readiness_qa.tsv",
 "stage4_layer_execution_status.tsv")
required_files <- vapply(file.path(readiness_dir,required_names),resolve_path,character(1L),base_dir=dirname(config_path))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
protected <- c(readiness_dir,dirname(required_files),config_path,package_root)
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

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# 普通执行使用独立输出根；源诊断图不等于冻结Fig. 4重绘。
output_root <- out_dir
out_dir <- file.path(output_root, "results")
fig_dir <- file.path(output_root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
`log_msg` <- bind_stage_helper(shared_log_msg_016, environment())
`assert_nonempty` <- bind_stage_helper(shared_assert_nonempty_021, environment())
`assert_cols` <- bind_stage_helper(shared_assert_cols_054, environment())
`write_tsv` <- bind_stage_helper(shared_write_tsv_049, environment())
`save_session_info` <- bind_stage_helper(shared_save_session_info_019, environment())
`wilcox_summary` <- bind_stage_helper(shared_wilcox_summary_056, environment())
`summarize_cross_cancer` <- bind_stage_helper(shared_summarize_cross_cancer_033, environment())
`save_figure_pair` <- bind_stage_helper(shared_save_figure_pair_028, environment())


continuous_endpoint <- function(dat, value_col, endpoint) {
  cancers <- sort(unique(dat$Cancer))
  results <- vector("list", length(cancers))
  for (i in seq_along(cancers)) {
    cancer <- cancers[i]
    z <- dat[Cancer == cancer & !is.na(get(value_col))]
    x <- z[C_subtype == "C1", get(value_col)]
    y <- z[C_subtype == "C4", get(value_col)]
    eligible <- length(x) >= 30L && length(y) >= 30L

    base <- data.table(
      endpoint = endpoint,
      Cancer = cancer,
      n_C1 = length(x),
      n_C4 = length(y),
      median_C1 = if (length(x)) median(x) else NA_real_,
      q1_C1 = if (length(x)) unname(quantile(x, 0.25)) else NA_real_,
      q3_C1 = if (length(x)) unname(quantile(x, 0.75)) else NA_real_,
      median_C4 = if (length(y)) median(y) else NA_real_,
      q1_C4 = if (length(y)) unname(quantile(y, 0.25)) else NA_real_,
      q3_C4 = if (length(y)) unname(quantile(y, 0.75)) else NA_real_,
      formal_eligible = eligible
    )
    if (eligible) {
      stats <- wilcox_summary(x, y)
      base[, names(stats) := as.list(stats)]
    } else {
      base[, `:=`(
        p_value = NA_real_,
        hl_shift_c1_minus_c4 = NA_real_,
        hl_ci_low = NA_real_,
        hl_ci_high = NA_real_,
        cliff_delta_c1_minus_c4 = NA_real_
      )]
    }
    base[, `:=`(
      effect_c1_minus_c4 = hl_shift_c1_minus_c4,
      significant = formal_eligible & !is.na(p_value) & p_value < 0.05,
      direction = fifelse(
        hl_shift_c1_minus_c4 > 0, "C1_higher",
        fifelse(hl_shift_c1_minus_c4 < 0, "C4_higher", "tie")
      )
    )]
    results[[i]] <- base
  }
  rbindlist(results, fill = TRUE)
}

make_continuous_forest <- function(x, title, x_label) {
  pdat <- x[formal_eligible == TRUE]
  pdat[, Cancer := factor(Cancer, levels = rev(Cancer[order(effect_c1_minus_c4)]))]
  ggplot(pdat, aes(
    x = effect_c1_minus_c4,
    y = Cancer,
    xmin = hl_ci_low,
    xmax = hl_ci_high,
    colour = significant
  )) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.45) +
    geom_errorbarh(height = 0.18, linewidth = 0.55) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D")) +
    labs(title = title, x = x_label, y = NULL, colour = "P < 0.05") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 11)
    )
}

tryCatch({
  log_msg("START", "tmb_msi", "Formal C1 vs C4 scalar genomic-instability endpoints")
  required <- c(
    file.path(readiness_dir, "tmb_mantis_c1c4_patient_data.tsv"),
    file.path(readiness_dir, "readiness_qa.tsv"),
    file.path(readiness_dir, "stage4_layer_execution_status.tsv")
  )
  invisible(lapply(required, assert_nonempty))
  qa <- fread(required[2L])
  if (!all(qa$pass)) stop("Readiness QA did not pass")

  dat <- fread(required[1L])
  assert_cols(
    dat,
    c("Patient_ID", "Cancer", "C_subtype", "TMB_Raw", "TMB_Log10", "MANTIS_Score", "MSI_Class"),
    "TMB/MANTIS patient data"
  )
  if (anyDuplicated(dat$Patient_ID)) stop("Duplicate patients in TMB/MANTIS input")
  if (any(!dat$C_subtype %chin% c("C1", "C4"))) stop("Non-C1/C4 labels in formal input")
  if (dat[!is.na(TMB_Raw), any(TMB_Raw < 0)]) stop("Negative TMB values")
  if (dat[!is.na(MANTIS_Score), any(!is.finite(MANTIS_Score))]) stop("Non-finite MANTIS values")

  tmb_results <- continuous_endpoint(dat, "TMB_Log10", "TMB_log10")
  mantis_results <- continuous_endpoint(dat, "MANTIS_Score", "MANTIS_continuous")
  write_tsv(tmb_results, file.path(out_dir, "tmb_log10_c1_vs_c4_by_cancer.tsv"))
  write_tsv(mantis_results, file.path(out_dir, "mantis_continuous_c1_vs_c4_by_cancer.tsv"))

  tmb_cross <- summarize_cross_cancer(
    tmb_results[formal_eligible == TRUE],
    "endpoint"
  )
  mantis_cross <- summarize_cross_cancer(
    mantis_results[formal_eligible == TRUE],
    "endpoint"
  )
  write_tsv(tmb_cross, file.path(out_dir, "tmb_cross_cancer_significant_counts.tsv"))
  write_tsv(mantis_cross, file.path(out_dir, "mantis_cross_cancer_significant_counts.tsv"))
  log_msg("CHECKPOINT", "continuous", "TMB and MANTIS continuous tests complete")

  msi_results <- dat[!is.na(MSI_Class), {
    n1 <- sum(C_subtype == "C1")
    n4 <- sum(C_subtype == "C4")
    h1 <- sum(C_subtype == "C1" & MSI_Class == "MSI-H")
    h4 <- sum(C_subtype == "C4" & MSI_Class == "MSI-H")
    eligible_n <- n1 >= 30L && n4 >= 30L
    eligible_event <- h1 + h4 >= 5L
    formal <- eligible_n && eligible_event
    if (formal) {
      ft <- fisher.test(matrix(c(h1, n1 - h1, h4, n4 - h4), nrow = 2L, byrow = TRUE))
      p <- unname(ft$p.value)
      odds <- unname(ft$estimate)
      ci_low <- unname(ft$conf.int[1L])
      ci_high <- unname(ft$conf.int[2L])
    } else {
      p <- odds <- ci_low <- ci_high <- NA_real_
    }
    .(
      endpoint = "MSI_binary",
      n_C1 = n1,
      n_C4 = n4,
      MSI_H_C1 = h1,
      MSI_H_C4 = h4,
      MSI_H_frequency_C1 = h1 / n1,
      MSI_H_frequency_C4 = h4 / n4,
      effect_c1_minus_c4 = h1 / n1 - h4 / n4,
      odds_ratio_C1_vs_C4 = odds,
      odds_ratio_ci_low = ci_low,
      odds_ratio_ci_high = ci_high,
      p_value = p,
      sample_size_eligible = eligible_n,
      event_count_eligible = eligible_event,
      formal_eligible = formal,
      status = if (formal) "TESTED" else if (!eligible_n) "HOLD_MIN_GROUP_N" else "HOLD_RARE_EVENT"
    )
  }, by = Cancer]
  msi_results[, significant := formal_eligible & !is.na(p_value) & p_value < 0.05]
  msi_results[, direction := fifelse(
    effect_c1_minus_c4 > 0, "C1_higher",
    fifelse(effect_c1_minus_c4 < 0, "C4_higher", "tie")
  )]
  write_tsv(msi_results, file.path(out_dir, "msi_binary_c1_vs_c4_by_cancer.tsv"))
  msi_cross <- summarize_cross_cancer(
    msi_results[formal_eligible == TRUE],
    "endpoint"
  )
  write_tsv(msi_cross, file.path(out_dir, "msi_cross_cancer_significant_counts.tsv"))
  log_msg("CHECKPOINT", "msi", "MSI binary tests complete")

  p_tmb <- make_continuous_forest(
    tmb_results,
    "TMB: C1 vs C4 within cancer",
    "Hodges-Lehmann shift in log10(TMB + 1), C1 - C4"
  )
  p_mantis <- make_continuous_forest(
    mantis_results,
    "MANTIS: C1 vs C4 within cancer",
    "Hodges-Lehmann shift in MANTIS score, C1 - C4"
  )
  save_figure_pair(
    p_tmb,
    file.path(fig_dir, "tmb_log10_c1_vs_c4_forest"),
    width = 6.8, height = 6.2
  )
  save_figure_pair(
    p_mantis,
    file.path(fig_dir, "mantis_continuous_c1_vs_c4_forest"),
    width = 6.8, height = 6.2
  )

  msi_plot_data <- msi_results[formal_eligible == TRUE]
  msi_plot_data[, Cancer := factor(
    Cancer,
    levels = rev(Cancer[order(effect_c1_minus_c4)])
  )]
  p_msi <- ggplot(msi_plot_data, aes(
    x = effect_c1_minus_c4,
    y = Cancer,
    colour = significant
  )) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.45) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D")) +
    labs(
      title = "MSI-H frequency: C1 vs C4 within cancer",
      x = "MSI-H frequency difference, C1 - C4",
      y = NULL,
      colour = "P < 0.05"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 11)
    )
  save_figure_pair(
    p_msi,
    file.path(fig_dir, "msi_binary_c1_vs_c4_forest"),
    width = 6.8, height = 5.8
  )

  qa_out <- data.table(
    check = c(
      "patient_unique", "only_c1_c4", "tmb_min_n_enforced",
      "mantis_min_n_enforced", "msi_min_n_enforced",
      "msi_rare_event_rule_enforced", "no_cross_cancer_bh"
    ),
    pass = c(
      !anyDuplicated(dat$Patient_ID),
      all(dat$C_subtype %chin% c("C1", "C4")),
      tmb_results[formal_eligible == TRUE, all(n_C1 >= 30L & n_C4 >= 30L)],
      mantis_results[formal_eligible == TRUE, all(n_C1 >= 30L & n_C4 >= 30L)],
      msi_results[formal_eligible == TRUE, all(n_C1 >= 30L & n_C4 >= 30L)],
      msi_results[formal_eligible == TRUE, all(MSI_H_C1 + MSI_H_C4 >= 5L)],
      !any(c("q_value", "padj") %chin% names(tmb_results))
    )
  )
  write_tsv(qa_out, file.path(out_dir, "tmb_msi_qa.tsv"))
  if (!all(qa_out$pass)) stop("TMB/MSI QA failed")
  save_session_info(file.path(out_dir, "sessionInfo.txt"))
  log_msg("DONE", "tmb_msi", "Formal scalar endpoint analyses and figures complete")
}, error = function(e) {
  log_msg("ERROR", "tmb_msi", conditionMessage(e))
  quit(status = 1L)
})
}

# Purpose: Compare copy-number burden across the reported state groups.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/genomics/03_run_cnv_burden_c1c4.R
stage_genomics_03_run_cnv_burden_c1c4 <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
readiness_dir <- unname(paths["readiness_dir"])
if(!dir.exists(readiness_dir) || file.access(readiness_dir,5L)!=0L) stop("Readiness input must be a readable searchable directory.")
required_names <- c("cnv_c1c4_sample_manifest.tsv", "readiness_qa.tsv")
required_files <- vapply(file.path(readiness_dir,required_names),resolve_path,character(1L),base_dir=dirname(config_path))

required_files <- c(required_files, unname(paths[c("continuous_cnv","thresholded_cnv")]))
if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
protected <- c(readiness_dir,dirname(required_files),config_path,package_root)
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

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# 独立输出根；保留完整主端点、敏感性与源诊断图。
output_root <- out_dir
out_dir <- file.path(output_root, "results")
fig_dir <- file.path(output_root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
continuous_path <- unname(paths["continuous_cnv"])
thresholded_path <- unname(paths["thresholded_cnv"])
`log_msg` <- bind_stage_helper(shared_log_msg_016, environment())
`assert_nonempty` <- bind_stage_helper(shared_assert_nonempty_021, environment())
`assert_cols` <- bind_stage_helper(shared_assert_cols_054, environment())
`write_tsv` <- bind_stage_helper(shared_write_tsv_049, environment())
`save_session_info` <- bind_stage_helper(shared_save_session_info_019, environment())
`wilcox_summary` <- bind_stage_helper(shared_wilcox_summary_056, environment())
`summarize_cross_cancer` <- bind_stage_helper(shared_summarize_cross_cancer_033, environment())
`save_figure_pair` <- bind_stage_helper(shared_save_figure_pair_028, environment())


tryCatch({
  log_msg("START", "cnv_burden", "Formal C1 vs C4 CNV burden")
  manifest_path <- file.path(readiness_dir, "cnv_c1c4_sample_manifest.tsv")
  qa_path <- file.path(readiness_dir, "readiness_qa.tsv")
  invisible(lapply(c(manifest_path, qa_path, continuous_path, thresholded_path), assert_nonempty))
  if (!all(fread(qa_path)$pass)) stop("Readiness QA did not pass")

  manifest <- fread(manifest_path)
  assert_cols(manifest, c("Patient_ID", "Sample15", "Cancer", "C_subtype"), "CNV manifest")
  if (anyDuplicated(manifest$Patient_ID)) stop("Duplicate patients in CNV manifest")
  if (any(!manifest$C_subtype %chin% c("C1", "C4"))) stop("Non-C1/C4 labels in CNV manifest")
  sample_cols <- manifest$Sample15

  log_msg("CHECKPOINT", "continuous_read", paste("Selected samples:", length(sample_cols)))
  cont <- fread(continuous_path, select = c("Sample", sample_cols))
  genes <- cont[[1L]]
  cont_mat <- as.matrix(cont[, -1L])
  storage.mode(cont_mat) <- "double"
  if (any(!is.finite(cont_mat), na.rm = TRUE)) stop("Non-finite continuous CNV values")
  burden <- data.table(
    Sample15 = colnames(cont_mat),
    mean_abs_continuous = colMeans(abs(cont_mat), na.rm = TRUE)
  )
  rm(cont, cont_mat)
  invisible(gc())

  log_msg("CHECKPOINT", "thresholded_read", "Reading thresholded selected-sample matrix")
  thr <- fread(thresholded_path, select = c("Sample", sample_cols))
  if (!identical(genes, thr[[1L]])) stop("Continuous and thresholded gene vectors differ")
  thr_mat <- as.matrix(thr[, -1L])
  storage.mode(thr_mat) <- "integer"
  observed_values <- sort(unique(as.vector(thr_mat)))
  if (any(!observed_values %in% -2:2)) {
    stop("Thresholded CNV values outside -2:-1:0:1:2")
  }
  sensitivity <- data.table(
    Sample15 = colnames(thr_mat),
    fraction_abs_ge1 = colMeans(abs(thr_mat) >= 1L, na.rm = TRUE),
    fraction_gain_ge1 = colMeans(thr_mat >= 1L, na.rm = TRUE),
    fraction_loss_le_minus1 = colMeans(thr_mat <= -1L, na.rm = TRUE),
    fraction_amp_plus2 = colMeans(thr_mat == 2L, na.rm = TRUE),
    fraction_deepdel_minus2 = colMeans(thr_mat == -2L, na.rm = TRUE)
  )
  rm(thr, thr_mat)
  invisible(gc())

  patient_burden <- merge(manifest, burden, by = "Sample15", all.x = TRUE)
  patient_burden <- merge(patient_burden, sensitivity, by = "Sample15", all.x = TRUE)
  metric_cols <- c(
    "mean_abs_continuous", "fraction_abs_ge1", "fraction_gain_ge1",
    "fraction_loss_le_minus1", "fraction_amp_plus2", "fraction_deepdel_minus2"
  )
  if (patient_burden[, anyNA(.SD), .SDcols = metric_cols]) {
    stop("Missing CNV burden metric after exact sample join")
  }
  write_tsv(patient_burden, file.path(out_dir, "cnv_burden_c1c4_patient_data.tsv"))
  log_msg("CHECKPOINT", "burden", "Patient-level burden metrics calculated")

  results <- rbindlist(lapply(metric_cols, function(metric) {
    patient_burden[, {
      x <- get(metric)[C_subtype == "C1"]
      y <- get(metric)[C_subtype == "C4"]
      eligible <- length(x) >= 30L && length(y) >= 30L
      stats <- if (eligible) wilcox_summary(x, y) else rep(NA_real_, 5L)
      names(stats) <- c(
        "p_value", "hl_shift_c1_minus_c4", "hl_ci_low",
        "hl_ci_high", "cliff_delta_c1_minus_c4"
      )
      c(list(
        n_C1 = length(x),
        n_C4 = length(y),
        median_C1 = median(x),
        q1_C1 = unname(quantile(x, 0.25)),
        q3_C1 = unname(quantile(x, 0.75)),
        median_C4 = median(y),
        q1_C4 = unname(quantile(y, 0.25)),
        q3_C4 = unname(quantile(y, 0.75)),
        formal_eligible = eligible
      ), as.list(stats))
    }, by = Cancer][, metric := metric]
  }), fill = TRUE)
  results[, `:=`(
    formal_primary = metric == "mean_abs_continuous",
    effect_c1_minus_c4 = hl_shift_c1_minus_c4,
    significant = formal_eligible & !is.na(p_value) & p_value < 0.05,
    direction = fifelse(
      hl_shift_c1_minus_c4 > 0, "C1_higher",
      fifelse(hl_shift_c1_minus_c4 < 0, "C4_higher", "tie")
    )
  )]
  setcolorder(results, c("metric", "formal_primary", "Cancer"))
  write_tsv(results, file.path(out_dir, "cnv_burden_c1_vs_c4_by_cancer.tsv"))

  primary <- results[formal_primary == TRUE & formal_eligible == TRUE]
  cross <- summarize_cross_cancer(primary, "metric")
  write_tsv(cross, file.path(out_dir, "cnv_burden_cross_cancer_significant_counts.tsv"))

  pdat <- copy(primary)
  pdat[, Cancer := factor(Cancer, levels = rev(Cancer[order(effect_c1_minus_c4)]))]
  p <- ggplot(pdat, aes(
    x = effect_c1_minus_c4,
    y = Cancer,
    xmin = hl_ci_low,
    xmax = hl_ci_high,
    colour = significant
  )) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.45) +
    geom_errorbarh(height = 0.18, linewidth = 0.55) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#4D4D4D")) +
    labs(
      title = "Gene-weighted CNV burden: C1 vs C4 within cancer",
      x = "Hodges-Lehmann shift in mean absolute continuous GISTIC score, C1 - C4",
      y = NULL,
      colour = "P < 0.05"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 11)
    )
  save_figure_pair(
    p,
    file.path(fig_dir, "cnv_mean_abs_continuous_c1_vs_c4_forest"),
    width = 7.2, height = 6.2
  )

  qa <- data.table(
    check = c(
      "patient_unique", "only_c1_c4", "all_metrics_complete",
      "threshold_value_domain", "min_group_n_enforced",
      "primary_endpoint_unique", "no_cross_cancer_bh"
    ),
    pass = c(
      !anyDuplicated(patient_burden$Patient_ID),
      all(patient_burden$C_subtype %chin% c("C1", "C4")),
      !patient_burden[, anyNA(.SD), .SDcols = metric_cols],
      all(observed_values %in% -2:2),
      results[formal_eligible == TRUE, all(n_C1 >= 30L & n_C4 >= 30L)],
      uniqueN(results[formal_primary == TRUE, metric]) == 1L,
      !"q_value" %chin% names(results)
    )
  )
  write_tsv(qa, file.path(out_dir, "cnv_burden_qa.tsv"))
  if (!all(qa$pass)) stop("CNV burden QA failed")
  save_session_info(file.path(out_dir, "sessionInfo.txt"))
  log_msg("DONE", "cnv_burden", "CNV burden analysis and figure complete")
}, error = function(e) {
  log_msg("ERROR", "cnv_burden", conditionMessage(e))
  quit(status = 1L)
})
}

# Purpose: Compare copy-number burden across the reported state groups.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/genomics/04_run_cnv_burden_four_state.R
stage_genomics_04_run_cnv_burden_four_state <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
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

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
# 四状态分支从标签与连续矩阵独立接入，不使用03的统计结果。
options(stringsAsFactors = FALSE, warn = 1)

script_name <- "04_run_cnv_burden_four_state.R"
start_time <- Sys.time()
log_event <- function(level, msg, ...) {
  extra <- list(...)
  extra_text <- if (length(extra) == 0L) "" else paste(
    vapply(names(extra), function(nm) {
      sprintf('%s="%s"', nm, gsub('"', "'", as.character(extra[[nm]]), fixed = TRUE))
    }, character(1)), collapse = " "
  )
  cat(sprintf(
    '[%s] %s script=%s msg="%s"%s%s\n',
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
    gsub('"', "'", msg, fixed = TRUE),
    if (nzchar(extra_text)) " " else "", extra_text
  ))
  flush.console()
}

sample_type_priority <- function(sample15) {
  priority_order <- c("01", "03", "06", "02", "05", "07", "08", "09")
  ans <- match(substr(sample15, 14L, 15L), priority_order)
  ans[is.na(ans)] <- length(priority_order) + 1L
  ans
}

main <- function() {
  canonical_input <- unname(paths["canonical_labels"])
  cnv_input <- unname(paths["continuous_cnv"])
  order_input <- unname(paths["cancer_order"])
  palette_input <- unname(paths["state_palette"])
  fig4_root <- out_dir
  source_dir <- file.path(fig4_root, "source_data")
  export_dir <- file.path(fig4_root, "panel_exports")
  manifest_source <- file.path(source_dir, "fig4_c_four_state_cnv_sample_manifest.tsv")
  patient_source <- file.path(source_dir, "fig4_c_four_state_cnv_patient_data.tsv")
  statistics_source <- file.path(source_dir, "fig4_c_four_state_cnv_omnibus_statistics.tsv")
  summary_source <- file.path(source_dir, "fig4_c_four_state_cnv_group_summary.tsv")
  figure_manifest <- file.path(source_dir, "fig4_c_four_state_cnv_figure_manifest.tsv")
  output_prefix <- file.path(export_dir, "fig4_c_four_state_cnv_burden_violin")

  group_order <- c("C1", "C2", "C3", "C4")
  min_group_n <- 20L
  width_mm <- 170
  height_mm <- 42

  required <- c(canonical_input, cnv_input, order_input, palette_input)
  if (any(!file.exists(required))) {
    stop("Missing input: ", paste(required[!file.exists(required)], collapse = "; "))
  }
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  log_event(
    "START", "Fig4 C four-state CNV extension begin",
    canonical = canonical_input, cnv = cnv_input, output = output_prefix,
    min_group_n = min_group_n
  )

  canonical <- fread(canonical_input)
  stopifnot(all(c("Sample", "Cancer", "C_subtype") %in% names(canonical)))
  canonical[, Sample15 := substr(as.character(Sample), 1L, 15L)]
  canonical[, Patient_ID := substr(Sample15, 1L, 12L)]
  if (any(!canonical$C_subtype %chin% group_order)) {
    stop("Canonical labels outside C1-C4")
  }
  audit <- canonical[, .(
    n_cancers = uniqueN(Cancer), n_subtypes = uniqueN(C_subtype)
  ), by = Patient_ID]
  discordant <- audit[n_cancers > 1L | n_subtypes > 1L, Patient_ID]

  matrix_samples <- names(fread(cnv_input, nrows = 0L))[-1L]
  candidates <- canonical[
    !Patient_ID %chin% discordant & Sample15 %chin% matrix_samples
  ]
  candidates[, sample_priority := sample_type_priority(Sample15)]
  setorder(candidates, Patient_ID, sample_priority, Sample15)
  selected <- candidates[, .SD[1L], by = Patient_ID]
  sample_manifest <- selected[, .(
    Patient_ID, Sample15, Cancer, C_subtype,
    sample_type = substr(Sample15, 14L, 15L),
    selection_rule = "available_exact_sample15_then_sample_type_priority_then_barcode"
  )]
  if (anyDuplicated(sample_manifest$Patient_ID)) stop("Duplicate patients in manifest")
  if (anyDuplicated(sample_manifest$Sample15)) stop("Duplicate samples in manifest")
  if (any(!sample_manifest$C_subtype %chin% group_order)) stop("Non-four-state manifest")
  fwrite(sample_manifest, manifest_source, sep = "\t", quote = FALSE, na = "NA")

  log_event(
    "CHECKPOINT", "Four-state matrix manifest frozen",
    n_patients = nrow(sample_manifest), n_cancers = uniqueN(sample_manifest$Cancer),
    n_discordant_excluded = length(discordant)
  )

  # 使用与 Stage 4 相同的连续 GISTIC 主端点定义，仅扩展到完整四状态患者。
  sample_cols <- sample_manifest$Sample15
  cont <- fread(cnv_input, select = c("Sample", sample_cols), showProgress = TRUE)
  cont_mat <- as.matrix(cont[, -1L])
  storage.mode(cont_mat) <- "double"
  if (any(!is.finite(cont_mat), na.rm = TRUE)) {
    stop("Non-finite continuous CNV values")
  }
  burden <- data.table(
    Sample15 = colnames(cont_mat),
    mean_abs_continuous = colMeans(abs(cont_mat), na.rm = TRUE)
  )
  rm(cont, cont_mat)
  invisible(gc())
  patient <- merge(sample_manifest, burden, by = "Sample15", all.x = TRUE)
  if (patient[, anyNA(mean_abs_continuous)]) stop("Missing burden after join")

  order_df <- fread(order_input)
  cancer_order_all <- order_df$Cancer[order_df$Cancer %chin% patient$Cancer]
  patient[, display_order := match(Cancer, cancer_order_all)]
  setorder(patient, display_order, C_subtype, Patient_ID)
  fwrite(patient, patient_source, sep = "\t", quote = FALSE, na = "NA")

  counts <- dcast(
    patient[, .N, by = .(Cancer, C_subtype)],
    Cancer ~ C_subtype, value.var = "N", fill = 0L
  )
  for (g in group_order) if (!g %chin% names(counts)) counts[, (g) := 0L]
  counts[, formal_eligible := C1 >= min_group_n & C2 >= min_group_n &
    C3 >= min_group_n & C4 >= min_group_n]

  group_summary <- patient[, .(
    n = .N,
    median = median(mean_abs_continuous),
    q1 = unname(quantile(mean_abs_continuous, 0.25)),
    q3 = unname(quantile(mean_abs_continuous, 0.75)),
    observed_max = max(mean_abs_continuous)
  ), by = .(Cancer, C_subtype)]

  eligible_cancers <- counts[formal_eligible == TRUE, Cancer]
  omnibus <- rbindlist(lapply(eligible_cancers, function(cancer) {
    dat <- patient[Cancer == cancer]
    test <- kruskal.test(mean_abs_continuous ~ C_subtype, data = dat)
    data.table(
      Cancer = cancer,
      kruskal_chisq = unname(test$statistic),
      df = unname(test$parameter),
      p_value = unname(test$p.value)
    )
  }))
  statistics <- merge(counts, omnibus, by = "Cancer", all.x = TRUE, sort = FALSE)
  statistics[, significant := formal_eligible & !is.na(p_value) & p_value < 0.05]
  statistics[, significance_code := fifelse(
    is.na(p_value) | p_value >= 0.05, "",
    fifelse(p_value < 0.0001, "****",
      fifelse(p_value < 0.001, "***", fifelse(p_value < 0.01, "**", "*")))
  )]
  statistics[, display_order := match(Cancer, cancer_order_all)]
  setorder(statistics, display_order)
  fwrite(statistics, statistics_source, sep = "\t", quote = FALSE, na = "NA")
  fwrite(group_summary, summary_source, sep = "\t", quote = FALSE, na = "NA")

  selected_cancers <- statistics[significant == TRUE, Cancer]
  if (length(selected_cancers) == 0L) {
    stop("No eligible cancer passed four-state Kruskal-Wallis P < 0.05")
  }
  selected_cancers <- cancer_order_all[cancer_order_all %chin% selected_cancers]
  plot_dat <- patient[Cancer %chin% selected_cancers]
  plot_dat[, Cancer := factor(Cancer, levels = selected_cancers)]
  plot_dat[, C_subtype := factor(C_subtype, levels = group_order)]
  marker <- group_summary[Cancer %chin% selected_cancers, .(
    marker_y = max(observed_max) + 0.018
  ), by = Cancer]
  marker <- merge(
    marker,
    statistics[, .(Cancer, significance_code, p_value)],
    by = "Cancer", all.x = TRUE, sort = FALSE
  )
  marker[, Cancer := factor(Cancer, levels = selected_cancers)]

  palette <- fread(palette_input)
  group_colors <- setNames(
    palette$color[match(group_order, palette$state)], group_order
  )
  # 为最高显著性标记保留最终尺寸下的顶部净空，避免星号贴边裁切。
  y_max <- max(marker$marker_y) + 0.040
  dodge <- position_dodge(width = 0.82)
  p <- ggplot(
    plot_dat,
    aes(
      x = Cancer, y = mean_abs_continuous, fill = C_subtype,
      group = interaction(Cancer, C_subtype)
    )
  ) +
    geom_violin(
      scale = "width", trim = TRUE, alpha = 0.84,
      position = dodge, width = 0.88,
      linewidth = 0.22, colour = "black"
    ) +
    geom_boxplot(
      width = 0.11, outlier.shape = NA, alpha = 0.48,
      position = dodge, linewidth = 0.22, colour = "black"
    ) +
    geom_text(
      data = marker,
      aes(x = Cancer, y = marker_y, label = significance_code),
      inherit.aes = FALSE, family = "Arial", size = 2.35,
      colour = "black", vjust = 0.30
    ) +
    scale_fill_manual(values = group_colors, breaks = group_order, name = NULL) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = "Mean absolute continuous CNV value") +
    guides(fill = guide_legend(direction = "horizontal", nrow = 1, byrow = TRUE)) +
    theme_classic(base_size = 6, base_family = "Arial") +
    theme(
      panel.grid.major.y = element_line(colour = "#E2E2E2", linewidth = 0.08),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.13),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.10),
      axis.ticks.length = grid::unit(1.0, "mm"),
      axis.text.x = element_text(size = 5.4, colour = "black", margin = margin(t = 1.0)),
      axis.text.y = element_text(size = 5.2, colour = "black"),
      axis.title.y = element_text(size = 5.5, colour = "black", margin = margin(r = 1.5)),
      legend.position = "bottom",
      legend.location = "plot",
      legend.justification = "center",
      legend.text = element_text(size = 5.0, colour = "black"),
      legend.key.size = grid::unit(2.5, "mm"),
      legend.spacing.x = grid::unit(0.8, "mm"),
      legend.margin = margin(t = 0.5, r = 0, b = 0, l = 0),
      plot.margin = margin(t = 1.5, r = 1.5, b = 0.8, l = 1.5)
    )

  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(
    paste0(output_prefix, ".pdf"), p,
    width = width_in, height = height_in, units = "in", device = cairo_pdf
  )
  ggsave(
    paste0(output_prefix, ".tiff"), p,
    width = width_in, height = height_in, units = "in",
    dpi = 600, device = "tiff", compression = "lzw"
  )
  ggsave(
    paste0(output_prefix, ".png"), p,
    width = width_in, height = height_in, units = "in", dpi = 300, bg = "white"
  )

  manifest <- data.table(
    field = c(
      "panel", "status", "analysis_class", "canonical_input", "cnv_input",
      "cancer_order_input", "palette_input", "sample_manifest", "patient_source",
      "statistics_source", "group_summary_source", "output_prefix",
      "width_mm", "height_mm", "endpoint", "group_order", "min_group_n",
      "formal_test", "selection_threshold", "multiple_testing",
      "n_selected_patients", "n_eligible_cancers", "n_significant_cancers",
      "significant_cancers_in_fig1_order", "evidence_boundary"
    ),
    value = c(
      "C", "public rerun candidate (unfrozen)", "Fig4-layer four-state extension",
      canonical_input, cnv_input, order_input, palette_input, manifest_source,
      patient_source, statistics_source, summary_source, output_prefix,
      width_mm, height_mm, "mean_abs_continuous", paste(group_order, collapse = ";"),
      min_group_n, "Cancer-wise Kruskal-Wallis omnibus test", "raw P < 0.05",
      "No secondary cross-cancer BH", nrow(sample_manifest),
      length(eligible_cancers), length(selected_cancers),
      paste(selected_cancers, collapse = ";"),
      "Four-state background association only; direct formal mainline comparisons remain frozen cancer-wise C1 versus C4"
    )
  )
  fwrite(manifest, figure_manifest, sep = "\t", quote = FALSE, na = "NA")

  expected <- c(
    manifest_source, patient_source, statistics_source, summary_source,
    figure_manifest, paste0(output_prefix, c(".pdf", ".tiff", ".png"))
  )
  missing <- expected[!file.exists(expected) | file.info(expected)$size <= 0]
  if (length(missing) > 0L) stop("Missing outputs: ", paste(missing, collapse = "; "))

  log_event(
    "CHECKPOINT", "Four-state omnibus and figure completed",
    n_patients = nrow(sample_manifest), n_eligible_cancers = length(eligible_cancers),
    n_significant_cancers = length(selected_cancers),
    significant_cancers = paste(selected_cancers, collapse = ";")
  )
  log_event(
    "DONE", "Fig4 C four-state CNV extension completed",
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2),
    output = output_prefix
  )
}

tryCatch(
  main(),
  error = function(e) {
    log_event("ERROR", conditionMessage(e))
    quit(status = 1)
  }
)
}

# Purpose: Aggregate significant gene-level copy-number differences into descriptive genomic regions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/genomics/05_run_cnv_alteration_regions.R
stage_genomics_05_run_cnv_alteration_regions <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
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

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})
# 四项输入与独立输出根参数化；科学和图形主体保持来源实现。
manifest_file <- unname(paths["manifest"])
cnv_file <- unname(paths["thresholded_cnv"])
probe_file <- unname(paths["gene_coordinates"])
cytoband_file <- unname(paths["cytoband"])
source_dir <- file.path(out_dir, "source_data")
export_dir <- file.path(out_dir, "panel_exports")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

MIN_N <- 30L
ALT_FREQ_MIN <- 0.01
FDR_MAX <- 0.05
DIFF_MIN <- 0.15
MAX_GAP <- 10e6
TOP_N <- 15L
STAR_GENES <- c("EGFR", "BRAF", "TP53", "KRAS", "MYC", "PIK3CA", "PTEN",
                "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")

message("[1/7] Read C1/C4 manifest and selected CNV matrix for an unfrozen public rerun")
manifest <- fread(manifest_file)
stopifnot(!anyDuplicated(manifest$Patient_ID), all(manifest$C_subtype %chin% c("C1", "C4")))
coverage <- dcast(manifest[, .N, by = .(Cancer, C_subtype)], Cancer ~ C_subtype,
                  value.var = "N", fill = 0L)
coverage[, formal_eligible := C1 >= MIN_N & C4 >= MIN_N]
eligible_cancers <- coverage[formal_eligible == TRUE, Cancer]
if (!length(eligible_cancers)) stop("No eligible cancers")

thr <- fread(cnv_file, select = c("Sample", manifest$Sample15), showProgress = TRUE)
genes <- as.character(thr[[1L]])
if (anyDuplicated(genes)) stop("Duplicate genes in CNV matrix")
mat <- as.matrix(thr[, -1L])
storage.mode(mat) <- "integer"
if (anyNA(mat) || any(!unique(as.vector(mat)) %in% -2:2)) stop("Invalid thresholded CNV values")
rm(thr)
invisible(gc())
manifest[, matrix_col := match(Sample15, colnames(mat))]
if (anyNA(manifest$matrix_col)) stop("Manifest samples missing from matrix")

# 与旧版相同的低变异预筛思想；在当前冻结C1/C4样本宇宙中计算。
alt_freq <- rowMeans(mat != 0)
keep_gene <- alt_freq > ALT_FREQ_MIN
mat <- mat[keep_gene, , drop = FALSE]
genes <- genes[keep_gene]
message("Genes after alt_freq > 0.01: ", length(genes))

fisher_vec <- function(a, b, n1, n4) {
  vapply(seq_along(a), function(i) {
    if (a[i] == 0L && b[i] == 0L) return(1)
    fisher.test(matrix(c(a[i], n1 - a[i], b[i], n4 - b[i]), nrow = 2L))$p.value
  }, numeric(1L))
}

message("[2/7] Run within-cancer C1 vs C4 Fisher tests")
all_stats <- vector("list", length(eligible_cancers))
for (ii in seq_along(eligible_cancers)) {
  cancer <- eligible_cancers[ii]
  idx1 <- manifest[Cancer == cancer & C_subtype == "C1", matrix_col]
  idx4 <- manifest[Cancer == cancer & C_subtype == "C4", matrix_col]
  n1 <- length(idx1); n4 <- length(idx4)
  c1_amp <- rowSums(mat[, idx1, drop = FALSE] > 0)
  c4_amp <- rowSums(mat[, idx4, drop = FALSE] > 0)
  c1_del <- rowSums(mat[, idx1, drop = FALSE] < 0)
  c4_del <- rowSums(mat[, idx4, drop = FALSE] < 0)
  amp <- data.table(Gene = genes, Cancer = cancer, Type = "Amp", n_C1 = n1, n_C4 = n4,
                    event_C1 = c1_amp, event_C4 = c4_amp,
                    Freq_C1 = c1_amp / n1, Freq_C4 = c4_amp / n4)
  del <- data.table(Gene = genes, Cancer = cancer, Type = "Del", n_C1 = n1, n_C4 = n4,
                    event_C1 = c1_del, event_C4 = c4_del,
                    Freq_C1 = c1_del / n1, Freq_C4 = c4_del / n4)
  amp[, `:=`(Diff = Freq_C1 - Freq_C4, P = fisher_vec(event_C1, event_C4, n1, n4))]
  del[, `:=`(Diff = Freq_C1 - Freq_C4, P = fisher_vec(event_C1, event_C4, n1, n4))]
  amp[, FDR := p.adjust(P, "BH")]
  del[, FDR := p.adjust(P, "BH")]
  z <- rbind(amp, del)
  z[, `:=`(Abs_Diff = abs(Diff), Bias = fifelse(Diff > 0, "C1 higher", "C4 higher"))]
  z[, significant := FDR < FDR_MAX & Abs_Diff > DIFF_MIN]
  all_stats[[ii]] <- z
  message(sprintf("  %s (%d/%d): Amp=%d, Del=%d", cancer, ii, length(eligible_cancers),
                  amp[FDR < FDR_MAX & abs(Diff) > DIFF_MIN, .N],
                  del[FDR < FDR_MAX & abs(Diff) > DIFF_MIN, .N]))
}
stats <- rbindlist(all_stats)
fwrite(stats, file.path(source_dir, "fig4_d_cnv_legacy_rule_all_gene_statistics.tsv.gz"),
       sep = "\t", compress = "gzip")
sig <- stats[significant == TRUE]
if (!nrow(sig)) stop("No significant events under legacy-adapted rules")

message("[3/7] Map significant genes to GENCODE v23 coordinates")
probe <- fread(probe_file)
probe <- probe[chrom %chin% paste0("chr", c(1:22, "X", "Y"))]
probe[, symbol_key := sub("\\|.*$", "", gene)]
coords <- probe[, .(Chr = unique(chrom)[1L], Start = min(chromStart), End = max(chromEnd),
                    n_chrom = uniqueN(chrom)), by = symbol_key]
coords <- coords[n_chrom == 1L]
setnames(coords, "symbol_key", "Gene")
mapped <- merge(sig, coords[, .(Gene, Chr, Start, End)], by = "Gene", all = FALSE)
mapping_summary <- data.table(metric = c("significant_gene_events", "mapped_gene_events", "mapping_fraction"),
                              value = c(nrow(sig), nrow(mapped), nrow(mapped) / nrow(sig)))
if (nrow(mapped) / nrow(sig) < 0.90) stop("Coordinate mapping below 90%")

message("[4/7] Merge significant genes into <=10 Mb alteration regions")
setorder(mapped, Cancer, Type, Bias, Chr, Start)
mapped[, new_cluster := as.integer(seq_len(.N) == 1L | Start - shift(End, fill = -1e9) > MAX_GAP),
       by = .(Cancer, Type, Bias, Chr)]
mapped[, cluster_id := cumsum(new_cluster), by = .(Cancer, Type, Bias, Chr)]
regions <- mapped[, {
  stars <- intersect(Gene, STAR_GENES)
  representative <- if (length(stars)) stars[1L] else Gene[which.max(Abs_Diff)]
  .(Start = min(Start), End = max(End), Mean_Diff = mean(Abs_Diff), Max_Diff = max(Abs_Diff),
    Gene_Count = .N, Representative_Gene = representative,
    Star_Genes = paste(stars, collapse = ";"), Is_Star = length(stars) > 0,
    Min_FDR = min(FDR))
}, by = .(Cancer, Type, Bias, Chr, cluster_id)]

message("[5/7] Apply old Top-15 plus star-region display rule")
display_regions <- regions[, {
  stars <- .SD[Is_Star == TRUE]
  ordinary <- .SD[Is_Star == FALSE][order(-Mean_Diff)][seq_len(min(.N, TOP_N))]
  unique(rbind(stars, ordinary, fill = TRUE))
}, by = Cancer]
setorder(display_regions, Cancer, Chr, Start, Type, Bias)

message("[6/7] Build two-column hg38 karyotype")
left_chroms <- paste0("chr", c(1:11, "X"))
right_chroms <- paste0("chr", c(22:12, "Y"))
chrom_levels <- c(left_chroms, right_chroms)
cyto <- fread(cytoband_file,
              col.names = c("chrom", "chromStart", "chromEnd", "band", "gieStain"))
cyto <- cyto[chrom %chin% chrom_levels]
sizes <- cyto[, .(chromLength = max(chromEnd)), by = chrom]
layout <- rbind(
  data.table(chrom = left_chroms, panel_side = "left", row_index = rev(seq_along(left_chroms))),
  data.table(chrom = right_chroms, panel_side = "right", row_index = rev(seq_along(right_chroms)))
)
layout <- merge(layout, sizes, by = "chrom", all.x = TRUE, sort = FALSE)
layout[, row_y := row_index * 0.52]
left_max <- max(layout[panel_side == "left", chromLength])
right_max <- max(layout[panel_side == "right", chromLength])
gap <- 70e6
right_end <- left_max + gap + right_max
to_x <- function(pos, side, len) ifelse(side == "left", pos, right_end - len + pos)
cyto <- merge(cyto, layout, by = "chrom", all.x = TRUE, sort = FALSE)
cyto[, `:=`(x_start = to_x(chromStart, panel_side, chromLength),
            x_end = to_x(chromEnd, panel_side, chromLength))]
outline <- copy(layout)
outline[, `:=`(x_start = fifelse(panel_side == "left", 0, right_end - chromLength),
               x_end = fifelse(panel_side == "left", chromLength, right_end))]
display_regions <- merge(display_regions, layout, by.x = "Chr", by.y = "chrom", all.x = TRUE)
display_regions[, `:=`(mid_pos = (Start + End) / 2,
                       x = to_x((Start + End) / 2, panel_side, chromLength),
                       y = row_y + fifelse(Bias == "C1 higher", 0.115, -0.115))]

labels <- display_regions[Is_Star == TRUE & nzchar(Star_Genes)]
labels[, Label := sub(";.*$", "", Star_Genes)]
setorder(labels, Label, -Mean_Diff)
labels <- labels[, .SD[1L], by = Label]

stain_cols <- c(gneg = "#FFFFFF", gpos25 = "#D9D9D9", gpos50 = "#A6A6A6",
                gpos75 = "#737373", gpos100 = "#111111", acen = "#111111",
                gvar = "#D9D9D9", stalk = "#D9D9D9")
cyto[!gieStain %chin% names(stain_cols), gieStain := "gneg"]
cancers <- sort(unique(display_regions$Cancer))
shape_pool <- c(0:14, 15:18)
if (length(cancers) > length(shape_pool)) stop("Too many cancers for stable shape map")
shape_values <- setNames(shape_pool[seq_along(cancers)], cancers)
axis_breaks <- sort(unique(layout$row_y))
left_axis_labels <- layout[panel_side == "left"][order(row_y), sub("^chr", "", chrom)]
right_axis_labels <- layout[panel_side == "right"][order(row_y), sub("^chr", "", chrom)]

p <- ggplot() +
  geom_segment(data = outline, aes(x = x_start, xend = x_end, y = row_y, yend = row_y),
               linewidth = 2.0, lineend = "round", color = "black") +
  geom_segment(data = cyto, aes(x = x_start, xend = x_end, y = row_y, yend = row_y, color = gieStain),
               linewidth = 1.45, lineend = "butt", show.legend = FALSE) +
  geom_point(data = display_regions, aes(x = x, y = y, shape = Cancer, color = Bias),
             size = 2.15, stroke = 0.62, fill = NA, alpha = 0.94) +
  geom_text_repel(data = labels, aes(x = x, y = y, label = Label),
                  size = 2.15, color = "black", box.padding = 0.24, point.padding = 0.12,
                  min.segment.length = 0, segment.size = 0.18, max.overlaps = Inf,
                  seed = 20260809) +
  scale_color_manual(values = c(stain_cols, "C1 higher" = "#FF2D2D", "C4 higher" = "#1764FF"),
                     breaks = c("C1 higher", "C4 higher"), name = "Enrichment bias") +
  scale_shape_manual(values = shape_values, name = NULL) +
  scale_y_continuous(breaks = axis_breaks, labels = left_axis_labels,
                     sec.axis = dup_axis(breaks = axis_breaks, labels = right_axis_labels, name = NULL),
                     expand = expansion(add = c(0.34, 0.10))) +
  scale_x_continuous(expand = expansion(add = c(8e6, 8e6))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Chromosome") +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.title.y = element_text(size = 7.0, margin = margin(r = 4)),
        axis.text.y = element_text(size = 6.0, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.3), axis.ticks.length.y = unit(1.4, "mm"),
        axis.line.y = element_line(linewidth = 0.3),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.position = "right", legend.box = "vertical", legend.text = element_text(size = 5.8),
        legend.title = element_text(size = 6.0), legend.key.height = unit(7.2, "pt"),
        plot.margin = margin(2.0, 2.5, 2.0, 5.0, unit = "mm")) +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.2)),
         shape = guide_legend(order = 2, override.aes = list(color = "black", size = 2.0)))

message("[7/7] Write unfrozen public rerun source data, QA and exports")
fwrite(display_regions, file.path(source_dir, "fig4_d_cnv_chromosome_map.tsv"), sep = "\t")
fwrite(labels, file.path(source_dir, "fig4_d_cnv_chromosome_labels.tsv"), sep = "\t")
fwrite(layout, file.path(source_dir, "fig4_d_cnv_chromosome_layout.tsv"), sep = "\t")
fwrite(coverage, file.path(source_dir, "fig4_d_cnv_legacy_rule_coverage.tsv"), sep = "\t")
fwrite(mapping_summary, file.path(source_dir, "fig4_d_cnv_legacy_rule_mapping_summary.tsv"), sep = "\t")
qa <- data.table(
  check = c("patient_unique", "only_c1_c4", "min_group_n", "within_cancer_bh",
            "formal_threshold", "mapping_ge90pct", "no_cross_cancer_pooling"),
  pass = c(!anyDuplicated(manifest$Patient_ID), all(manifest$C_subtype %chin% c("C1", "C4")),
           all(stats$n_C1 >= 30 & stats$n_C4 >= 30), TRUE,
           all(sig$FDR < 0.05 & sig$Abs_Diff > 0.15), nrow(mapped) / nrow(sig) >= 0.90, TRUE)
)
fwrite(qa, file.path(source_dir, "fig4_d_cnv_legacy_rule_qa.tsv"), sep = "\t")
if (!all(qa$pass)) stop("QA failed")
manifest_out <- data.table(
  field = c("panel", "method", "comparison", "eligible_cancers", "significant_gene_events",
            "mapped_gene_events", "display_regions", "display_cancers", "label_genes",
            "width_mm", "height_mm", "evidence_boundary"),
  value = c("D", "legacy CNV landscape rules adapted to canonical C1 vs C4", "within-cancer C1 vs C4",
            length(eligible_cancers), nrow(sig), nrow(mapped), nrow(display_regions), length(cancers),
            paste(labels$Label, collapse = ";"), 180, 72,
            "10 Mb alteration clusters; not GISTIC focal peaks or arm-level events")
)
fwrite(manifest_out, file.path(source_dir, "fig4_d_cnv_chromosome_map_manifest.tsv"), sep = "\t")

stem <- file.path(export_dir, "fig4_d_cnv_chromosome_map")
ggsave(paste0(stem, ".pdf"), p, width = 180, height = 72, units = "mm", device = cairo_pdf)
ggsave(paste0(stem, ".tiff"), p, width = 180, height = 72, units = "mm", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(paste0(stem, ".png"), p, width = 180, height = 72, units = "mm", dpi = 300,
       device = "png", bg = "white")
writeLines(capture.output(sessionInfo()), file.path(source_dir, "fig4_d_cnv_legacy_rule_sessionInfo.txt"))
message(sprintf("DONE: %d significant gene-events, %d display regions, %d cancers, %d labels",
                nrow(sig), nrow(display_regions), length(cancers), nrow(labels)))
}
