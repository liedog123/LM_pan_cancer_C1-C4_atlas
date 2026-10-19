#!/usr/bin/env Rscript
# Purpose: Compare mutation burden and MSI measures using the adopted eligibility rules.
# Task: tmb_msi; step: s02
# Inputs (DCF fields): readiness_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task tmb_msi --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 02_run_tmb_msi_c1c4.R --config PATH [--check-config]\n       Rscript 02_run_tmb_msi_c1c4.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("readiness_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
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
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

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

log_msg <- function(level, step, ...) {
  detail <- paste(..., collapse = " ")
  cat(sprintf(
    "%s\t%s\t%s\t%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level, step, detail
  ))
  flush.console()
}


assert_nonempty <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }
  size <- file.info(path)$size
  if (is.na(size) || size <= 0) {
    stop("Empty input: ", path)
  }
  invisible(size)
}


assert_cols <- function(x, required, object_name) {
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) {
    stop(object_name, " missing columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}


write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}


save_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), path)
}


wilcox_summary <- function(x, y) {
  test <- suppressWarnings(wilcox.test(
    x, y,
    alternative = "two.sided",
    exact = FALSE,
    conf.int = TRUE,
    conf.level = 0.95
  ))
  combined <- c(x, y)
  ranks <- rank(combined, ties.method = "average")
  n_x <- length(x)
  n_y <- length(y)
  u_x <- sum(ranks[seq_len(n_x)]) - n_x * (n_x + 1) / 2
  cliff_delta <- 2 * u_x / (n_x * n_y) - 1

  c(
    p_value = unname(test$p.value),
    hl_shift_c1_minus_c4 = if (length(test$estimate)) unname(test$estimate) else NA_real_,
    hl_ci_low = if (!is.null(test$conf.int)) unname(test$conf.int[1L]) else NA_real_,
    hl_ci_high = if (!is.null(test$conf.int)) unname(test$conf.int[2L]) else NA_real_,
    cliff_delta_c1_minus_c4 = cliff_delta
  )
}


summarize_cross_cancer <- function(
  results,
  feature_cols,
  sig_col = "significant",
  effect_col = "effect_c1_minus_c4"
) {
  x <- copy(results)
  x[, direction := fifelse(
    get(effect_col) > 0, "C1_higher",
    fifelse(get(effect_col) < 0, "C4_higher", "tie")
  )]
  x[, sig_direction := fifelse(
    get(sig_col) & direction == "C1_higher", "C1_higher",
    fifelse(get(sig_col) & direction == "C4_higher", "C4_higher", "not_significant")
  )]

  out <- x[, .(
    n_testable_cancers = uniqueN(Cancer),
    n_sig_total = sum(get(sig_col), na.rm = TRUE),
    n_sig_C1_higher = sum(sig_direction == "C1_higher", na.rm = TRUE),
    n_sig_C4_higher = sum(sig_direction == "C4_higher", na.rm = TRUE),
    median_effect_c1_minus_c4 = median(get(effect_col), na.rm = TRUE),
    sig_cancers_C1_higher = paste(sort(Cancer[sig_direction == "C1_higher"]), collapse = ";"),
    sig_cancers_C4_higher = paste(sort(Cancer[sig_direction == "C4_higher"]), collapse = ";")
  ), by = feature_cols]
  out[, dominant_n_sig := pmax(n_sig_C1_higher, n_sig_C4_higher)]
  out[, dominant_direction := fifelse(
    n_sig_C1_higher > n_sig_C4_higher, "C1_higher",
    fifelse(n_sig_C4_higher > n_sig_C1_higher, "C4_higher", "tie")
  )]
  out[, context_dependent := n_sig_C1_higher > 0L & n_sig_C4_higher > 0L]
  out[, recurrent_ge3 := dominant_n_sig >= 3L]
  out
}


save_figure_pair <- function(plot, stem, width, height, dpi = 600L) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot = plot,
    width = width, height = height, units = "in",
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"), plot = plot,
    width = width, height = height, units = "in",
    dpi = dpi, compression = "lzw"
  )
  ggplot2::ggsave(
    paste0(stem, "_preview.png"), plot = plot,
    width = width, height = height, units = "in",
    dpi = 180
  )
}


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
