#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 03_run_cnv_burden_c1c4.R --config PATH [--check-config]\n       Rscript 03_run_cnv_burden_c1c4.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("readiness_dir", "continuous_cnv", "thresholded_cnv", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
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
# 独立输出根；保留完整主端点、敏感性与源诊断图。
output_root <- out_dir
out_dir <- file.path(output_root, "results")
fig_dir <- file.path(output_root, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
continuous_path <- unname(paths["continuous_cnv"])
thresholded_path <- unname(paths["thresholded_cnv"])

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
