#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 04_run_cnv_burden_four_state.R --config PATH [--check-config]\n       Rscript 04_run_cnv_burden_four_state.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("canonical_labels", "continuous_cnv", "cancer_order", "state_palette", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
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
