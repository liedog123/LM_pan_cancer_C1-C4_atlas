#!/usr/bin/env Rscript
# Purpose: Compare immune and stromal features with the specified state contrasts and testing families.
# Task: immune_ecology; step: s02
# Inputs (DCF fields): input_results
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task immune_ecology --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 02_finalize_stage2_evidence.R --config PATH [--check-config]\n",
      "       Rscript 02_finalize_stage2_evidence.R --help\n\n",
      "DCF fields: input_results, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary invocation classifies evidence and produces source diagnostics and plots.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (!(length(args) %in% c(2L, 3L)) || args[1L] != "--config" || !nzchar(args[2L]) ||
    (length(args) == 3L && args[3L] != "--check-config")) {
  usage()
  stop("Expected --config PATH, optionally followed by --check-config.", call. = FALSE)
}
flags <- if (length(args) == 3L) args[3L] else character()

# Resolve DCF-relative paths through existing ancestors and symlinks.
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
if (!file_test("-f", config_path) || file.access(config_path, mode = 4L) != 0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size <= 0) {
  stop("Configuration must be a readable nonempty regular file.", call. = FALSE)
}
cfg <- read.dcf(config_path)
required_fields <- c("input_results", "output_dir")
if (nrow(cfg) != 1L || !setequal(colnames(cfg), required_fields) ||
    ncol(cfg) != length(required_fields)) {
  stop("Configuration must contain one DCF record with exactly: ",
       paste(required_fields, collapse = ", "), call. = FALSE)
}
cfg <- setNames(as.character(cfg[1L, required_fields]), required_fields)
if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME", cfg, fixed = TRUE))) {
  stop("Replace all example placeholders with explicit paths.", call. = FALSE)
}
paths <- vapply(cfg, resolve_path, character(1L), base_dir = dirname(config_path))
input_result_dir <- unname(paths["input_results"])
out_dir <- unname(paths["output_dir"])
if (!dir.exists(input_result_dir) || file.access(input_result_dir, mode = 5L) != 0L) {
  stop("input_results must be a readable searchable directory.", call. = FALSE)
}
input_stems <- c("immune_random_effects_meta",
                 "immune_leave_one_cancer_out",
                 "immune_pairwise_effects_by_cancer",
                 "immune_continuous_axis_associations_by_cancer",
                 "immune_omnibus_by_cancer",
                 "immune_adjusted_models_by_cancer",
                 "immune_input_and_feature_registry",
                 "cibersort_confidence_retention_audit",
                 "immune_state_profile_equal_weight",
                 "stage2_sample_analysis_manifest",
                 "immune_coverage_and_evaluability_by_cancer")
required_files <- file.path(input_result_dir, paste0(input_stems, "_v0.1_20260723.tsv"))
# Resolve each input's final target before checking output independence.
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)

# Protect inputs and existing results; checks create no directories.
protected_paths <- c(dirname(required_files), input_result_dir, config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
if ("--check-config" %in% flags) {
  cat("Configuration paths passed metadata checks. No inputs loaded; no analysis or outputs created.\n")
  quit(status = 0L)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
script_name <- "02_finalize_stage2_evidence.R"
version_tag <- "v0.1_20260723"
smoke_mode <- FALSE
# 输入统计表只读；新分类、诊断及绘图写入独立输出目录。
result_dir <- file.path(out_dir, "results")
figure_dir <- file.path(out_dir, "figures")
doc_dir <- file.path(out_dir, "docs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)

log_event <- function(level, msg, ...) {
  extras <- list(...)
  extra_text <- if (length(extras) == 0L) "" else paste(
    vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)),
    collapse = " "
  )
  cat(sprintf("[%s] %s script=%s msg=%s%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg),
              ifelse(nzchar(extra_text), " ", ""), extra_text))
  flush.console()
}

fail <- function(msg) {
  log_event("ERROR", msg)
  stop(msg, call. = FALSE)
}

read_result <- function(stem) {
  path <- file.path(input_result_dir, paste0(stem, "_", version_tag, ".tsv"))
  if (!file.exists(path) || file.info(path)$size <= 0) fail(paste("missing result", path))
  fread(path, na.strings = c("NA", "NaN", ""))
}

write_result <- function(x, stem) {
  path <- file.path(result_dir, paste0(stem, "_", version_tag, ".tsv"))
  fwrite(x, path, sep = "\t", na = "NA", quote = TRUE)
  path
}

md5_safe <- function(path) unname(tools::md5sum(path))
mm_to_in <- function(x) x / 25.4

result_theme <- function(base_size = 6.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.3),
      strip.background = element_blank(), strip.text = element_text(size = base_size),
      panel.grid = element_blank(), plot.title = element_blank()
    )
}

save_figure <- function(p, stem, width_mm, height_mm) {
  prefix <- file.path(figure_dir, paste0(stem, "_", version_tag))
  ggsave(paste0(prefix, ".pdf"), p, width = mm_to_in(width_mm), height = mm_to_in(height_mm),
         units = "in", device = cairo_pdf)
  ggsave(paste0(prefix, ".tiff"), p, width = mm_to_in(width_mm), height = mm_to_in(height_mm),
         units = "in", dpi = 600, compression = "lzw")
  c(pdf = paste0(prefix, ".pdf"), tiff = paste0(prefix, ".tiff"))
}

add_sensitivity <- function(base, sens, prefix) {
  keys <- c("Method", "Feature", "Broad_Compartment", "Design_ID", "Design_Type")
  keep <- sens[, c(keys, "Beta_REML", "Meta_FDR", "I2", "LOO_Direction_Stable"), with = FALSE]
  setnames(keep, c("Beta_REML", "Meta_FDR", "I2", "LOO_Direction_Stable"),
           paste0(prefix, c("_Beta", "_FDR", "_I2", "_LOO_Direction_Stable")))
  merge(base, keep, by = keys, all.x = TRUE, sort = FALSE)
}

direction_status <- function(pos, neg) {
  if (pos > 0L && neg == 0L) return("POSITIVE")
  if (neg > 0L && pos == 0L) return("NEGATIVE")
  if (pos > 0L && neg > 0L) return("MIXED")
  "NONE"
}

main <- function() {
  start_time <- Sys.time()
  log_event("START", "Stage 2 evidence finalization begin", input = input_result_dir,
            output = paste(figure_dir, doc_dir, sep = ";"), parameters = ifelse(smoke_mode, "SMOKE", "FORMAL"))

  meta <- read_result("immune_random_effects_meta")
  loo <- read_result("immune_leave_one_cancer_out")
  pairwise <- read_result("immune_pairwise_effects_by_cancer")
  continuous <- read_result("immune_continuous_axis_associations_by_cancer")
  omnibus <- read_result("immune_omnibus_by_cancer")
  adjusted <- read_result("immune_adjusted_models_by_cancer")
  registry <- read_result("immune_input_and_feature_registry")
  retention <- read_result("cibersort_confidence_retention_audit")
  state_profile <- read_result("immune_state_profile_equal_weight")
  sample_manifest <- read_result("stage2_sample_analysis_manifest")
  coverage <- read_result("immune_coverage_and_evaluability_by_cancer")

  main_set <- "REPRESENTATIVE_MAIN__ALL_MATCHED"
  conf_set <- "REPRESENTATIVE_MAIN__CONFIDENT_P_LT_0.05"
  evidence <- copy(meta[Model_Type == "UNADJUSTED" & Analysis_Set == main_set])
  if (nrow(evidence) == 0L) fail("no primary meta results")

  age_sens <- meta[Model_Type == "CORE_AGE_SEX" & Analysis_Set == main_set]
  purity_sens <- meta[Model_Type == "CORE_PURITY" & Analysis_Set == main_set]
  conf_sens <- meta[Method == "CIBERSORT_ABS" & Model_Type == "UNADJUSTED" & Analysis_Set == conf_set]
  evidence <- add_sensitivity(evidence, age_sens, "AgeSex")
  evidence <- add_sensitivity(evidence, purity_sens, "Purity")
  evidence <- add_sensitivity(evidence, conf_sens, "CIBConfidence")

  evidence[, AgeSex_Direction_Compatible := fifelse(is.na(AgeSex_Beta), NA, sign(Beta_REML) == sign(AgeSex_Beta))]
  evidence[, Purity_Direction_Compatible := fifelse(is.na(Purity_Beta), NA, sign(Beta_REML) == sign(Purity_Beta))]
  evidence[, CIBConfidence_Direction_Compatible := fifelse(is.na(CIBConfidence_Beta), NA,
                                                            sign(Beta_REML) == sign(CIBConfidence_Beta))]
  evidence[, Base_Cross_Cancer_Support := Meta_FDR < 0.05 & Direction_Consistency >= 0.70 & LOO_Direction_Stable]
  evidence[, Sensitivity_Discordant :=
             (!is.na(AgeSex_Direction_Compatible) & !AgeSex_Direction_Compatible) |
             (!is.na(Purity_Direction_Compatible) & !Purity_Direction_Compatible) |
             (!is.na(CIBConfidence_Direction_Compatible) & !CIBConfidence_Direction_Compatible)]
  evidence[, Evidence_Class := fifelse(
    Base_Cross_Cancer_Support & (I2 > 75 | Sensitivity_Discordant), "CONTEXT_DEPENDENT",
    fifelse(Base_Cross_Cancer_Support, "CROSS_CANCER_SUPPORTED",
            fifelse(P_Value < 0.05, "DESCRIPTIVE_NOMINAL", "NO_STABLE_SIGNAL"))
  )]
  evidence[, Direction := fifelse(Beta_REML > 0, "POSITIVE", fifelse(Beta_REML < 0, "NEGATIVE", "ZERO"))]

  method_direction <- evidence[Evidence_Class == "CROSS_CANCER_SUPPORTED",
    .(N_Positive = sum(Beta_REML > 0), N_Negative = sum(Beta_REML < 0),
      Supported_Features = paste(sort(Feature), collapse = ";")),
    by = .(Broad_Compartment, Design_ID, Method)]
  method_direction[, Method_Direction := mapply(direction_status, N_Positive, N_Negative)]

  all_compartments <- unique(evidence[, .(Broad_Compartment, Design_ID)])
  all_methods <- CJ(Method = c("ESTIMATE", "CIBERSORT_ABS", "MCPCOUNTER"), unique = TRUE)
  # 使用显式临时键构造笛卡尔积，避免 data.table 版本差异导致 by = NULL 失败。
  all_compartments[, Join_Key := 1L]
  all_methods[, Join_Key := 1L]
  convergence <- merge(all_compartments, all_methods, by = "Join_Key", allow.cartesian = TRUE)
  convergence[, Join_Key := NULL]
  convergence <- merge(convergence, method_direction,
                       by = c("Broad_Compartment", "Design_ID", "Method"), all.x = TRUE)
  convergence[is.na(Method_Direction), `:=`(Method_Direction = "NONE", N_Positive = 0L,
                                             N_Negative = 0L, Supported_Features = "")]
  convergence_summary <- convergence[, .(
    Positive_Methods = sum(Method_Direction == "POSITIVE"),
    Negative_Methods = sum(Method_Direction == "NEGATIVE"),
    Mixed_Methods = sum(Method_Direction == "MIXED"),
    Positive_Method_List = paste(Method[Method_Direction == "POSITIVE"], collapse = ";"),
    Negative_Method_List = paste(Method[Method_Direction == "NEGATIVE"], collapse = ";")),
    by = .(Broad_Compartment, Design_ID)]
  convergence_summary[, Convergent_Direction := fifelse(
    Positive_Methods >= 2L & Negative_Methods == 0L, "POSITIVE",
    fifelse(Negative_Methods >= 2L & Positive_Methods == 0L, "NEGATIVE",
            fifelse(Positive_Methods > 0L & Negative_Methods > 0L, "CONFLICT", "NONE")))]
  convergence_summary[, Method_Convergent := Convergent_Direction %in% c("POSITIVE", "NEGATIVE")]
  convergence <- merge(convergence, convergence_summary,
                       by = c("Broad_Compartment", "Design_ID"), all.x = TRUE)

  conv_key <- convergence_summary[Method_Convergent == TRUE,
                                  .(Broad_Compartment, Design_ID, Convergent_Direction)]
  evidence <- merge(evidence, conv_key, by = c("Broad_Compartment", "Design_ID"), all.x = TRUE)
  evidence[Evidence_Class == "CROSS_CANCER_SUPPORTED" &
             ((Beta_REML > 0 & Convergent_Direction == "POSITIVE") |
              (Beta_REML < 0 & Convergent_Direction == "NEGATIVE")),
           Evidence_Class := "METHOD_CONVERGENT"]

  cancer_specific <- rbindlist(list(
    pairwise[Analysis_Set == main_set & Status == "ESTIMATED" & FDR < 0.05,
             .(Cancer, Method, Feature, Broad_Compartment, Design_ID, Design_Type,
               Effect = Hedges_g, Effect_Metric = "HEDGES_G", P_Value, FDR)],
    continuous[Analysis_Set == main_set & Status == "ESTIMATED" & FDR < 0.05,
               .(Cancer, Method, Feature, Broad_Compartment, Design_ID, Design_Type,
                 Effect = Spearman_rho, Effect_Metric = "SPEARMAN_RHO", P_Value, FDR)]
  ), fill = TRUE)
  cancer_specific[, Evidence_Class := "CANCER_SPECIFIC_SUPPORTED"]

  gate_checks <- data.table(
    Check = c("THREE_METHOD_INPUTS", "CANONICAL_PATIENT_INTERFACE", "CANCER_LEVEL_ANALYSIS_COMPLETE",
              "META_AND_LOO_COMPLETE", "DEFINITION_CIRCULARITY_COMPLETE", "METHOD_CONVERGENT_CORE",
              "CIBERSORT_DUAL_ANALYSIS_COMPLETE", "POSITIVE_NEGATIVE_CONTEXT_REPORTED"),
    Passed = c(
      all(c("ESTIMATE", "CIBERSORT_ABS", "MCPCOUNTER") %in% unique(evidence$Method)),
      nrow(sample_manifest) > 0L && !anyDuplicated(sample_manifest$Patient_ID),
      all(c("ESTIMATE", "CIBERSORT_ABS", "MCPCOUNTER") %in% unique(coverage$Method)),
      nrow(meta) > 0L && nrow(loo) > 0L,
      nrow(registry) == 37L && all(registry$Circularity_Class == "FIXED_SIGNATURE_OVERLAP_AUDITED_NO_RESCORING"),
      any(convergence_summary$Method_Convergent),
      any(meta$Analysis_Set == conf_set),
      all(c("CONTEXT_DEPENDENT", "NO_STABLE_SIGNAL") %in% unique(evidence$Evidence_Class)) ||
        any(evidence$Evidence_Class %in% c("METHOD_CONVERGENT", "CROSS_CANCER_SUPPORTED"))
    ),
    Detail = c(
      paste(sort(unique(evidence$Method)), collapse = ";"),
      sprintf("n_patients=%d", nrow(sample_manifest)),
      sprintf("n_coverage_rows=%d", nrow(coverage)),
      sprintf("n_meta=%d;n_loo=%d", nrow(meta), nrow(loo)),
      sprintf("n_signature_rows=%d", nrow(registry)),
      sprintf("n_convergent_compartment_design=%d", sum(convergence_summary$Method_Convergent)),
      sprintf("n_confident_meta=%d", nrow(meta[Analysis_Set == conf_set])),
      paste(sort(unique(evidence$Evidence_Class)), collapse = ";"))
  )
  gate_decision <- if (all(gate_checks$Passed)) {
    "PASS_STAGE2_COMPLETE_GO_STAGE3"
  } else if (all(gate_checks$Passed[gate_checks$Check != "METHOD_CONVERGENT_CORE"])) {
    "HOLD_STAGE2_COMPLETE_NO_METHOD_CONVERGENT_CORE"
  } else {
    "HOLD_STAGE2_CONTRACT_INCOMPLETE"
  }
  gate <- rbind(gate_checks, data.table(Check = "GATE2_DECISION",
                                        Passed = startsWith(gate_decision, "PASS"), Detail = gate_decision))

  evidence_path <- write_result(evidence, "immune_evidence_classification")
  cancer_path <- write_result(cancer_specific, "cancer_specific_immune_evidence_classification")
  convergence_path <- write_result(convergence, "immune_cross_method_convergence_matrix")
  gate_path <- write_result(gate, "stage2_gate")

  figure_rows <- list()
  design_order <- c("AG1_IG1_vs_AG2_IG1", "AG2_IG1_vs_AG2_IG2",
                    "AG2_IG2_vs_AG3_IG2", "C1_vs_C4", "CONTINUOUS_ACID", "CONTINUOUS_IG")
  heat_palette <- c("-1" = "#2166AC", "0" = "#F7F7F7", "1" = "#B2182B")
  for (method in c("ESTIMATE", "CIBERSORT_ABS", "MCPCOUNTER")) {
    src <- evidence[Method == method]
    if (nrow(src) == 0L) next
    src[, Design_ID := factor(Design_ID, levels = design_order)]
    src[, Feature := factor(Feature, levels = rev(unique(Feature)))]
    src[, Effect_Capped := pmax(-1, pmin(1, Beta_REML))]
    p <- ggplot(src, aes(x = Design_ID, y = Feature, fill = Effect_Capped)) +
      geom_tile(colour = "white", linewidth = 0.15) +
      geom_point(data = src[Meta_FDR < 0.05], shape = 21, size = 1.15,
                 stroke = 0.25, fill = "black", colour = "black") +
      scale_fill_gradient2(low = heat_palette[["-1"]], mid = heat_palette[["0"]],
                           high = heat_palette[["1"]], midpoint = 0, limits = c(-1, 1),
                           name = "Meta effect") +
      labs(x = NULL, y = NULL) +
      result_theme(5.8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
    size <- if (method == "CIBERSORT_ABS") c(120, 120) else if (method == "MCPCOUNTER") c(120, 75) else c(90, 60)
    stem <- paste0("immune_cross_cancer_meta_heatmap_", tolower(method))
    files <- save_figure(p, stem, size[1], size[2])
    source_path <- file.path(figure_dir, paste0(stem, "_source_data_", version_tag, ".tsv"))
    fwrite(src, source_path, sep = "\t", na = "NA", quote = TRUE)
    figure_rows[[length(figure_rows) + 1L]] <- data.table(
      Figure_ID = stem, Evidence_Role = "immune_ecology_association", Output_Mode = "independent_single_plot",
      Figure_Type = "cross_cancer_meta_heatmap", Width_mm = size[1], Height_mm = size[2],
      PDF = files[["pdf"]], TIFF = files[["tiff"]], Source_Data = source_path,
      Statistics = "REML meta; BH within method-design-analysis-model",
      Evidence_Boundary = "association_not_mechanism")
  }

  conv_src <- convergence_summary[Method_Convergent == TRUE]
  if (nrow(conv_src) > 0L) {
    conv_src[, Design_ID := factor(Design_ID, levels = design_order)]
    conv_src[, Broad_Compartment := factor(Broad_Compartment,
                                            levels = rev(sort(unique(Broad_Compartment))))]
    p_conv <- ggplot(conv_src, aes(x = Design_ID, y = Broad_Compartment,
                                   colour = Convergent_Direction,
                                   size = pmax(Positive_Methods, Negative_Methods))) +
      geom_point(alpha = 0.9) +
      scale_colour_manual(values = c(POSITIVE = "#B2182B", NEGATIVE = "#2166AC"), name = "Direction") +
      scale_size_continuous(range = c(2, 4), breaks = c(2, 3), name = "Methods") +
      labs(x = NULL, y = NULL) + result_theme(6) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
    stem <- "immune_cross_method_convergence_dotplot"
    files <- save_figure(p_conv, stem, 120, 75)
    source_path <- file.path(figure_dir, paste0(stem, "_source_data_", version_tag, ".tsv"))
    fwrite(conv_src, source_path, sep = "\t", na = "NA", quote = TRUE)
    figure_rows[[length(figure_rows) + 1L]] <- data.table(
      Figure_ID = stem, Evidence_Role = "orthogonal_method_convergence", Output_Mode = "independent_single_plot",
      Figure_Type = "dotplot", Width_mm = 120, Height_mm = 75,
      PDF = files[["pdf"]], TIFF = files[["tiff"]], Source_Data = source_path,
      Statistics = "two-or-more-method directional convergence",
      Evidence_Boundary = "method_convergence_not_causality")
  }

  retention_src <- retention[Patient_Set == "REPRESENTATIVE_MAIN"]
  retention_src[, C_subtype := factor(C_subtype, levels = c("C1", "C2", "C3", "C4"))]
  state_colors <- c(C1 = "#D55E00", C2 = "#E69F00", C3 = "#0072B2", C4 = "#7A7A7A")
  p_ret <- ggplot(retention_src, aes(x = C_subtype, y = Retention_Rate, colour = C_subtype)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.65) +
    scale_colour_manual(values = state_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "CIBERSORT P < 0.05 retention") + result_theme(6.2)
  stem <- "cibersort_confidence_retention_by_state"
  files <- save_figure(p_ret, stem, 90, 60)
  source_path <- file.path(figure_dir, paste0(stem, "_source_data_", version_tag, ".tsv"))
  fwrite(retention_src, source_path, sep = "\t", na = "NA", quote = TRUE)
  figure_rows[[length(figure_rows) + 1L]] <- data.table(
    Figure_ID = stem, Evidence_Role = "QC_background", Output_Mode = "independent_single_plot",
    Figure_Type = "distribution", Width_mm = 90, Height_mm = 60,
    PDF = files[["pdf"]], TIFF = files[["tiff"]], Source_Data = source_path,
    Statistics = "state-specific retention; per-cancer dependency audit",
    Evidence_Boundary = "algorithm_confidence_QC_not_biology")

  figure_manifest <- rbindlist(figure_rows, fill = TRUE)
  figure_manifest_path <- write_result(figure_manifest, "figure_manifest")

  class_counts <- evidence[, .N, by = Evidence_Class][order(Evidence_Class)]
  top_supported <- evidence[Evidence_Class %in% c("METHOD_CONVERGENT", "CROSS_CANCER_SUPPORTED")][
    order(Meta_FDR, -abs(Beta_REML))]
  top_supported <- head(top_supported, 25)
  top_lines <- if (nrow(top_supported) == 0L) {
    "- 无跨癌稳定特征。"
  } else {
    paste0("- ", top_supported$Method, " / ", top_supported$Feature, " / ", top_supported$Design_ID,
           ": effect=", sprintf("%.3f", top_supported$Effect_REML_Interpretable),
           ", FDR=", format(top_supported$Meta_FDR, digits = 3),
           ", I2=", sprintf("%.1f%%", top_supported$I2),
           ", class=", top_supported$Evidence_Class)
  }
  conv_lines <- if (nrow(conv_src) == 0L) {
    "- 无满足两种方法同向的compartment × design单元。"
  } else {
    paste0("- ", conv_src$Broad_Compartment, " / ", conv_src$Design_ID,
           ": ", conv_src$Convergent_Direction,
           "，methods=", pmax(conv_src$Positive_Methods, conv_src$Negative_Methods))
  }
  report_path <- file.path(doc_dir, paste0("stage2_immune_ecology_report_", version_tag, ".md"))
  report_lines <- c(
    "# Stage 2 Immune Ecology Report V0.1", "", "日期：2026-07-23", "",
    paste0("状态：`STAGE2_EXECUTION_COMPLETE / ", gate_decision, "`"), "",
    "## 1. 执行范围", "",
    sprintf("本Stage复用ESTIMATE、CIBERSORT absolute和MCPcounter样本级分数，正式患者接口包含%d名患者、%d个癌种。所有检验先在癌种内完成；未重新运行去卷积，未读取旧Subtype作为标签。",
            nrow(sample_manifest), uniqueN(sample_manifest$Cancer)), "",
    "固定比较包括C1-C4总体、三条原子局部边、C1 vs C4和continuous Acid/IG；CIBERSORT同时完成ALL_MATCHED主分析和P<0.05置信样本敏感性。", "",
    "## 2. 分析规模与证据分类", "",
    sprintf("癌种级结果包括%d条omnibus、%d条pairwise、%d条continuous和%d条调整模型记录；跨癌产生%d条REML meta和%d条LOO记录。",
            nrow(omnibus), nrow(pairwise), nrow(continuous), nrow(adjusted), nrow(meta), nrow(loo)), "",
    paste0("证据分类：", paste(class_counts$Evidence_Class, class_counts$N, sep = "=", collapse = "；"), "。"), "",
    "## 3. 跨癌稳定免疫生态特征", "", top_lines, "",
    "上述结果只支持bulk免疫丰度/生态关联；固定signature与Acid/LM12定义基因的重叠已审计，不将其提升为因果机制。", "",
    "## 4. 多方法一致性", "", conv_lines, "",
    "只有至少两种互补方法同向的单元标记为METHOD_CONVERGENT；同一方法内多个细胞亚群不计作多方法互证。", "",
    "## 5. CIBERSORT与纯度敏感性", "",
    sprintf("CIBERSORT置信样本meta记录%d条；代表样本中状态依赖保留率审计覆盖%d个癌种。任何ALL_MATCHED与P<0.05方向不兼容的结果均降级。",
            nrow(meta[Analysis_Set == conf_set]), uniqueN(retention_src$Cancer)), "",
    "CIBERSORT和MCPcounter均完成年龄+性别及可用TumorPurity调整；ESTIMATE不进行同源派生指标间的循环调整。", "",
    "## 6. Gate 2", "",
    paste0("Gate 2决定：`", gate_decision, "`。"), "",
    paste0("- ", gate_checks$Check, ": ", ifelse(gate_checks$Passed, "PASS", "FAIL"), "（", gate_checks$Detail, "）"), "",
    "Gate PASS只表示免疫生态层按合同完成并形成可审计的正向、阴性和异质结果，不等于建立免疫机制或临床标志物。", "",
    "## 7. 主要输出", "",
    paste0("- 证据分类：`../results/", basename(evidence_path), "`；"),
    paste0("- 多方法矩阵：`../results/", basename(convergence_path), "`；"),
    paste0("- Gate 2：`../results/", basename(gate_path), "`；"),
    paste0("- 图件manifest：`../results/", basename(figure_manifest_path), "`。")
  )
  writeLines(report_lines, report_path, useBytes = TRUE)

  all_outputs <- c(evidence_path, cancer_path, convergence_path, gate_path,
                   figure_manifest_path, report_path,
                   unlist(figure_manifest[, .(PDF, TIFF, Source_Data)], use.names = FALSE))
  all_outputs <- unique(all_outputs[file.exists(all_outputs)])
  output_manifest <- data.table(
    File = all_outputs, Exists = file.exists(all_outputs), Size_Bytes = file.info(all_outputs)$size,
    MD5 = vapply(all_outputs, md5_safe, character(1)), Stage = "STAGE2_FINALIZATION",
    Smoke_Mode = smoke_mode)
  manifest_path <- write_result(output_manifest, "stage2_output_manifest")

  log_event("DONE", "Stage 2 evidence finalization complete",
            runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
            gate = gate_decision, n_supported = nrow(top_supported),
            n_convergent = sum(convergence_summary$Method_Convergent),
            report = report_path, manifest = manifest_path)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1)
})
