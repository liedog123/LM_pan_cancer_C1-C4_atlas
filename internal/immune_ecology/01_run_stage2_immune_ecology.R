#!/usr/bin/env Rscript
# Purpose: Compare immune and stromal features with the specified state contrasts and testing families.
# Task: immune_ecology; step: s01
# Inputs (DCF fields): canonical_labels,patient_manifest,estimate_scores,cibersort_scores,mcpcounter_scores,acid_genes,lm12_genes
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task immune_ecology --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 01_run_stage2_immune_ecology.R --config PATH [--check-config]\n",
      "       Rscript 01_run_stage2_immune_ecology.R --help\n\n",
      "DCF fields: canonical_labels, patient_manifest, estimate_scores, cibersort_scores, mcpcounter_scores, acid_genes, lm12_genes, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary invocation runs the full source statistical family; no deconvolution is rerun.\n",
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
required_fields <- c("canonical_labels", "patient_manifest", "estimate_scores", "cibersort_scores", "mcpcounter_scores", "acid_genes", "lm12_genes", "output_dir")
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
canonical_file <- unname(paths["canonical_labels"])
patient_file <- unname(paths["patient_manifest"])
estimate_file <- unname(paths["estimate_scores"])
cibersort_file <- unname(paths["cibersort_scores"])
mcpcounter_file <- unname(paths["mcpcounter_scores"])
acid_gene_file <- unname(paths["acid_genes"])
lm12_gene_file <- unname(paths["lm12_genes"])
out_dir <- unname(paths["output_dir"])
required_files <- c(canonical_file, patient_file, estimate_file, cibersort_file,
                    mcpcounter_file, acid_gene_file, lm12_gene_file)
# Resolve each input's final target before checking output independence.
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)

# Protect inputs and existing results; checks create no directories.
protected_paths <- c(dirname(required_files), config_path)
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
  library(IOBR)
  library(sandwich)
  library(lmtest)
})
script_name <- "01_run_stage2_immune_ecology.R"
version_tag <- "v0.1_20260723"
# 正式入口固定完整分析范围；不提供按癌种截取的运行模式。
smoke_mode <- FALSE
selected_cancers <- character()
output_dir <- file.path(out_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

standardize_safe <- function(x) {
  x <- safe_numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

rint_safe <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n < 3L) return(out)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- qnorm((r - 3 / 8) / (n + 1 / 4))
  out
}

md5_safe <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

design_specs <- list(
  C1C4_OMNIBUS = list(type = "omnibus", variable = "C_subtype",
                      levels = c("C4", "C1", "C2", "C3"), target = NA_character_, reference = "C4"),
  AG1_IG1_vs_AG2_IG1 = list(type = "pair", variable = "AGIG",
                            levels = c("AG2_IG1", "AG1_IG1"), target = "AG1_IG1", reference = "AG2_IG1"),
  AG2_IG1_vs_AG2_IG2 = list(type = "pair", variable = "AGIG",
                            levels = c("AG2_IG2", "AG2_IG1"), target = "AG2_IG1", reference = "AG2_IG2"),
  AG2_IG2_vs_AG3_IG2 = list(type = "pair", variable = "AGIG",
                            levels = c("AG3_IG2", "AG2_IG2"), target = "AG2_IG2", reference = "AG3_IG2"),
  C1_vs_C4 = list(type = "pair", variable = "C_subtype",
                  levels = c("C4", "C1"), target = "C1", reference = "C4"),
  CONTINUOUS_ACID = list(type = "continuous", variable = "Acid_z",
                         levels = NULL, target = "per_1SD", reference = NA_character_),
  CONTINUOUS_IG = list(type = "continuous", variable = "IG_z",
                       levels = NULL, target = "per_1SD", reference = NA_character_)
)

estimate_features <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
cibersort_features <- c(
  "B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8",
  "T_cells_CD4_naive", "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated",
  "T_cells_follicular_helper", "T_cells_regulatory_(Tregs)", "T_cells_gamma_delta",
  "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0",
  "Macrophages_M1", "Macrophages_M2", "Dendritic_cells_resting",
  "Dendritic_cells_activated", "Mast_cells_resting", "Mast_cells_activated",
  "Eosinophils", "Neutrophils", "Absolute_score_(sig_score)_CIBERSORT"
)
mcpcounter_features <- c(
  "T_cells", "CD8_T_cells", "Cytotoxic_lymphocytes", "B_lineage", "NK_cells",
  "Monocytic_lineage", "Myeloid_dendritic_cells", "Neutrophils",
  "Endothelial_cells", "Fibroblasts"
)

broad_compartment <- function(method, feature) {
  if (method == "ESTIMATE") {
    return(c(StromalScore = "STROMAL", ImmuneScore = "GLOBAL_IMMUNE",
             ESTIMATEScore = "GLOBAL_TME", TumorPurity = "PURITY")[[feature]])
  }
  if (method == "MCPCOUNTER") {
    return(c(T_cells = "T_CELL", CD8_T_cells = "CD8_T",
             Cytotoxic_lymphocytes = "CYTOTOXIC_LYMPHOCYTE", B_lineage = "B_LINEAGE",
             NK_cells = "NK", Monocytic_lineage = "MONOCYTIC",
             Myeloid_dendritic_cells = "DENDRITIC", Neutrophils = "NEUTROPHILS",
             Endothelial_cells = "ENDOTHELIAL", Fibroblasts = "FIBROBLAST")[[feature]])
  }
  if (feature == "Absolute_score_(sig_score)_CIBERSORT") return("GLOBAL_IMMUNE")
  if (grepl("^B_cells|^Plasma", feature)) return("B_LINEAGE")
  if (feature == "T_cells_CD8") return("CD8_T")
  if (grepl("^T_cells_", feature)) return("T_CELL")
  if (grepl("^NK_cells", feature)) return("NK")
  if (feature == "Monocytes") return("MONOCYTIC")
  if (grepl("^Macrophages", feature)) return("MACROPHAGE")
  if (grepl("^Dendritic", feature)) return("DENDRITIC")
  if (grepl("^Mast", feature)) return("MAST")
  if (feature == "Eosinophils") return("EOSINOPHILS")
  if (feature == "Neutrophils") return("NEUTROPHILS")
  "OTHER"
}

make_signature_registry <- function() {
  acid_env <- new.env(parent = emptyenv())
  loaded <- load(acid_gene_file, envir = acid_env)
  if (!("final_selected_genes" %in% loaded)) fail("final_selected_genes not found")
  acid_genes <- unique(toupper(as.character(acid_env$final_selected_genes$Gene)))
  lm12_dt <- fread(lm12_gene_file)
  gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol"), names(lm12_dt))[1]
  if (is.na(gene_col)) fail("LM12 gene column not found")
  lm12_genes <- unique(toupper(as.character(lm12_dt[[gene_col]])))

  data("lm22", package = "IOBR", envir = environment())
  data("mcp_genes", package = "IOBR", envir = environment())
  data("SI_geneset", package = "IOBR", envir = environment())

  signature_lists <- list()
  lm22_names <- gsub(" ", "_", colnames(lm22), fixed = TRUE)
  for (i in seq_along(lm22_names)) {
    signature_lists[[paste("CIBERSORT_ABS", lm22_names[i], sep = "::")]] <-
      toupper(rownames(lm22)[is.finite(lm22[[i]]) & lm22[[i]] != 0])
  }
  signature_lists[["CIBERSORT_ABS::Absolute_score_(sig_score)_CIBERSORT"]] <- toupper(rownames(lm22))

  mcp_map <- c(
    "T cells" = "T_cells", "CD8 T cells" = "CD8_T_cells",
    "Cytotoxic lymphocytes" = "Cytotoxic_lymphocytes", "B lineage" = "B_lineage",
    "NK cells" = "NK_cells", "Monocytic lineage" = "Monocytic_lineage",
    "Myeloid dendritic cells" = "Myeloid_dendritic_cells", "Neutrophils" = "Neutrophils",
    "Endothelial cells" = "Endothelial_cells", "Fibroblasts" = "Fibroblasts"
  )
  for (pop in names(mcp_map)) {
    genes <- toupper(as.character(mcp_genes[mcp_genes[["Cell population"]] == pop, "HUGO symbols"]))
    signature_lists[[paste("MCPCOUNTER", mcp_map[[pop]], sep = "::")]] <- genes
  }

  si_stromal <- toupper(as.character(unlist(SI_geneset["StromalSignature", ], use.names = FALSE)))
  si_immune <- toupper(as.character(unlist(SI_geneset["ImmuneSignature", ], use.names = FALSE)))
  si_stromal <- setdiff(si_stromal, c("", "NA", "ESTIMATE"))
  si_immune <- setdiff(si_immune, c("", "NA", "ESTIMATE"))
  signature_lists[["ESTIMATE::StromalScore"]] <- si_stromal
  signature_lists[["ESTIMATE::ImmuneScore"]] <- si_immune
  signature_lists[["ESTIMATE::ESTIMATEScore"]] <- union(si_stromal, si_immune)
  signature_lists[["ESTIMATE::TumorPurity"]] <- union(si_stromal, si_immune)

  rows <- lapply(names(signature_lists), function(key) {
    bits <- strsplit(key, "::", fixed = TRUE)[[1]]
    genes <- unique(signature_lists[[key]][nzchar(signature_lists[[key]])])
    acid_overlap <- intersect(genes, acid_genes)
    lm12_overlap <- intersect(genes, lm12_genes)
    data.table(
      Method = bits[1], Feature = bits[2], Broad_Compartment = broad_compartment(bits[1], bits[2]),
      Signature_Source = ifelse(bits[1] == "CIBERSORT_ABS", "IOBR::lm22",
                                ifelse(bits[1] == "MCPCOUNTER", "IOBR::mcp_genes", "IOBR::SI_geneset")),
      Signature_Version = as.character(packageVersion("IOBR")), Signature_Gene_Count = length(genes),
      Acid_Definition_Gene_Count = length(acid_genes), LM12_Definition_Gene_Count = length(lm12_genes),
      Acid_Overlap_N = length(acid_overlap), Acid_Overlap_Proportion = length(acid_overlap) / max(1, length(genes)),
      LM12_Overlap_N = length(lm12_overlap), LM12_Overlap_Proportion = length(lm12_overlap) / max(1, length(genes)),
      Acid_Overlap_Genes = paste(sort(acid_overlap), collapse = ";"),
      LM12_Overlap_Genes = paste(sort(lm12_overlap), collapse = ";"),
      Fixed_Public_Signature = TRUE, Non_Definition_Rescore = FALSE,
      Circularity_Class = "FIXED_SIGNATURE_OVERLAP_AUDITED_NO_RESCORING"
    )
  })
  rbindlist(rows, fill = TRUE)
}

read_method_scores <- function(method, path, features, patient_dt) {
  if (!file.exists(path)) fail(paste("missing input", path))
  x <- fread(path, na.strings = c("NA", "NaN", ""))
  if ("Subtype" %in% names(x)) x[, Subtype := NULL]
  required <- c("Sample", "Cancer", features)
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) fail(paste(method, "missing columns", paste(missing_cols, collapse = ",")))
  x[, Sample := substr(as.character(Sample), 1, 15)]
  if (anyDuplicated(x$Sample)) fail(paste(method, "score table is not unique by Sample"))
  setnames(x, "Cancer", "Cancer_Source")
  keep <- c("Sample", "Cancer_Source", features)
  if (method == "CIBERSORT_ABS") keep <- c(keep, "P-value", "Correlation", "RMSE")
  x <- x[, ..keep]
  out <- merge(patient_dt, x, by = "Sample", all = FALSE, sort = FALSE)
  if (any(out$Cancer != out$Cancer_Source)) fail(paste(method, "Cancer mismatch after join"))
  out[, Cancer_Source := NULL]
  for (f in features) out[, (f) := safe_numeric(get(f))]
  if (method == "CIBERSORT_ABS") {
    out[, `P-value` := safe_numeric(`P-value`)]
    out[, Correlation := safe_numeric(Correlation)]
    out[, RMSE := safe_numeric(RMSE)]
  }
  out
}

prepare_design <- function(dt, spec, feature) {
  x <- copy(dt[is.finite(get(feature))])
  x[, y := get(feature)]
  x[, y_rint := rint_safe(y)]
  if (spec$type == "continuous") {
    x[, exposure := safe_numeric(get(spec$variable))]
    x <- x[is.finite(exposure)]
  } else {
    x <- x[get(spec$variable) %in% spec$levels]
    x[, exposure := factor(get(spec$variable), levels = spec$levels)]
  }
  x
}

check_threshold <- function(x, spec) {
  if (spec$type == "omnibus") {
    tab <- table(x$exposure)
    ok <- length(tab) == 4L && all(tab >= 10L) && nrow(x) >= 60L
    detail <- paste0("group_n=", paste(names(tab), tab, sep = ":", collapse = "|"), ";total_n=", nrow(x))
  } else if (spec$type == "pair") {
    tab <- table(x$exposure)
    ok <- length(tab) == 2L && all(tab >= 10L) && nrow(x) >= 30L
    detail <- paste0("group_n=", paste(names(tab), tab, sep = ":", collapse = "|"), ";total_n=", nrow(x))
  } else {
    ok <- nrow(x) >= 50L && is.finite(sd(x$exposure)) && sd(x$exposure) > 0 &&
      is.finite(sd(x$y)) && sd(x$y) > 0
    detail <- sprintf("n=%d;exposure_sd=%.5f;feature_sd=%.5f", nrow(x), sd(x$exposure), sd(x$y))
  }
  list(ok = isTRUE(ok), detail = detail)
}

hedges_g <- function(target, reference) {
  target <- target[is.finite(target)]
  reference <- reference[is.finite(reference)]
  n1 <- length(target); n0 <- length(reference)
  if (n1 < 2L || n0 < 2L) return(c(g = NA_real_, se = NA_real_))
  sp2 <- ((n1 - 1) * var(target) + (n0 - 1) * var(reference)) / (n1 + n0 - 2)
  if (!is.finite(sp2) || sp2 <= 0) return(c(g = NA_real_, se = NA_real_))
  d <- (mean(target) - mean(reference)) / sqrt(sp2)
  j <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- j * d
  se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2)))
  c(g = g, se = se)
}

rank_biserial <- function(target, reference) {
  n1 <- length(target); n0 <- length(reference)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(c(target, reference), ties.method = "average")
  u <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u / (n1 * n0) - 1
}

run_unadjusted <- function(x, cancer, method, feature, design_id, spec, analysis_set, threshold) {
  base <- list(Cancer = cancer, Method = method, Feature = feature,
               Broad_Compartment = broad_compartment(method, feature), Design_ID = design_id,
               Design_Type = spec$type, Analysis_Set = analysis_set, Model_Type = "UNADJUSTED",
               Target = spec$target, Reference = spec$reference, N = nrow(x), Threshold_Detail = threshold$detail)
  if (!threshold$ok) {
    return(list(type = spec$type, row = as.data.table(c(base, list(Status = "NOT_ESTIMABLE_BY_CONTRACT")))))
  }
  if (spec$type == "omnibus") {
    fit <- tryCatch(kruskal.test(y ~ exposure, data = x), error = function(e) NULL)
    if (is.null(fit)) return(list(type = spec$type, row = as.data.table(c(base, list(Status = "MODEL_FAILED")))))
    h <- unname(fit$statistic); k <- nlevels(x$exposure)
    eps <- max(0, (h - k + 1) / (nrow(x) - k))
    centers <- x[, .(Median = median(y), Mean = mean(y)), by = exposure]
    return(list(type = spec$type, row = as.data.table(c(base, list(
      Status = "ESTIMATED", Statistic = h, DF = unname(fit$parameter),
      Epsilon_Squared = eps, P_Value = fit$p.value,
      State_Centers = paste(centers$exposure, sprintf("%.6g", centers$Median), sep = ":", collapse = "|")
    )))))
  }
  if (spec$type == "pair") {
    target <- x[exposure == spec$target, y]
    reference <- x[exposure == spec$reference, y]
    fit <- tryCatch(wilcox.test(target, reference, exact = FALSE, correct = TRUE), error = function(e) NULL)
    effect <- hedges_g(target, reference)
    if (is.null(fit) || !is.finite(effect[["g"]])) {
      return(list(type = spec$type, row = as.data.table(c(base, list(Status = "MODEL_FAILED")))))
    }
    return(list(type = spec$type, row = as.data.table(c(base, list(
      Status = "ESTIMATED", N_Target = length(target), N_Reference = length(reference),
      Median_Target = median(target), Median_Reference = median(reference),
      Median_Difference = median(target) - median(reference), Rank_Biserial = rank_biserial(target, reference),
      Hedges_g = effect[["g"]], Hedges_g_SE = effect[["se"]], P_Value = fit$p.value
    )))))
  }
  fit <- tryCatch(suppressWarnings(cor.test(x$exposure, x$y, method = "spearman", exact = FALSE)),
                  error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$estimate)) {
    return(list(type = spec$type, row = as.data.table(c(base, list(Status = "MODEL_FAILED")))))
  }
  rho <- max(-0.999999, min(0.999999, unname(fit$estimate)))
  list(type = spec$type, row = as.data.table(c(base, list(
    Status = "ESTIMATED", Spearman_rho = rho, Fisher_z = atanh(rho),
    Fisher_z_SE = 1 / sqrt(nrow(x) - 3), P_Value = fit$p.value
  ))))
}

empty_adjusted <- function(cancer, method, feature, design_id, spec, analysis_set, model_type,
                           x, status, reason, complete_rate = NA_real_) {
  data.table(Cancer = cancer, Method = method, Feature = feature,
             Broad_Compartment = broad_compartment(method, feature), Design_ID = design_id,
             Design_Type = spec$type, Analysis_Set = analysis_set, Model_Type = model_type,
             Effect_Level = NA_character_, Target = spec$target, Reference = spec$reference,
             N_Input = nrow(x), N_Used = NA_integer_, Complete_Rate = complete_rate,
             Estimate = NA_real_, SE_HC3 = NA_real_, P_Value = NA_real_, Global_P = NA_real_,
             Status = status, Failure_Reason = reason)
}

run_adjusted <- function(x, cancer, method, feature, design_id, spec, analysis_set, model_type) {
  covars <- c("age10", "sex")
  if (model_type == "CORE_PURITY") covars <- c(covars, "purity_z")
  needed <- c("y_rint", "exposure", covars)
  comp <- x[complete.cases(x[, ..needed])]
  complete_rate <- if (nrow(x) == 0L) 0 else nrow(comp) / nrow(x)
  if (complete_rate < 0.70) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "NOT_ESTIMABLE_BY_CONTRACT", "COMPLETE_RATE_LT_0.70", complete_rate))
  }
  comp[, sex := droplevels(factor(sex))]
  if (nlevels(comp$sex) < 2L) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "NOT_ESTIMABLE_BY_CONTRACT", "SEX_SINGLE_LEVEL", complete_rate))
  }
  if (spec$type != "continuous") {
    comp[, exposure := droplevels(exposure)]
    tab <- table(comp$exposure)
    if (length(tab) != length(spec$levels) || any(tab < 10L)) {
      return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                            x, "NOT_ESTIMABLE_BY_CONTRACT", "ADJUSTED_GROUP_N_LT10", complete_rate))
    }
  } else if (nrow(comp) < 50L) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "NOT_ESTIMABLE_BY_CONTRACT", "ADJUSTED_CONTINUOUS_N_LT50", complete_rate))
  }
  form <- as.formula(paste("y_rint ~", paste(c("exposure", covars), collapse = " + ")))
  mm <- tryCatch(model.matrix(form, data = comp), error = function(e) NULL)
  if (is.null(mm) || nrow(comp) < max(40L, 10L * ncol(mm))) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "NOT_ESTIMABLE_BY_CONTRACT", "N_LT_MAX40_OR_10_PER_DF", complete_rate))
  }
  fit <- tryCatch(lm(form, data = comp), error = function(e) NULL)
  if (is.null(fit)) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "MODEL_FAILED", "LM_FAILED", complete_rate))
  }
  vc <- tryCatch(suppressWarnings(sandwich::vcovHC(fit, type = "HC3")), error = function(e) NULL)
  if (is.null(vc)) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "MODEL_FAILED", "HC3_FAILED", complete_rate))
  }
  cf <- coef(fit)
  effect_names <- grep("^exposure", names(cf), value = TRUE)
  effect_names <- effect_names[is.finite(cf[effect_names])]
  if (length(effect_names) == 0L) {
    return(empty_adjusted(cancer, method, feature, design_id, spec, analysis_set, model_type,
                          x, "MODEL_FAILED", "NO_EXPOSURE_COEFFICIENT", complete_rate))
  }
  vv <- vc[effect_names, effect_names, drop = FALSE]
  global_stat <- tryCatch(as.numeric(t(cf[effect_names]) %*% solve(vv, cf[effect_names])), error = function(e) NA_real_)
  global_p <- if (is.finite(global_stat)) pchisq(global_stat, df = length(effect_names), lower.tail = FALSE) else NA_real_
  out <- rbindlist(lapply(effect_names, function(nm) {
    se <- sqrt(vc[nm, nm])
    data.table(Cancer = cancer, Method = method, Feature = feature,
               Broad_Compartment = broad_compartment(method, feature), Design_ID = design_id,
               Design_Type = spec$type, Analysis_Set = analysis_set, Model_Type = model_type,
               Effect_Level = if (spec$type == "continuous") "per_1SD" else sub("^exposure", "", nm),
               Target = spec$target, Reference = spec$reference, N_Input = nrow(x), N_Used = nrow(comp),
               Complete_Rate = complete_rate, Estimate = unname(cf[nm]), SE_HC3 = se,
               CI_Lower = unname(cf[nm]) - 1.96 * se, CI_Upper = unname(cf[nm]) + 1.96 * se,
               P_Value = 2 * pnorm(abs(cf[nm] / se), lower.tail = FALSE), Global_P = global_p,
               Status = "ESTIMATED", Failure_Reason = NA_character_)
  }), fill = TRUE)
  # 高杠杆小分层可能使 HC3 标准误不再有限；显式降级，禁止进入后续 meta。
  out[!is.finite(SE_HC3) | SE_HC3 <= 0 | !is.finite(P_Value),
      `:=`(Status = "NUMERICALLY_UNSTABLE_HC3",
           Failure_Reason = "NONFINITE_OR_NONPOSITIVE_HC3")]
  out
}

reml_meta <- function(yi, sei) {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]; vi <- sei[ok]^2
  k <- length(yi)
  if (k < 2L) return(NULL)
  objective <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  upper <- max(10, var(yi) * 20, max(vi) * 20)
  opt <- optimize(objective, interval = c(0, upper))
  tau2 <- if (objective(0) <= opt$objective + 1e-10) 0 else opt$minimum
  w <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  se_mu <- sqrt(1 / sum(w))
  w_fixed <- 1 / vi
  mu_fixed <- sum(w_fixed * yi) / sum(w_fixed)
  q <- sum(w_fixed * (yi - mu_fixed)^2)
  q_df <- k - 1L
  list(k = k, beta = mu, se = se_mu, p = 2 * pnorm(abs(mu / se_mu), lower.tail = FALSE),
       tau2 = tau2, Q = q, Q_df = q_df, Q_p = pchisq(q, q_df, lower.tail = FALSE),
       I2 = ifelse(q <= 0, 0, max(0, (q - q_df) / q) * 100),
       pred_low = mu - 1.96 * sqrt(tau2 + se_mu^2), pred_high = mu + 1.96 * sqrt(tau2 + se_mu^2))
}

run_meta <- function(meta_input) {
  groups <- unique(meta_input[, .(Method, Feature, Broad_Compartment, Design_ID, Design_Type,
                                  Analysis_Set, Model_Type, Metric)])
  meta_rows <- list(); loo_rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- meta_input[g, on = names(g)]
    if (uniqueN(z$Cancer) < 5L) next
    m <- reml_meta(z$Beta, z$SE)
    if (is.null(m) || m$k < 5L) next
    meta_rows[[length(meta_rows) + 1L]] <- data.table(
      Method = g$Method, Feature = g$Feature, Broad_Compartment = g$Broad_Compartment,
      Design_ID = g$Design_ID, Design_Type = g$Design_Type, Analysis_Set = g$Analysis_Set,
      Model_Type = g$Model_Type, Metric = g$Metric, K_Cancers = m$k, Beta_REML = m$beta,
      SE_REML = m$se, CI_Lower = m$beta - 1.96 * m$se, CI_Upper = m$beta + 1.96 * m$se,
      P_Value = m$p, Tau2 = m$tau2, Q = m$Q, Q_DF = m$Q_df, Q_P = m$Q_p, I2 = m$I2,
      Prediction_Lower = m$pred_low, Prediction_Upper = m$pred_high,
      Direction_Consistency = mean(sign(z$Beta) == sign(m$beta)),
      Cancers = paste(sort(unique(z$Cancer)), collapse = ";"))
    for (left_out in unique(z$Cancer)) {
      ml <- reml_meta(z[Cancer != left_out]$Beta, z[Cancer != left_out]$SE)
      if (is.null(ml)) next
      loo_rows[[length(loo_rows) + 1L]] <- data.table(
        Method = g$Method, Feature = g$Feature, Broad_Compartment = g$Broad_Compartment,
        Design_ID = g$Design_ID, Design_Type = g$Design_Type, Analysis_Set = g$Analysis_Set,
        Model_Type = g$Model_Type, Metric = g$Metric, Left_Out_Cancer = left_out,
        K_Cancers = ml$k, Beta_REML = ml$beta, SE_REML = ml$se,
        CI_Lower = ml$beta - 1.96 * ml$se, CI_Upper = ml$beta + 1.96 * ml$se,
        P_Value = ml$p, I2 = ml$I2)
    }
  }
  meta <- rbindlist(meta_rows, fill = TRUE)
  loo <- rbindlist(loo_rows, fill = TRUE)
  if (nrow(meta) > 0L) {
    meta[, Meta_FDR := p.adjust(P_Value, method = "BH"),
         by = .(Method, Design_ID, Analysis_Set, Model_Type)]
    meta[, LOO_Direction_Stable := {
      ll <- loo[.BY, on = .(Method, Feature, Design_ID, Analysis_Set, Model_Type, Metric)]
      nrow(ll) > 0L && all(sign(ll$Beta_REML) == sign(Beta_REML))
    }, by = .(Method, Feature, Design_ID, Analysis_Set, Model_Type, Metric)]
    meta[, LOO_All_Nominal_Significant := {
      ll <- loo[.BY, on = .(Method, Feature, Design_ID, Analysis_Set, Model_Type, Metric)]
      nrow(ll) > 0L && all(ll$P_Value < 0.05)
    }, by = .(Method, Feature, Design_ID, Analysis_Set, Model_Type, Metric)]
    meta[Metric == "FISHER_Z", `:=`(
      Effect_REML_Interpretable = tanh(Beta_REML),
      Effect_CI_Lower_Interpretable = tanh(CI_Lower),
      Effect_CI_Upper_Interpretable = tanh(CI_Upper))]
    meta[Metric != "FISHER_Z", `:=`(
      Effect_REML_Interpretable = Beta_REML,
      Effect_CI_Lower_Interpretable = CI_Lower,
      Effect_CI_Upper_Interpretable = CI_Upper)]
  }
  list(meta = meta, loo = loo)
}

write_table <- function(x, stem) {
  path <- file.path(output_dir, paste0(stem, "_", version_tag, ".tsv"))
  fwrite(x, path, sep = "\t", na = "NA", quote = TRUE)
  path
}

main <- function() {
  start_time <- Sys.time()
  log_event("START", "Stage 2 immune ecology formal analysis begin",
            input = "canonical+stage1_patient+reused_scores", output = output_dir,
            parameters = ifelse(smoke_mode, "SMOKE", "FORMAL"))
  required_files <- c(canonical_file, patient_file, estimate_file, cibersort_file,
                      mcpcounter_file, acid_gene_file, lm12_gene_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) fail(paste("missing files", paste(missing_files, collapse = ";")))

  patient <- fread(patient_file, na.strings = c("NA", "", "NaN"))
  required_patient <- c("Patient_ID", "Sample", "Cancer", "C_subtype", "AGIG", "Acid_z", "IG_z",
                        "discordant_multi_sample_state", "age10", "sex", "purity_z")
  miss_patient <- setdiff(required_patient, names(patient))
  if (length(miss_patient) > 0L) fail(paste("patient manifest missing", paste(miss_patient, collapse = ",")))
  patient[, Sample := substr(as.character(Sample), 1, 15)]
  if (anyDuplicated(patient$Patient_ID) || anyDuplicated(patient$Sample)) fail("patient manifest is not unique")
  patient[, discordant_multi_sample_state := as.logical(discordant_multi_sample_state)]
  patient[, age10 := safe_numeric(age10)]
  patient[, purity_z := safe_numeric(purity_z)]
  patient[, Acid_z := safe_numeric(Acid_z)]
  patient[, IG_z := safe_numeric(IG_z)]
  patient[, sex := factor(sex)]
  if (length(selected_cancers) > 0L) patient <- patient[Cancer %in% selected_cancers]
  if (nrow(patient) == 0L) fail("no patients after cancer selection")

  canonical <- fread(canonical_file, select = c("Sample", "Cancer", "C_subtype", "AGIG"))
  canonical[, Sample := substr(as.character(Sample), 1, 15)]
  if (anyDuplicated(canonical$Sample)) fail("canonical labels are not unique")
  cc <- canonical[patient, on = "Sample", nomatch = 0L]
  if (nrow(cc) != nrow(patient) || any(cc$Cancer != cc$i.Cancer) ||
      any(cc$C_subtype != cc$i.C_subtype) || any(cc$AGIG != cc$i.AGIG)) {
    fail("patient interface does not match canonical labels")
  }

  patient_keep <- patient[, .(Patient_ID, Sample, Cancer, C_subtype, AGIG, Acid_z, IG_z,
                              discordant_multi_sample_state, age10, sex, purity_z)]
  scores <- list(
    ESTIMATE = read_method_scores("ESTIMATE", estimate_file, estimate_features, patient_keep),
    CIBERSORT_ABS = read_method_scores("CIBERSORT_ABS", cibersort_file, cibersort_features, patient_keep),
    MCPCOUNTER = read_method_scores("MCPCOUNTER", mcpcounter_file, mcpcounter_features, patient_keep)
  )
  feature_map <- list(ESTIMATE = estimate_features, CIBERSORT_ABS = cibersort_features, MCPCOUNTER = mcpcounter_features)
  log_event("CHECKPOINT", "canonical score interfaces validated",
            n_patients = nrow(patient), n_cancers = uniqueN(patient$Cancer),
            estimate_n = nrow(scores$ESTIMATE), cibersort_n = nrow(scores$CIBERSORT_ABS),
            mcpcounter_n = nrow(scores$MCPCOUNTER))

  input_registry <- data.table(
    Input_ID = c("CANONICAL_LABEL", "STAGE1_PATIENT", "ESTIMATE", "CIBERSORT_ABS", "MCPCOUNTER",
                 "ACID_GENES", "LM12_GENES"),
    Path = required_files, Exists = file.exists(required_files),
    Size_Bytes = file.info(required_files)$size, MD5 = vapply(required_files, md5_safe, character(1)),
    Access_Mode = "READ_ONLY", Reuse_Boundary = c(
      "CANONICAL_ONLY", "PATIENT_INTERFACE_ONLY", "DROP_OLD_LABELS_REUSE_SCORES",
      "REUSE_LABEL_INDEPENDENT_SCORES", "DROP_OLD_LABELS_REUSE_SCORES",
      "DEFINITION_AUDIT_ONLY", "DEFINITION_AUDIT_ONLY"))

  feature_registry <- make_signature_registry()

  sample_manifest <- copy(patient_keep)
  for (method in names(scores)) {
    key <- scores[[method]][, .(Sample, Available = TRUE)]
    setnames(key, "Available", paste0("Available_", method))
    sample_manifest <- merge(sample_manifest, key, by = "Sample", all.x = TRUE, sort = FALSE)
    sample_manifest[is.na(get(paste0("Available_", method))), (paste0("Available_", method)) := FALSE]
  }
  cib_qc <- scores$CIBERSORT_ABS[, .(Sample, CIBERSORT_P = `P-value`, CIBERSORT_Correlation = Correlation,
                                     CIBERSORT_RMSE = RMSE)]
  sample_manifest <- merge(sample_manifest, cib_qc, by = "Sample", all.x = TRUE, sort = FALSE)

  retention_rows <- list()
  for (patient_set in c("REPRESENTATIVE_MAIN", "EXCLUDE_DISCORDANT")) {
    z <- copy(scores$CIBERSORT_ABS)
    if (patient_set == "EXCLUDE_DISCORDANT") z <- z[discordant_multi_sample_state == FALSE]
    for (cancer in sort(unique(z$Cancer))) {
      zz <- z[Cancer == cancer]
      tab <- table(zz$C_subtype, zz$`P-value` < 0.05)
      p <- tryCatch({
        cs <- suppressWarnings(chisq.test(tab))
        if (any(cs$expected < 5)) fisher.test(tab)$p.value else cs$p.value
      }, error = function(e) NA_real_)
      by_state <- zz[, .(N_Total = .N, N_Confident = sum(`P-value` < 0.05, na.rm = TRUE)), by = C_subtype]
      by_state[, `:=`(Cancer = cancer, Patient_Set = patient_set,
                      Retention_Rate = N_Confident / N_Total, State_Dependency_P = p)]
      retention_rows[[length(retention_rows) + 1L]] <- by_state
    }
  }
  retention <- rbindlist(retention_rows, fill = TRUE)
  retention_p <- unique(retention[, .(Cancer, Patient_Set, State_Dependency_P)])
  retention_p[, State_Dependency_FDR := p.adjust(State_Dependency_P, method = "BH"), by = Patient_Set]
  retention <- merge(retention, retention_p, by = c("Cancer", "Patient_Set", "State_Dependency_P"), all.x = TRUE)

  omnibus_rows <- list(); pair_rows <- list(); continuous_rows <- list(); adjusted_rows <- list(); coverage_rows <- list()
  state_center_rows <- list()
  for (method in names(scores)) {
    method_dt <- scores[[method]]
    algorithm_sets <- if (method == "CIBERSORT_ABS") c("ALL_MATCHED", "CONFIDENT_P_LT_0.05") else "ALL_MATCHED"
    for (patient_set in c("REPRESENTATIVE_MAIN", "EXCLUDE_DISCORDANT")) {
      for (algorithm_set in algorithm_sets) {
        variant <- copy(method_dt)
        if (patient_set == "EXCLUDE_DISCORDANT") variant <- variant[discordant_multi_sample_state == FALSE]
        if (algorithm_set == "CONFIDENT_P_LT_0.05") variant <- variant[`P-value` < 0.05]
        analysis_set <- paste(patient_set, algorithm_set, sep = "__")
        for (cancer in sort(unique(variant$Cancer))) {
          cancer_dt <- variant[Cancer == cancer]
          for (feature in feature_map[[method]]) {
            if (patient_set == "REPRESENTATIVE_MAIN" && algorithm_set == "ALL_MATCHED") {
              centers <- cancer_dt[is.finite(get(feature)), .(
                N = .N, Median_Raw = median(get(feature)), Mean_Raw = mean(get(feature))), by = C_subtype]
              if (nrow(centers) > 0L) {
                tmp <- cancer_dt[is.finite(get(feature)), .(C_subtype, y = get(feature))]
                tmp[, y_rint := rint_safe(y)]
                rc <- tmp[, .(Median_RINT = median(y_rint), Mean_RINT = mean(y_rint)), by = C_subtype]
                centers <- merge(centers, rc, by = "C_subtype", all = TRUE)
                centers[, `:=`(Cancer = cancer, Method = method, Feature = feature,
                               Broad_Compartment = broad_compartment(method, feature))]
                state_center_rows[[length(state_center_rows) + 1L]] <- centers
              }
            }
            for (design_id in names(design_specs)) {
              spec <- design_specs[[design_id]]
              x <- prepare_design(cancer_dt, spec, feature)
              threshold <- check_threshold(x, spec)
              coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
                Cancer = cancer, Method = method, Feature = feature,
                Broad_Compartment = broad_compartment(method, feature), Design_ID = design_id,
                Design_Type = spec$type, Analysis_Set = analysis_set, N = nrow(x),
                Status = ifelse(threshold$ok, "ESTIMABLE", "NOT_ESTIMABLE_BY_CONTRACT"),
                Threshold_Detail = threshold$detail)
              u <- run_unadjusted(x, cancer, method, feature, design_id, spec, analysis_set, threshold)
              if (u$type == "omnibus") omnibus_rows[[length(omnibus_rows) + 1L]] <- u$row
              if (u$type == "pair") pair_rows[[length(pair_rows) + 1L]] <- u$row
              if (u$type == "continuous") continuous_rows[[length(continuous_rows) + 1L]] <- u$row
              model_types <- if (method == "ESTIMATE") "CORE_AGE_SEX" else c("CORE_AGE_SEX", "CORE_PURITY")
              for (model_type in model_types) {
                adjusted_rows[[length(adjusted_rows) + 1L]] <- run_adjusted(
                  x, cancer, method, feature, design_id, spec, analysis_set, model_type)
              }
            }
          }
        }
      }
    }
    log_event("CHECKPOINT", "method cancer-level models complete", method = method)
  }

  omnibus <- rbindlist(omnibus_rows, fill = TRUE)
  pairwise <- rbindlist(pair_rows, fill = TRUE)
  continuous <- rbindlist(continuous_rows, fill = TRUE)
  adjusted <- rbindlist(adjusted_rows, fill = TRUE)
  coverage <- rbindlist(coverage_rows, fill = TRUE)
  state_centers <- rbindlist(state_center_rows, fill = TRUE)

  omnibus[Status == "ESTIMATED", FDR := p.adjust(P_Value, method = "BH"),
          by = .(Cancer, Method, Design_ID, Analysis_Set)]
  pairwise[Status == "ESTIMATED", FDR := p.adjust(P_Value, method = "BH"),
           by = .(Cancer, Method, Design_ID, Analysis_Set)]
  continuous[Status == "ESTIMATED", FDR := p.adjust(P_Value, method = "BH"),
             by = .(Cancer, Method, Design_ID, Analysis_Set)]
  adjusted[Status == "ESTIMATED", FDR := p.adjust(P_Value, method = "BH"),
           by = .(Cancer, Method, Design_ID, Analysis_Set, Model_Type, Effect_Level)]
  adjusted_global <- unique(adjusted[Status == "ESTIMATED" & is.finite(Global_P),
                                     .(Cancer, Method, Feature, Design_ID, Analysis_Set, Model_Type, Global_P)])
  adjusted_global[, Global_FDR := p.adjust(Global_P, method = "BH"),
                  by = .(Cancer, Method, Design_ID, Analysis_Set, Model_Type)]
  adjusted <- merge(adjusted, adjusted_global,
                    by = c("Cancer", "Method", "Feature", "Design_ID", "Analysis_Set", "Model_Type", "Global_P"),
                    all.x = TRUE, sort = FALSE)

  state_profile <- state_centers[, .(
    K_Cancers = uniqueN(Cancer), Equal_Weight_Mean_Median_RINT = mean(Median_RINT, na.rm = TRUE),
    Equal_Weight_Median_Median_RINT = median(Median_RINT, na.rm = TRUE)),
    by = .(Method, Feature, Broad_Compartment, C_subtype)]

  meta_inputs <- list(
    pairwise[Status == "ESTIMATED" & is.finite(Hedges_g) & is.finite(Hedges_g_SE),
             .(Cancer, Method, Feature, Broad_Compartment, Design_ID, Design_Type, Analysis_Set,
               Model_Type = "UNADJUSTED", Metric = "HEDGES_G", Beta = Hedges_g, SE = Hedges_g_SE)],
    continuous[Status == "ESTIMATED" & is.finite(Fisher_z) & is.finite(Fisher_z_SE),
               .(Cancer, Method, Feature, Broad_Compartment, Design_ID, Design_Type, Analysis_Set,
                 Model_Type = "UNADJUSTED", Metric = "FISHER_Z", Beta = Fisher_z, SE = Fisher_z_SE)],
    adjusted[Status == "ESTIMATED" & Design_Type %in% c("pair", "continuous") &
               is.finite(Estimate) & is.finite(SE_HC3),
             .(Cancer, Method, Feature, Broad_Compartment, Design_ID, Design_Type, Analysis_Set,
               Model_Type, Metric = "RINT_BETA_HC3", Beta = Estimate, SE = SE_HC3)]
  )
  meta_input <- rbindlist(meta_inputs, fill = TRUE)
  meta_out <- run_meta(meta_input)

  output_paths <- c(
    write_table(input_registry, "immune_input_manifest"),
    write_table(feature_registry, "immune_input_and_feature_registry"),
    write_table(feature_registry, "immune_definition_circularity_audit"),
    write_table(sample_manifest, "stage2_sample_analysis_manifest"),
    write_table(coverage, "immune_coverage_and_evaluability_by_cancer"),
    write_table(retention, "cibersort_confidence_retention_audit"),
    write_table(omnibus, "immune_omnibus_by_cancer"),
    write_table(pairwise, "immune_pairwise_effects_by_cancer"),
    write_table(continuous, "immune_continuous_axis_associations_by_cancer"),
    write_table(adjusted, "immune_adjusted_models_by_cancer"),
    write_table(meta_out$meta, "immune_random_effects_meta"),
    write_table(meta_out$loo, "immune_leave_one_cancer_out"),
    write_table(state_centers, "immune_state_centers_by_cancer"),
    write_table(state_profile, "immune_state_profile_equal_weight")
  )
  manifest <- data.table(
    File = output_paths, Exists = file.exists(output_paths), Size_Bytes = file.info(output_paths)$size,
    MD5 = vapply(output_paths, md5_safe, character(1)), Stage = "STAGE2_CORE_ANALYSIS",
    Smoke_Mode = smoke_mode)
  write_table(manifest, "stage2_core_output_manifest")

  log_event("DONE", "Stage 2 immune ecology core analysis complete",
            runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
            n_omnibus = nrow(omnibus), n_pairwise = nrow(pairwise), n_continuous = nrow(continuous),
            n_adjusted = nrow(adjusted), n_meta = nrow(meta_out$meta), output = output_dir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1)
})
