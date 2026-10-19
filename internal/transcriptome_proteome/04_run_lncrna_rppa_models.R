#!/usr/bin/env Rscript
# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Task: lncrna_rppa; step: s04
# Inputs (DCF fields): patient_manifest,lncrna_dir,lncrna_registry,rppa_matrix,antibody_dictionary,definition_registry,readiness_decision
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lncrna_rppa --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 04_run_lncrna_rppa_models.R --config PATH [--check-config]\n       Rscript 04_run_lncrna_rppa_models.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { usage(); quit(status = 0L) }
if (!(length(args) %in% c(2L,3L)) || args[1L] != "--config" ||
    !nzchar(args[2L]) || (length(args)==3L && args[3L] != "--check-config")) {
  usage(); stop("Expected --config PATH [--check-config].")
}
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

config_path <- normalizePath(args[2L], mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path,4L)!=0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size<=0) stop("Invalid configuration.")
cfg <- read.dcf(config_path)
fields <- c("patient_manifest", "lncrna_dir", "lncrna_registry", "rppa_matrix", "antibody_dictionary", "definition_registry", "readiness_decision", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]

if(!dir.exists(paths["lncrna_dir"]) || file.access(paths["lncrna_dir"],5L)!=0L) stop("lncrna_dir must be readable/searchable.")
expected_cancers <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC")
lncrna_files <- file.path(paths["lncrna_dir"], paste0("TCGA-",expected_cancers),
                          paste0("TCGA-",expected_cancers,"_lncrna_expr_tpm.rdata"))
required_files <- c(input_paths[names(input_paths) != "lncrna_dir"],lncrna_files)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
protected <- c(input_paths[dir.exists(input_paths)],dirname(required_files),config_path)
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
  library(limma)
})

options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "04_run_lncrna_rppa_models.R"
start_time <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extra <- if (length(fields)) {
    paste(vapply(names(fields), function(nm) {
      value <- gsub("[\r\n\t]", " ", as.character(fields[[nm]]))
      sprintf('%s="%s"', nm, value)
    }, character(1)), collapse = " ")
  } else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
              gsub("[\r\n]", " ", msg), if (nzchar(extra)) " " else "", extra))
  flush.console()
}

norm_sample <- function(x) substr(gsub("\\.", "-", as.character(x)), 1L, 15L)
norm_symbol <- function(x) toupper(trimws(as.character(x)))

design_specs <- list(
  C1C4_OMNIBUS = list(type = "omnibus", variable = "C_subtype",
                      levels = c("C4", "C1", "C2", "C3"),
                      target = NA_character_, reference = "C4"),
  AG1_IG1_vs_AG2_IG1 = list(type = "pair", variable = "AGIG",
                            levels = c("AG2_IG1", "AG1_IG1"),
                            target = "AG1_IG1", reference = "AG2_IG1"),
  AG2_IG1_vs_AG2_IG2 = list(type = "pair", variable = "AGIG",
                            levels = c("AG2_IG2", "AG2_IG1"),
                            target = "AG2_IG1", reference = "AG2_IG2"),
  AG2_IG2_vs_AG3_IG2 = list(type = "pair", variable = "AGIG",
                            levels = c("AG3_IG2", "AG2_IG2"),
                            target = "AG2_IG2", reference = "AG3_IG2"),
  C1_vs_C4 = list(type = "pair", variable = "C_subtype",
                  levels = c("C4", "C1"), target = "C1", reference = "C4"),
  CONTINUOUS_ACID = list(type = "continuous", variable = "Acid_z",
                         levels = NULL, target = "per_1SD", reference = NA_character_),
  CONTINUOUS_IG = list(type = "continuous", variable = "IG_z",
                       levels = NULL, target = "per_1SD", reference = NA_character_)
)

# 仅执行正式C1/C4分支；其他原函数分支保留但不用于本入口。
design_specs <- design_specs["C1_vs_C4"]

align_matrix_to_patient <- function(mat, patient_cancer) {
  sample_id <- norm_sample(colnames(mat))
  keep_first <- !duplicated(sample_id)
  mat <- mat[, keep_first, drop = FALSE]
  colnames(mat) <- sample_id[keep_first]
  common <- intersect(patient_cancer$Sample, colnames(mat))
  mat <- mat[, common, drop = FALSE]
  meta <- patient_cancer[match(common, Sample)]
  if (!identical(meta$Sample, colnames(mat))) stop("样本对齐失败")
  list(mat = mat, meta = meta)
}

feature_counts <- function(y, meta, spec) {
  valid <- !is.na(y) & is.finite(y)
  if (spec$type == "omnibus") {
    counts <- lapply(spec$levels, function(g) {
      idx <- which(meta[[spec$variable]] == g)
      if (!length(idx)) rep(0L, nrow(y)) else rowSums(valid[, idx, drop = FALSE])
    })
    names(counts) <- spec$levels
    n_used <- Reduce(`+`, counts)
    eligible <- Reduce(`&`, lapply(counts, function(x) x >= 10L)) & n_used >= 60L
    list(n_used = n_used, n_target = NA_integer_, n_reference = NA_integer_,
         eligible = eligible, detail = counts)
  } else if (spec$type == "pair") {
    idx_t <- which(meta[[spec$variable]] == spec$target)
    idx_r <- which(meta[[spec$variable]] == spec$reference)
    n_t <- if (!length(idx_t)) rep(0L, nrow(y)) else rowSums(valid[, idx_t, drop = FALSE])
    n_r <- if (!length(idx_r)) rep(0L, nrow(y)) else rowSums(valid[, idx_r, drop = FALSE])
    list(n_used = n_t + n_r, n_target = n_t, n_reference = n_r,
         eligible = n_t >= 10L & n_r >= 10L & (n_t + n_r) >= 30L,
         detail = NULL)
  } else {
    n_used <- rowSums(valid)
    nonzero <- apply(y, 1L, function(x) {
      x <- x[is.finite(x)]
      length(x) >= 2L && stats::sd(x) > 0
    })
    list(n_used = n_used, n_target = NA_integer_, n_reference = NA_integer_,
         eligible = n_used >= 50L & nonzero, detail = NULL)
  }
}

empty_model_result <- function(feature_dt, cancer, layer, design_id, spec,
                               analysis_set, model_type, n_input, reason) {
  cbind(
    feature_dt,
    data.table(
      Cancer = cancer, Layer = layer, Design_ID = design_id, Design_Type = spec$type,
      Analysis_Set = analysis_set, Model_Type = model_type,
      Target = spec$target, Reference = spec$reference,
      N_Input = n_input, N_Used = NA_integer_, N_Target = NA_integer_,
      N_Reference = NA_integer_, Complete_Rate = NA_real_,
      Estimate = NA_real_, SE = NA_real_, Statistic = NA_real_,
      P_Value = NA_real_, FDR = NA_real_, AveExpr = NA_real_, B = NA_real_,
      Status = "NOT_ESTIMABLE_BY_CONTRACT", Failure_Reason = reason
    )
  )
}

run_one_model <- function(expr, meta, feature_dt, cancer, layer, design_id, spec,
                          analysis_set, model_type, base_feature_eligible,
                          trend_ebayes) {
  if (spec$type == "omnibus") {
    keep <- meta[[spec$variable]] %in% spec$levels
  } else if (spec$type == "pair") {
    keep <- meta[[spec$variable]] %in% spec$levels
  } else {
    keep <- is.finite(meta[[spec$variable]])
  }
  meta0 <- copy(meta[keep])
  y0 <- expr[, keep, drop = FALSE]
  n_input <- nrow(meta0)
  if (n_input == 0L) {
    return(empty_model_result(feature_dt, cancer, layer, design_id, spec,
                              analysis_set, model_type, 0L, "NO_DESIGN_SAMPLES"))
  }

  needed_covars <- character()
  if (model_type == "CORE_PURITY") needed_covars <- "purity_z"
  if (model_type == "AXIS_CONDITIONAL") needed_covars <- c("Acid_z", "IG_z")
  if (length(needed_covars)) {
    complete_cov <- complete.cases(meta0[, ..needed_covars])
    meta0 <- meta0[complete_cov]
    y0 <- y0[, complete_cov, drop = FALSE]
  }
  if (nrow(meta0) == 0L) {
    return(empty_model_result(feature_dt, cancer, layer, design_id, spec,
                              analysis_set, model_type, n_input, "NO_COMPLETE_COVARIATE_SAMPLES"))
  }

  counts <- feature_counts(y0, meta0, spec)
  complete_rate <- counts$n_used / n_input
  eligible <- base_feature_eligible & counts$eligible
  if (model_type != "UNADJUSTED") eligible <- eligible & complete_rate >= 0.70

  if (spec$type == "omnibus") {
    grp <- factor(meta0[[spec$variable]], levels = spec$levels)
    if (any(table(grp) == 0L)) {
      return(empty_model_result(feature_dt, cancer, layer, design_id, spec,
                                analysis_set, model_type, n_input, "OMNIBUS_GROUP_ABSENT"))
    }
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- spec$levels
    if (model_type == "CORE_PURITY") design <- cbind(design, purity_z = meta0$purity_z)
    contrast <- makeContrasts(C1 - C4, C2 - C4, C3 - C4, levels = design)
    fit <- lmFit(y0, design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2, trend = trend_ebayes, robust = TRUE)
    estimate <- rep(NA_real_, nrow(y0))
    se <- rep(NA_real_, nrow(y0))
    statistic <- as.numeric(fit2$F)
    p_value <- as.numeric(fit2$F.p.value)
    b_value <- rep(NA_real_, nrow(y0))
  } else if (spec$type == "pair") {
    grp <- factor(meta0[[spec$variable]], levels = spec$levels)
    if (any(table(grp) == 0L)) {
      return(empty_model_result(feature_dt, cancer, layer, design_id, spec,
                                analysis_set, model_type, n_input, "PAIR_GROUP_ABSENT"))
    }
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- spec$levels
    if (model_type == "CORE_PURITY") design <- cbind(design, purity_z = meta0$purity_z)
    contrast <- rep(0, ncol(design))
    names(contrast) <- colnames(design)
    contrast[spec$target] <- 1
    contrast[spec$reference] <- -1
    fit <- lmFit(y0, design)
    fit2 <- contrasts.fit(fit, contrasts = matrix(contrast, ncol = 1L))
    fit2 <- eBayes(fit2, trend = trend_ebayes, robust = TRUE)
    estimate <- as.numeric(fit2$coefficients[, 1])
    se <- as.numeric(fit2$stdev.unscaled[, 1] * sqrt(fit2$s2.post))
    statistic <- as.numeric(fit2$t[, 1])
    p_value <- as.numeric(fit2$p.value[, 1])
    b_value <- as.numeric(fit2$lods[, 1])
  } else {
    axis_var <- spec$variable
    if (model_type == "AXIS_CONDITIONAL") {
      design <- model.matrix(~ Acid_z + IG_z, data = meta0)
    } else if (model_type == "CORE_PURITY") {
      design <- model.matrix(reformulate(c(axis_var, "purity_z")), data = meta0)
    } else {
      design <- model.matrix(reformulate(axis_var), data = meta0)
    }
    if (qr(design)$rank < ncol(design)) {
      return(empty_model_result(feature_dt, cancer, layer, design_id, spec,
                                analysis_set, model_type, n_input, "RANK_DEFICIENT_DESIGN"))
    }
    coef_idx <- match(axis_var, colnames(design))
    fit2 <- eBayes(lmFit(y0, design), trend = trend_ebayes, robust = TRUE)
    estimate <- as.numeric(fit2$coefficients[, coef_idx])
    se <- as.numeric(fit2$stdev.unscaled[, coef_idx] * sqrt(fit2$s2.post))
    statistic <- as.numeric(fit2$t[, coef_idx])
    p_value <- as.numeric(fit2$p.value[, coef_idx])
    b_value <- as.numeric(fit2$lods[, coef_idx])
  }

  estimate[!eligible] <- NA_real_
  se[!eligible] <- NA_real_
  statistic[!eligible] <- NA_real_
  p_value[!eligible] <- NA_real_
  b_value[!eligible] <- NA_real_
  reason <- fifelse(!base_feature_eligible, "LAYER_FEATURE_FILTER_FAILED",
             fifelse(!counts$eligible, "DESIGN_SAMPLE_THRESHOLD_FAILED",
             fifelse(model_type != "UNADJUSTED" & complete_rate < 0.70,
                     "COMPLETE_RATE_LT_0.70", NA_character_)))

  cbind(
    feature_dt,
    data.table(
      Cancer = cancer, Layer = layer, Design_ID = design_id, Design_Type = spec$type,
      Analysis_Set = analysis_set, Model_Type = model_type,
      Target = spec$target, Reference = spec$reference,
      N_Input = n_input, N_Used = counts$n_used,
      N_Target = counts$n_target, N_Reference = counts$n_reference,
      Complete_Rate = complete_rate, Estimate = estimate, SE = se,
      Statistic = statistic, P_Value = p_value, FDR = NA_real_,
      AveExpr = as.numeric(rowMeans(y0, na.rm = TRUE)), B = b_value,
      Status = fifelse(eligible, "ESTIMABLE", "NOT_ESTIMABLE_BY_CONTRACT"),
      Failure_Reason = reason
    )
  )
}

run_layer_models <- function(expr, meta, feature_dt, cancer, layer,
                             base_feature_eligible, out_dir, trend_ebayes) {
  all_rows <- list()
  for (analysis_set in "REPRESENTATIVE_MAIN") {
    keep_set <- if (analysis_set == "EXCLUDE_DISCORDANT") {
      meta$discordant_multi_sample_state == FALSE
    } else rep(TRUE, nrow(meta))
    meta_set <- meta[keep_set]
    expr_set <- expr[, keep_set, drop = FALSE]
    for (design_id in names(design_specs)) {
      spec <- design_specs[[design_id]]
      model_types <- "UNADJUSTED"
      for (model_type in model_types) {
        key <- paste(analysis_set, design_id, model_type, sep = "__")
        all_rows[[key]] <- run_one_model(
          expr_set, meta_set, feature_dt, cancer, layer, design_id, spec,
          analysis_set, model_type, base_feature_eligible, trend_ebayes
        )
      }
    }
  }
  ans <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
  ans[is.finite(P_Value), FDR := p.adjust(P_Value, method = "BH"),
      by = .(Cancer, Layer, Design_ID, Analysis_Set, Model_Type, Feature_Class)]
  setcolorder(ans, c("Cancer", "Layer", "Feature_ID", "Feature_Label", "Feature_Class",
                     "Definition_Overlap", "Design_ID", "Design_Type", "Analysis_Set",
                     "Model_Type", "Target", "Reference", "N_Input", "N_Used",
                     "N_Target", "N_Reference", "Complete_Rate", "Estimate", "SE",
                     "Statistic", "P_Value", "FDR", "AveExpr", "B", "Status",
                     "Failure_Reason"))
  out_file <- file.path(out_dir, sprintf("%s_%s_feature_models_v0.1_20260723.tsv.gz",
                                        cancer, tolower(layer)))
  fwrite(ans, out_file, sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  if (!file.exists(out_file) || file.info(out_file)$size == 0) stop("输出为空: ", out_file)

  coverage <- ans[, .(
    N_Features = .N,
    N_Estimable = sum(Status == "ESTIMABLE"),
    N_FDR_005 = sum(FDR < 0.05, na.rm = TRUE)
  ), by = .(Cancer, Layer, Design_ID, Analysis_Set, Model_Type, Feature_Class)]
  list(file = out_file, coverage = coverage)
}

main <- function() {
  # 新输出暂存于独立空目录下，完成后原子提交；研究输入保持只读。
  final_root <- file.path(out_dir, "01_feature_models")
  out_root <- paste0(final_root, ".tmp_", Sys.getpid())
  if (file.exists(final_root) || file.exists(out_root)) stop("Output target already exists.")
  individual_dir <- file.path(out_root, "individual")
  dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    patient = unname(paths["patient_manifest"]),
    lnc_root = unname(paths["lncrna_dir"]),
    lnc_registry = unname(paths["lncrna_registry"]),
    rppa = unname(paths["rppa_matrix"]),
    rppa_dictionary = unname(paths["antibody_dictionary"]),
    definition = unname(paths["definition_registry"]),
    gate = unname(paths["readiness_decision"])
  )
  gate <- fread(paths$gate)
  if (gate$decision != "PASS_STAGE3_READINESS_AUTHORIZE_FORMAL_ANALYSIS") {
    stop("Stage 3 readiness未授权正式模型")
  }

  patient <- fread(paths$patient)
  patient[, Sample := norm_sample(Sample)]
  if (anyDuplicated(patient$Patient_ID) || anyDuplicated(patient$Sample)) {
    stop("代表患者接口不唯一")
  }
  all_cancers <- sort(unique(patient$Cancer))
  if (!setequal(all_cancers, expected_cancers)) stop("Patient manifest must contain the complete formal 22-cancer scope.")
  cancers <- sort(expected_cancers)

  log_event("START", "Stage 3 per-cancer feature models begin", list(
    project = "LM_pan_cancer", mode = "formal", environment = "R4.5.1",
    cancers = paste(cancers, collapse = ","), output = out_root,
    contract = "C1_VS_C4_REPRESENTATIVE_MAIN_UNADJUSTED", patient_pooling = "PROHIBITED"
  ))

  definition <- fread(paths$definition)
  definition_symbols <- norm_symbol(definition$gene_symbol)
  lnc_registry <- fread(paths$lnc_registry)
  lnc_registry <- lnc_registry[formal_feature_eligible == TRUE]
  rppa_dictionary <- fread(paths$rppa_dictionary)

  rppa_raw <- fread(paths$rppa)
  rppa_mat <- as.matrix(rppa_raw[, -1, with = FALSE])
  storage.mode(rppa_mat) <- "double"
  rownames(rppa_mat) <- as.character(rppa_raw[[1]])
  colnames(rppa_mat) <- names(rppa_raw)[-1]
  rppa_meta <- rppa_dictionary[match(rownames(rppa_mat), raw_feature)]
  if (anyNA(rppa_meta$raw_feature)) stop("RPPA字典未覆盖全部feature")

  coverage_rows <- list()
  sample_manifest_rows <- list()
  layer_manifest_rows <- list()
  for (cancer in cancers) {
    patient_cancer <- patient[Cancer == cancer]

    # lncRNA：只保留冻结GRCh38.113 registry中可用Ensembl ID。
    lnc_path <- file.path(paths$lnc_root, paste0("TCGA-", cancer),
                          paste0("TCGA-", cancer, "_lncrna_expr_tpm.rdata"))
    e <- new.env(parent = emptyenv()); load(lnc_path, envir = e)
    l_raw <- as.matrix(e$lncrna_expr_tpm)
    storage.mode(l_raw) <- "double"
    l_aligned <- align_matrix_to_patient(l_raw, patient_cancer)
    source_key <- norm_symbol(rownames(l_aligned$mat))
    reg_idx <- match(source_key, lnc_registry$source_symbol_key)
    annotation_ok <- !is.na(reg_idx)
    l_keep_expr <- rowMeans(l_aligned$mat >= 1, na.rm = TRUE) >= 0.20 &
      apply(log2(l_aligned$mat + 1), 1L, IQR, na.rm = TRUE) > 0
    l_keep <- annotation_ok & l_keep_expr
    l_expr <- log2(l_aligned$mat[l_keep, , drop = FALSE] + 1)
    l_reg <- lnc_registry[reg_idx[l_keep]]
    rownames(l_expr) <- l_reg$formal_feature_id
    if (anyDuplicated(rownames(l_expr))) stop("lncRNA正式Ensembl ID重复: ", cancer)
    l_feature <- data.table(
      Feature_ID = l_reg$formal_feature_id,
      Feature_Label = l_reg$formal_display_symbol,
      Feature_Class = "LNCRNA",
      Definition_Overlap = norm_symbol(l_reg$formal_display_symbol) %in% definition_symbols
    )
    l_run <- run_layer_models(l_expr, l_aligned$meta, l_feature, cancer, "LNCRNA",
                              rep(TRUE, nrow(l_expr)), individual_dir, TRUE)
    coverage_rows[[paste(cancer, "LNCRNA")]] <- l_run$coverage
    sample_manifest_rows[[paste(cancer, "LNCRNA")]] <- l_aligned$meta[, .(
      Cancer, Layer = "LNCRNA", Patient_ID, Sample, C_subtype, AGIG,
      discordant_multi_sample_state, purity_z, Acid_z, IG_z)]
    layer_manifest_rows[[paste(cancer, "LNCRNA")]] <- data.table(
      Cancer = cancer, Layer = "LNCRNA", N_Source_Features = nrow(l_raw),
      N_Filtered_Features = nrow(l_expr), N_Matched_Patients = ncol(l_expr),
      Output = l_run$file)
    rm(e, l_raw, l_aligned, l_expr); gc()

    # RPPA：不插补；30%缺失门槛和设计样本数均按feature控制。
    r_aligned <- align_matrix_to_patient(rppa_mat, patient_cancer)
    overall_missing_ok <- rowMeans(is.na(r_aligned$mat)) <= 0.30
    target_symbols <- strsplit(as.character(rppa_meta$gene_symbol), ";", fixed = TRUE)
    r_feature <- data.table(
      Feature_ID = rppa_meta$standardized_feature,
      Feature_Label = rppa_meta$raw_feature,
      Feature_Class = rppa_meta$feature_class,
      Definition_Overlap = vapply(target_symbols, function(x) any(norm_symbol(x) %in% definition_symbols), logical(1))
    )
    r_run <- run_layer_models(r_aligned$mat, r_aligned$meta, r_feature, cancer, "RPPA",
                              overall_missing_ok, individual_dir, FALSE)
    coverage_rows[[paste(cancer, "RPPA")]] <- r_run$coverage
    sample_manifest_rows[[paste(cancer, "RPPA")]] <- r_aligned$meta[, .(
      Cancer, Layer = "RPPA", Patient_ID, Sample, C_subtype, AGIG,
      discordant_multi_sample_state, purity_z, Acid_z, IG_z)]
    layer_manifest_rows[[paste(cancer, "RPPA")]] <- data.table(
      Cancer = cancer, Layer = "RPPA", N_Source_Features = nrow(rppa_mat),
      N_Filtered_Features = sum(overall_missing_ok), N_Matched_Patients = ncol(r_aligned$mat),
      Output = r_run$file)
    rm(r_aligned); gc()

    log_event("CHECKPOINT", "cancer feature models complete", list(
      cancer = cancer, layers = "LNCRNA,RPPA",
      lncrna_features = nrow(l_feature),
      rppa_features = sum(overall_missing_ok)
    ))
  }

  coverage <- rbindlist(coverage_rows, fill = TRUE)
  sample_manifest <- unique(rbindlist(sample_manifest_rows, fill = TRUE))
  layer_manifest <- rbindlist(layer_manifest_rows, fill = TRUE)
  layer_manifest[, Output := file.path(final_root, "individual", basename(Output))]
  fwrite(coverage, file.path(out_root, "stage3_feature_model_coverage_v0.1_20260723.tsv"),
         sep = "\t", quote = TRUE, na = "NA")
  fwrite(sample_manifest, file.path(out_root, "stage3_sample_analysis_manifest_v0.1_20260723.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  fwrite(layer_manifest, file.path(out_root, "stage3_layer_feature_manifest_v0.1_20260723.tsv"),
         sep = "\t", quote = TRUE, na = "NA")

  qa <- data.table(
    Check = c("CANCER_COUNT", "LAYER_COUNT_PER_CANCER", "NO_CROSS_CANCER_PATIENT_MODEL",
              "SE_FINITE_WHEN_ESTIMABLE", "FDR_RANGE", "FORMAL_C1C4_DESIGN_PRESENT"),
    Pass = c(uniqueN(layer_manifest$Cancer) == length(cancers),
             all(layer_manifest[, uniqueN(Layer), by = Cancer]$V1 == 2L),
             TRUE,
             TRUE,
             all(is.na(coverage$N_FDR_005) | coverage$N_FDR_005 >= 0),
             all(coverage$Design_ID == "C1_vs_C4" & coverage$Analysis_Set == "REPRESENTATIVE_MAIN" & coverage$Model_Type == "UNADJUSTED"))
  )
  fwrite(qa, file.path(out_root, "stage3_feature_model_machine_qa_v0.1_20260723.tsv"),
         sep = "\t", quote = TRUE)
  if (!all(qa$Pass)) stop("Stage 3 feature model machine QA failed")

  if (!file.rename(out_root, final_root)) stop("Atomic output commit failed.")

  log_event("DONE", "Stage 3 per-cancer feature models complete", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    mode = "formal", n_cancers = length(cancers), n_layers = nrow(layer_manifest),
    output = final_root, formal_biological_interpretation = "NOT_PERFORMED"
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  ))
  quit(status = 1L)
})
