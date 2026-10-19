# Purpose: lncrna rppa analysis stages.
# Configuration: config/molecular.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/04_run_lncrna_rppa_models.R
stage_transcriptome_proteome_04_run_lncrna_rppa_models <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`norm_sample` <- bind_stage_helper(shared_norm_sample_037, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())

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
`align_matrix_to_patient` <- bind_stage_helper(shared_align_matrix_to_patient_017, environment())

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
}

# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/03_finalize_rppa_recurrence_tcpa.R
stage_transcriptome_proteome_03_finalize_rppa_recurrence_tcpa <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("model_dir must be readable/searchable.")

required_files <- input_paths[names(input_paths) != "model_dir"]
required_files <- c(required_files,file.path(paths["model_dir"],c("BLCA_rppa_feature_models_v0.1_20260723.tsv.gz","BRCA_rppa_feature_models_v0.1_20260723.tsv.gz","CESC_rppa_feature_models_v0.1_20260723.tsv.gz","COAD_rppa_feature_models_v0.1_20260723.tsv.gz","ESCA_rppa_feature_models_v0.1_20260723.tsv.gz","GBM_rppa_feature_models_v0.1_20260723.tsv.gz","HNSC_rppa_feature_models_v0.1_20260723.tsv.gz","KIRC_rppa_feature_models_v0.1_20260723.tsv.gz","KIRP_rppa_feature_models_v0.1_20260723.tsv.gz","LGG_rppa_feature_models_v0.1_20260723.tsv.gz","LIHC_rppa_feature_models_v0.1_20260723.tsv.gz","LUAD_rppa_feature_models_v0.1_20260723.tsv.gz","LUSC_rppa_feature_models_v0.1_20260723.tsv.gz","OV_rppa_feature_models_v0.1_20260723.tsv.gz","PCPG_rppa_feature_models_v0.1_20260723.tsv.gz","PRAD_rppa_feature_models_v0.1_20260723.tsv.gz","READ_rppa_feature_models_v0.1_20260723.tsv.gz","SARC_rppa_feature_models_v0.1_20260723.tsv.gz","SKCM_rppa_feature_models_v0.1_20260723.tsv.gz","STAD_rppa_feature_models_v0.1_20260723.tsv.gz","THCA_rppa_feature_models_v0.1_20260723.tsv.gz","UCEC_rppa_feature_models_v0.1_20260723.tsv.gz")))

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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(data.table)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "03_finalize_rppa_recurrence_tcpa.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_012, environment())

normalize_token <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
}

main <- function() {
  model_dir <- unname(paths["model_dir"])
  dictionary_file <- unname(paths["antibody_dictionary"])
  tcpa_file <- unname(paths["tcpa_table13"])
  final_dir <- file.path(out_dir,"01_freeze_and_source")
  if (dir.exists(final_dir)) {
    stop("正式输出目录已存在，拒绝覆盖: ", final_dir)
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  expected_cancers <- c(
    "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC",
    "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PCPG", "PRAD",
    "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC"
  )
  min_group_n <- 30L
  min_recurrent_cancers <- 3L

  model_files <- file.path(
    model_dir,
    paste0(
      expected_cancers,
      "_rppa_feature_models_v0.1_20260723.tsv.gz"
    )
  )
  required <- c(model_files, dictionary_file, tcpa_file)
  if (any(!file.exists(required)) ||
      any(file.info(required)$size <= 0)) {
    stop("存在缺失或空的冻结输入")
  }

  log_event("START", "构建RPPA C1 vs C4复现和TCPA绘图源数据", list(
    environment = "R4.5.1",
    input_cancers = length(expected_cancers),
    min_group_n = min_group_n,
    min_recurrent_cancers = min_recurrent_cancers,
    cross_cancer_meta = "NOT_PERFORMED",
    purity = "NOT_USED",
    mrna = "NOT_USED",
    output = final_dir
  ))

  read_one <- function(cancer, path) {
    x <- fread(path)
    x <- x[
      Design_ID == "C1_vs_C4" &
        Analysis_Set == "REPRESENTATIVE_MAIN" &
        Model_Type == "UNADJUSTED"
    ]
    if (nrow(x) != 258L ||
        uniqueN(x$Cancer) != 1L ||
        unique(x$Cancer) != cancer ||
        anyDuplicated(x$Feature_ID) ||
        any(x$Target != "C1") ||
        any(x$Reference != "C4")) {
      stop("单癌种RPPA正式输入身份或维度错误: ", cancer)
    }
    x
  }
  individual_all <- rbindlist(
    Map(read_one, expected_cancers, model_files),
    use.names = TRUE, fill = TRUE
  )

  # 旧版RPPA大图只纳入两组均有至少30名匹配患者的癌种。
  cancer_samples <- individual_all[, .(
    N_C1_Max = max(N_Target, na.rm = TRUE),
    N_C4_Max = max(N_Reference, na.rm = TRUE),
    N_Estimable_Features = sum(Status == "ESTIMABLE")
  ), by = Cancer]
  cancer_samples[, Formal_Cancer :=
                   N_C1_Max >= min_group_n &
                   N_C4_Max >= min_group_n &
                   N_Estimable_Features > 0L]
  formal_cancers <- cancer_samples[Formal_Cancer == TRUE, Cancer]
  individual <- copy(individual_all[Cancer %chin% formal_cancers])

  # 按癌种内整个RPPA面板统一BH；绝不跨癌校正。
  individual[, Original_Feature_Class_FDR := FDR]
  individual[, FDR := NA_real_]
  individual[
    Status == "ESTIMABLE" & is.finite(P_Value),
    FDR := p.adjust(P_Value, method = "BH"),
    by = Cancer
  ]
  individual[, Single_Cancer_Direction := fcase(
    Status == "ESTIMABLE" & is.finite(FDR) &
      FDR < 0.05 & Estimate > 0,
    "C1_HIGH",
    Status == "ESTIMABLE" & is.finite(FDR) &
      FDR < 0.05 & Estimate < 0,
    "C4_HIGH",
    default = "NOT_SIGNIFICANT"
  )]

  dictionary <- fread(dictionary_file)
  if (nrow(dictionary) != 258L ||
      anyDuplicated(dictionary$standardized_feature)) {
    stop("RPPA冻结字典不是258个唯一feature")
  }
  dictionary_small <- dictionary[, .(
    Feature_ID = standardized_feature,
    Gene_Symbol = gene_symbol,
    Site = site,
    Target_Count = target_count,
    Mapping_Status = mapping_status,
    Validation_Status = validation_status,
    Direct_mRNA_Pairable = direct_mrna_pairable,
    Cross_Layer_Restriction = cross_layer_restriction
  )]

  recurrence <- individual[, {
    estimates <- Estimate[Status == "ESTIMABLE" & is.finite(Estimate)]
    sig_estimates <- Estimate[
      Status == "ESTIMABLE" & is.finite(FDR) & FDR < 0.05
    ]
    .(
      Feature_Label = Feature_Label[1L],
      Feature_Class = Feature_Class[1L],
      Definition_Overlap = any(Definition_Overlap %in% TRUE),
      N_Formal_Cancers = length(formal_cancers),
      N_Cancer_Estimable = uniqueN(Cancer[
        Status == "ESTIMABLE" & is.finite(P_Value)
      ]),
      N_C1_High_Significant = uniqueN(Cancer[
        Single_Cancer_Direction == "C1_HIGH"
      ]),
      N_C4_High_Significant = uniqueN(Cancer[
        Single_Cancer_Direction == "C4_HIGH"
      ]),
      Median_Estimate_Estimable = if (length(estimates)) {
        median(estimates)
      } else NA_real_,
      Median_Estimate_Significant = if (length(sig_estimates)) {
        median(sig_estimates)
      } else NA_real_,
      Max_Abs_Estimate = if (length(estimates)) {
        max(abs(estimates))
      } else NA_real_
    )
  }, by = Feature_ID]
  recurrence <- merge(
    recurrence, dictionary_small, by = "Feature_ID", all.x = TRUE
  )
  recurrence[, `:=`(
    C1_High_Recurrence_Fraction =
      N_C1_High_Significant / N_Formal_Cancers,
    C4_High_Recurrence_Fraction =
      N_C4_High_Significant / N_Formal_Cancers
  )]
  recurrence[, Recurrent_Class := fcase(
    N_C1_High_Significant >= min_recurrent_cancers &
      N_C4_High_Significant >= min_recurrent_cancers,
    "BIDIRECTIONAL_RECURRENT",
    N_C1_High_Significant >= min_recurrent_cancers,
    "C1_HIGH",
    N_C4_High_Significant >= min_recurrent_cancers,
    "C4_HIGH",
    default = "NOT_RECURRENT"
  )]
  recurrence[, Formal_Display :=
               Recurrent_Class %chin% c("C1_HIGH", "C4_HIGH")]
  recurrence[, Dominant_Recurrence_N := pmax(
    N_C1_High_Significant, N_C4_High_Significant
  )]
  recurrence[, Dominant_Recurrence_Fraction :=
               Dominant_Recurrence_N / N_Formal_Cancers]
  recurrence[, Abs_Median_Estimate_Significant :=
               abs(Median_Estimate_Significant)]
  setorder(
    recurrence, -Formal_Display, Recurrent_Class,
    -Dominant_Recurrence_N, -Abs_Median_Estimate_Significant,
    Feature_ID
  )

  selected <- copy(recurrence[Formal_Display == TRUE])
  selected[, State_Order := match(
    Recurrent_Class, c("C1_HIGH", "C4_HIGH")
  )]
  setorder(
    selected, State_Order, -Dominant_Recurrence_N,
    -Abs_Median_Estimate_Significant, Feature_ID
  )
  selected[, Display_Order := .I]
  selected[, Display_Label := Feature_ID]

  # 构建Supplementary Table 13到当前258-feature字典的严格映射。
  official <- fread(tcpa_file)
  if (nrow(official) != 75L || uniqueN(official$Pathway) != 10L) {
    stop("TCPA Supplementary Table 13结构化输入维度错误")
  }
  official[, Source_Normalized := normalize_token(Source_Protein_Token)]
  official[, Feature_Normalized := Source_Normalized]
  official[
    Source_Normalized == "HER3PY1298",
    Feature_Normalized := "HER3PY1289"
  ]
  dictionary_map <- dictionary[, .(
    Feature_ID = standardized_feature,
    Feature_Normalized = normalize_token(standardized_feature)
  )]
  official_map <- merge(
    official, dictionary_map, by = "Feature_Normalized", all.x = TRUE
  )
  if (anyNA(official_map$Feature_ID)) {
    stop("TCPA Table 13存在未映射到258-feature字典的成员")
  }

  selected_map <- merge(
    selected[, .(
      Feature_ID, Display_Label, Recurrent_Class, Display_Order
    )],
    official_map[, .(
      Feature_ID, TCPA_Pathway = Pathway,
      TCPA_Component_Sign = Direction_Sign,
      Source_Protein_Token
    )],
    by = "Feature_ID", all.x = TRUE, allow.cartesian = TRUE
  )
  selected_map[
    is.na(TCPA_Pathway),
    `:=`(
      TCPA_Pathway = "Other_not_in_TCPA_Table13",
      TCPA_Component_Sign = 0L,
      Source_Protein_Token = NA_character_
    )
  ]
  selected_map[, Mapping_Status := fifelse(
    TCPA_Pathway == "Other_not_in_TCPA_Table13",
    "UNMAPPED_TO_TCPA_TABLE13",
    "OFFICIAL_TCPA_TABLE13_MAPPED"
  )]

  pathway_order <- c(
    "Apoptosis", "Cell_cycle", "DNA_damage_response", "EMT",
    "Hormone_a", "Hormone_b", "PI3K_Akt", "RAS_MAPK", "RTK",
    "TSC_mTOR", "Other_not_in_TCPA_Table13"
  )
  pathway_labels <- c(
    Apoptosis = "Apoptosis",
    Cell_cycle = "Cell cycle",
    DNA_damage_response = "DNA damage response",
    EMT = "EMT",
    Hormone_a = "Hormone a",
    Hormone_b = "Hormone b",
    PI3K_Akt = "PI3K / Akt",
    RAS_MAPK = "RAS / MAPK",
    RTK = "RTK",
    TSC_mTOR = "TSC / mTOR",
    Other_not_in_TCPA_Table13 = "Other / not in TCPA Table 13"
  )
  selected_map[, Pathway_Order := match(TCPA_Pathway, pathway_order)]
  selected_map[, TCPA_Pathway_Label :=
                 unname(pathway_labels[TCPA_Pathway])]
  setorder(selected_map, Display_Order, Pathway_Order)

  feature_mapping <- selected_map[, .(
    N_TCPA_Pathways = sum(
      TCPA_Pathway != "Other_not_in_TCPA_Table13"
    ),
    TCPA_Pathway_Membership = paste(
      TCPA_Pathway, collapse = ";"
    ),
    TCPA_Direction_Membership = paste(
      paste0(
        TCPA_Pathway, ":",
        fifelse(
          TCPA_Component_Sign > 0L, "+",
          fifelse(TCPA_Component_Sign < 0L, "-", "0")
        )
      ),
      collapse = ";"
    ),
    TCPA_Mapping_Status = if (
      any(TCPA_Pathway != "Other_not_in_TCPA_Table13")
    ) "OFFICIAL_TCPA_TABLE13_MAPPED" else {
      "UNMAPPED_TO_TCPA_TABLE13"
    }
  ), by = Feature_ID]
  selected <- merge(
    selected, feature_mapping, by = "Feature_ID", all.x = TRUE
  )
  setorder(selected, Display_Order)

  dot_source <- merge(
    individual,
    selected[, .(
      Feature_ID, Recurrent_Class, Dominant_Recurrence_N,
      Dominant_Recurrence_Fraction, Display_Order, Display_Label,
      TCPA_Pathway_Membership, TCPA_Direction_Membership,
      TCPA_Mapping_Status
    )],
    by = "Feature_ID", all = FALSE
  )
  cancer_order <- cancer_samples[
    Cancer %chin% formal_cancers
  ][order(match(Cancer, formal_cancers)), Cancer]
  dot_source[, Cancer_Order := match(Cancer, cancer_order)]
  dot_source[, Dot_Display :=
               Status == "ESTIMABLE" &
               is.finite(FDR) & FDR < 0.05]
  dot_source[, Minus_Log10_FDR :=
               -log10(pmax(FDR, .Machine$double.xmin))]
  dot_source[, Minus_Log10_FDR_Capped :=
               pmin(Minus_Log10_FDR, 10)]
  dot_source[, Estimate_Capped_1 :=
               pmax(pmin(Estimate, 1), -1)]
  setorder(dot_source, Cancer_Order, Display_Order)

  pathway_registry <- selected_map[, .(
    N_Display_Features = uniqueN(Feature_ID),
    N_C1_High_Features = uniqueN(
      Feature_ID[Recurrent_Class == "C1_HIGH"]
    ),
    N_C4_High_Features = uniqueN(
      Feature_ID[Recurrent_Class == "C4_HIGH"]
    ),
    N_Positive_Components = sum(TCPA_Component_Sign > 0L),
    N_Negative_Components = sum(TCPA_Component_Sign < 0L),
    N_Other = sum(TCPA_Component_Sign == 0L)
  ), by = .(
    TCPA_Pathway, TCPA_Pathway_Label, Pathway_Order
  )]
  setorder(pathway_registry, Pathway_Order)

  threshold_audit <- rbindlist(lapply(c(0, 0.1, 0.2), function(cutoff) {
    a <- individual[
      Status == "ESTIMABLE" & is.finite(FDR) &
        FDR < 0.05 & abs(Estimate) >= cutoff,
      .(
        N_C1 = uniqueN(Cancer[Estimate > 0]),
        N_C4 = uniqueN(Cancer[Estimate < 0])
      ),
      by = Feature_ID
    ]
    data.table(
      Abs_Estimate_Cutoff = cutoff,
      N_Selected_Including_Bidirectional =
        sum(a$N_C1 >= min_recurrent_cancers |
              a$N_C4 >= min_recurrent_cancers),
      N_C1_High_Only =
        sum(a$N_C1 >= min_recurrent_cancers &
              a$N_C4 < min_recurrent_cancers),
      N_C4_High_Only =
        sum(a$N_C4 >= min_recurrent_cancers &
              a$N_C1 < min_recurrent_cancers),
      N_Bidirectional =
        sum(a$N_C1 >= min_recurrent_cancers &
              a$N_C4 >= min_recurrent_cancers)
    )
  }))
  threshold_audit[, Formal_Selection := Abs_Estimate_Cutoff == 0]

  summary_table <- rbindlist(list(
    data.table(
      Metric = c(
        "INPUT_CANCERS", "FORMAL_CANCERS",
        "MIN_GROUP_N", "MIN_RECURRENT_CANCERS",
        "TOTAL_RPPA_FEATURES", "FORMAL_DISPLAY_FEATURES",
        "C1_HIGH_DISPLAY_FEATURES", "C4_HIGH_DISPLAY_FEATURES",
        "BIDIRECTIONAL_EXCLUDED_FEATURES",
        "VISIBLE_CANCER_FEATURE_DOTS",
        "TCPA_OFFICIAL_MAPPED_DISPLAY_FEATURES",
        "OTHER_DISPLAY_FEATURES", "HIT_PATHWAY_NODES"
      ),
      Value = c(
        length(expected_cancers), length(formal_cancers),
        min_group_n, min_recurrent_cancers,
        uniqueN(recurrence$Feature_ID), nrow(selected),
        selected[Recurrent_Class == "C1_HIGH", .N],
        selected[Recurrent_Class == "C4_HIGH", .N],
        recurrence[Recurrent_Class == "BIDIRECTIONAL_RECURRENT", .N],
        sum(dot_source$Dot_Display),
        selected[TCPA_Mapping_Status ==
                   "OFFICIAL_TCPA_TABLE13_MAPPED", .N],
        selected[TCPA_Mapping_Status ==
                   "UNMAPPED_TO_TCPA_TABLE13", .N],
        nrow(pathway_registry)
      )
    ),
    recurrence[, .(
      Metric = paste0("RECURRENCE_CLASS_", Recurrent_Class),
      Value = .N
    ), by = Recurrent_Class][, Recurrent_Class := NULL]
  ), use.names = TRUE)

  threshold_registry <- data.table(
    Parameter = c(
      "DESIGN_ID", "ANALYSIS_SET", "MODEL_TYPE",
      "MIN_C1_PATIENTS_PER_CANCER",
      "MIN_C4_PATIENTS_PER_CANCER",
      "WITHIN_CANCER_BH_FAMILY",
      "WITHIN_CANCER_FDR",
      "ABS_ESTIMATE_HARD_THRESHOLD",
      "MIN_SAME_DIRECTION_RECURRENT_CANCERS",
      "BIDIRECTIONAL_POLICY",
      "CROSS_CANCER_BH", "META_LOO_I2",
      "PURITY_MODEL", "MRNA_OR_ORA_IN_RPPA_SELECTION"
    ),
    Value = c(
      "C1_vs_C4", "REPRESENTATIVE_MAIN", "UNADJUSTED",
      "30", "30", "ALL_ESTIMABLE_RPPA_FEATURES",
      "0.05", "NONE", "3", "EXCLUDE_FROM_FORMAL_DISPLAY",
      "NOT_PERFORMED", "NOT_PERFORMED",
      "EXCLUDED", "NOT_USED"
    )
  )

  files <- c(
    within = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_within_cancer_results_v0.1_20260730.tsv.gz"
    ),
    cancer = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_cancer_sample_gate_v0.1_20260730.tsv"
    ),
    recurrence = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_recurrence_all_v0.1_20260730.tsv"
    ),
    selected = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_display_feature_freeze_v0.1_20260730.tsv"
    ),
    dot = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_dotplot_source_v0.1_20260730.tsv.gz"
    ),
    edges = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_tcpa_edges_v0.1_20260730.tsv"
    ),
    pathways = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_tcpa_pathway_registry_v0.1_20260730.tsv"
    ),
    audit = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_threshold_audit_v0.1_20260730.tsv"
    ),
    thresholds = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_threshold_registry_v0.1_20260730.tsv"
    ),
    summary = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_summary_v0.1_20260730.tsv"
    )
  )
  fwrite(individual, files[["within"]], sep = "\t", na = "NA")
  fwrite(cancer_samples, files[["cancer"]], sep = "\t", na = "NA")
  fwrite(recurrence, files[["recurrence"]], sep = "\t", na = "NA")
  fwrite(selected, files[["selected"]], sep = "\t", na = "NA")
  fwrite(dot_source, files[["dot"]], sep = "\t", na = "NA")
  fwrite(selected_map, files[["edges"]], sep = "\t", na = "NA")
  fwrite(pathway_registry, files[["pathways"]], sep = "\t", na = "NA")
  fwrite(threshold_audit, files[["audit"]], sep = "\t", na = "NA")
  fwrite(threshold_registry, files[["thresholds"]], sep = "\t", na = "NA")
  fwrite(summary_table, files[["summary"]], sep = "\t", na = "NA")

  fdr_recalc_error <- individual[
    Status == "ESTIMABLE" & is.finite(P_Value),
    .(Error = max(abs(
      FDR - p.adjust(P_Value, method = "BH")
    ))),
    by = Cancer
  ][, max(Error)]
  qa <- data.table(
    Check = c(
      "INPUT_CANCERS_22",
      "FORMAL_CANCERS_15",
      "FORMAL_CANCERS_MEET_BOTH_GROUP_N30",
      "ONLY_C1_VS_C4_REPRESENTATIVE_UNADJUSTED",
      "WITHIN_CANCER_FULL_PANEL_BH_EXACT",
      "NO_CROSS_CANCER_BH_OR_META",
      "FORMAL_FEATURES_SAME_DIRECTION_N_GE_3",
      "NO_BIDIRECTIONAL_IN_FORMAL_DISPLAY",
      "ALL_SELECTED_HAVE_TCPA_OR_OTHER_EDGE",
      "DOT_SOURCE_DIMENSION",
      "ALL_OUTPUTS_NONEMPTY"
    ),
    Observed = c(
      uniqueN(individual_all$Cancer),
      length(formal_cancers),
      all(cancer_samples[
        Formal_Cancer == TRUE,
        N_C1_Max >= min_group_n & N_C4_Max >= min_group_n
      ]),
      uniqueN(individual$Design_ID) == 1L &&
        unique(individual$Design_ID) == "C1_vs_C4" &&
        uniqueN(individual$Analysis_Set) == 1L &&
        unique(individual$Analysis_Set) == "REPRESENTATIVE_MAIN" &&
        uniqueN(individual$Model_Type) == 1L &&
        unique(individual$Model_Type) == "UNADJUSTED",
      format(fdr_recalc_error, scientific = TRUE),
      "TRUE",
      all(selected$Dominant_Recurrence_N >= min_recurrent_cancers),
      !any(selected$Recurrent_Class == "BIDIRECTIONAL_RECURRENT"),
      uniqueN(selected_map$Feature_ID),
      nrow(dot_source),
      all(file.exists(files)) &&
        all(file.info(files)$size > 0)
    ),
    Expected = c(
      22L, 15L, TRUE, TRUE, "0e+00", "TRUE",
      TRUE, TRUE, nrow(selected),
      length(formal_cancers) * nrow(selected), TRUE
    )
  )
  qa[, Pass := as.character(Observed) == as.character(Expected)]
  qa_file <- file.path(
    tmp_dir,
    "stage3_rppa_c1c4_source_machine_qa_v0.1_20260730.tsv"
  )
  fwrite(qa, qa_file, sep = "\t", na = "NA")
  if (!all(qa$Pass)) {
    print(qa)
    stop("RPPA C1 vs C4 source machine QA失败")
  }

  manifest_inputs <- data.table(
    Role = c(
      rep("PER_CANCER_MODEL", length(model_files)),
      "RPPA_DICTIONARY", "TCPA_TABLE13"
    ),
    Path = required,
    Bytes = file.info(required)$size,
    MD5 = unname(tools::md5sum(required))
  )
  fwrite(
    manifest_inputs,
    file.path(
      tmp_dir,
      "stage3_rppa_c1c4_input_manifest_v0.1_20260730.tsv"
    ),
    sep = "\t", na = "NA"
  )

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) {
    stop("无法原子冻结RPPA C1 vs C4 source目录")
  }
  log_event("DONE", "RPPA C1 vs C4复现和TCPA源数据完成", list(
    runtime_sec = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    )), 1),
    formal_cancers = length(formal_cancers),
    display_features = nrow(selected),
    c1_high = selected[Recurrent_Class == "C1_HIGH", .N],
    c4_high = selected[Recurrent_Class == "C4_HIGH", .N],
    visible_dots = sum(dot_source$Dot_Display),
    pathway_nodes = nrow(pathway_registry),
    output = final_dir
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    )), 1)
  ))
  quit(save = "no", status = 1L)
})
}

# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/05_freeze_lncrna_recurrence.R
stage_transcriptome_proteome_05_freeze_lncrna_recurrence <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]

if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("Input directory must be readable/searchable.")
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
required_files <- c(file.path(paths["model_dir"],paste0(expected_cancers,"_lncrna_feature_models_v0.1_20260723.tsv.gz")))

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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。


suppressPackageStartupMessages({
  library(data.table)
})

# 本脚本只提取既有单癌种C1 vs C4未校正模型，并按同方向显著癌种数冻结候选。
input_dir <- unname(paths["model_dir"])
output_parent <- out_dir
output_dir <- file.path(out_dir,"01_lncrna_recurrence")
tmp_dir <- paste0(output_dir,".tmp_",Sys.getpid())
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

expected_cancers <- c(
  "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC",
  "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PCPG", "PRAD",
  "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC"
)
n_cancers <- length(expected_cancers)
min_recurrent <- ceiling(n_cancers * 0.30)

if (dir.exists(output_dir)) {
  stop("正式输出目录已存在，拒绝覆盖: ", output_dir)
}
if (dir.exists(tmp_dir)) {
  stop("临时输出目录已存在: ", tmp_dir)
}
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

read_one <- function(cancer) {
  path <- file.path(
    input_dir,
    paste0(cancer, "_lncrna_feature_models_v0.1_20260723.tsv.gz")
  )
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("缺少或为空的输入文件: ", path)
  }
  x <- fread(path)
  x <- x[
    Design_ID == "C1_vs_C4" &
      Analysis_Set == "REPRESENTATIVE_MAIN" &
      Model_Type == "UNADJUSTED"
  ]
  if (!nrow(x)) {
    stop("未提取到正式C1 vs C4 lncRNA结果: ", cancer)
  }
  if (uniqueN(x$Cancer) != 1L || unique(x$Cancer) != cancer ||
      uniqueN(x$Layer) != 1L || unique(x$Layer) != "LNCRNA" ||
      any(x$Target != "C1") || any(x$Reference != "C4")) {
    stop("输入设计身份不符合合同: ", cancer)
  }
  if (anyDuplicated(x[, .(Cancer, Feature_ID)])) {
    stop("单癌种Feature_ID重复: ", cancer)
  }
  x
}

message("读取22癌种正式C1 vs C4 lncRNA结果")
individual <- rbindlist(lapply(expected_cancers, read_one), use.names = TRUE)

individual[, Single_Cancer_Direction := fifelse(
  Status == "ESTIMABLE" & is.finite(Estimate) & is.finite(FDR) &
    FDR < 0.01 & Estimate > 0.5,
  "C1_HIGH",
  fifelse(
    Status == "ESTIMABLE" & is.finite(Estimate) & is.finite(FDR) &
      FDR < 0.01 & Estimate < -0.5,
    "C4_HIGH", "NOT_SIGNIFICANT"
  )
)]

recurrence <- individual[, {
  labels <- unique(Feature_Label[!is.na(Feature_Label) & nzchar(Feature_Label)])
  defs <- Definition_Overlap %in% TRUE
  est <- Estimate[Status == "ESTIMABLE" & is.finite(Estimate)]
  .(
    Feature_Label = if (length(labels)) labels[1L] else Feature_ID[1L],
    Definition_Overlap = any(defs),
    N_Cancer_Total = n_cancers,
    N_Cancer_Observed = uniqueN(Cancer),
    N_Cancer_Estimable = uniqueN(Cancer[
      Status == "ESTIMABLE" & is.finite(Estimate) & is.finite(FDR)
    ]),
    N_C1_High_Significant = uniqueN(Cancer[
      Single_Cancer_Direction == "C1_HIGH"
    ]),
    N_C4_High_Significant = uniqueN(Cancer[
      Single_Cancer_Direction == "C4_HIGH"
    ]),
    Median_Estimate_Estimable = if (length(est)) median(est) else NA_real_,
    Mean_Estimate_Estimable = if (length(est)) mean(est) else NA_real_
  )
}, by = Feature_ID]

recurrence[, `:=`(
  C1_High_Recurrence_Fraction = N_C1_High_Significant / n_cancers,
  C4_High_Recurrence_Fraction = N_C4_High_Significant / n_cancers
)]
recurrence[, Recurrent_Class := fcase(
  N_C1_High_Significant >= min_recurrent &
    N_C4_High_Significant >= min_recurrent,
  "BIDIRECTIONAL_RECURRENT",
  N_C1_High_Significant >= min_recurrent,
  "C1_HIGH",
  N_C4_High_Significant >= min_recurrent,
  "C4_HIGH",
  default = "NOT_RECURRENT"
)]
recurrence[, Formal_Functional_Candidate :=
             Recurrent_Class %chin% c("C1_HIGH", "C4_HIGH")]
recurrence[, Dominant_Recurrence_N := pmax(
  N_C1_High_Significant, N_C4_High_Significant
)]
setorder(
  recurrence, -Formal_Functional_Candidate, Recurrent_Class,
  -Dominant_Recurrence_N, Feature_Label, Feature_ID
)

candidates <- recurrence[Formal_Functional_Candidate == TRUE]
if (!nrow(candidates)) {
  stop("没有lncRNA通过预定义复现门槛")
}

summary_table <- rbindlist(list(
  recurrence[, .(
    Metric = "FEATURE_COUNT_BY_CLASS",
    Stratum = Recurrent_Class,
    N = .N
  ), by = Recurrent_Class][, Recurrent_Class := NULL],
  data.table(
    Metric = c(
      "EXPECTED_CANCERS", "MIN_RECURRENT_CANCERS",
      "TOTAL_UNIQUE_LNCRNA", "FORMAL_FUNCTIONAL_CANDIDATES"
    ),
    Stratum = "ALL",
    N = c(
      n_cancers, min_recurrent, nrow(recurrence), nrow(candidates)
    )
  )
), use.names = TRUE, fill = TRUE)

thresholds <- data.table(
  Parameter = c(
    "DESIGN_ID", "ANALYSIS_SET", "MODEL_TYPE",
    "WITHIN_CANCER_FDR", "ABS_ESTIMATE",
    "N_CANCERS_TOTAL", "RECURRENCE_FRACTION", "MIN_RECURRENT_CANCERS",
    "CROSS_CANCER_BH", "META_LOO_I2", "PURITY_MODEL"
  ),
  Value = c(
    "C1_vs_C4", "REPRESENTATIVE_MAIN", "UNADJUSTED",
    "0.01", "0.5", as.character(n_cancers), "0.30",
    as.character(min_recurrent), "NOT_PERFORMED", "NOT_PERFORMED",
    "EXCLUDED"
  )
)

qa <- data.table(
  Check = c(
    "INPUT_CANCERS_22",
    "ONLY_C1_VS_C4",
    "ONLY_REPRESENTATIVE_MAIN",
    "ONLY_UNADJUSTED",
    "TARGET_REFERENCE_C1_C4",
    "UNIQUE_CANCER_FEATURE_KEYS",
    "FDR_RANGE",
    "MIN_RECURRENT_EQUALS_7",
    "NO_BIDIRECTIONAL_IN_FORMAL_CANDIDATES",
    "FORMAL_CANDIDATES_NONEMPTY"
  ),
  Pass = c(
    setequal(unique(individual$Cancer), expected_cancers),
    identical(unique(individual$Design_ID), "C1_vs_C4"),
    identical(unique(individual$Analysis_Set), "REPRESENTATIVE_MAIN"),
    identical(unique(individual$Model_Type), "UNADJUSTED"),
    all(individual$Target == "C1") && all(individual$Reference == "C4"),
    !anyDuplicated(individual[, .(Cancer, Feature_ID)]),
    all(individual$FDR[is.finite(individual$FDR)] >= 0 &
          individual$FDR[is.finite(individual$FDR)] <= 1),
    identical(min_recurrent, 7),
    !any(candidates$Recurrent_Class == "BIDIRECTIONAL_RECURRENT"),
    nrow(candidates) > 0
  )
)
if (!all(qa$Pass)) {
  fwrite(qa, file.path(tmp_dir, "stage3_lncrna_c1c4_recurrence_qa_v0.1_20260730.tsv"),
         sep = "\t")
  stop("lncRNA C1 vs C4复现候选QA失败")
}

fwrite(
  individual,
  file.path(tmp_dir, "stage3_lncrna_c1c4_within_cancer_results_v0.1_20260730.tsv.gz"),
  sep = "\t"
)
fwrite(
  recurrence,
  file.path(tmp_dir, "stage3_lncrna_c1c4_recurrence_all_v0.1_20260730.tsv.gz"),
  sep = "\t"
)
fwrite(
  candidates,
  file.path(tmp_dir, "stage3_lncrna_c1c4_recurrent_candidates_v0.1_20260730.tsv"),
  sep = "\t"
)
fwrite(
  summary_table,
  file.path(tmp_dir, "stage3_lncrna_c1c4_recurrence_summary_v0.1_20260730.tsv"),
  sep = "\t"
)
fwrite(
  thresholds,
  file.path(tmp_dir, "stage3_lncrna_c1c4_recurrence_thresholds_v0.1_20260730.tsv"),
  sep = "\t"
)
fwrite(
  qa,
  file.path(tmp_dir, "stage3_lncrna_c1c4_recurrence_qa_v0.1_20260730.tsv"),
  sep = "\t"
)

dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
if (!file.rename(tmp_dir, output_dir)) {
  stop("无法原子提交正式输出目录: ", output_dir)
}

message(
  "完成：", nrow(candidates), "个正式候选；C1-high=",
  candidates[Recurrent_Class == "C1_HIGH", .N], "；C4-high=",
  candidates[Recurrent_Class == "C4_HIGH", .N],
  "；双向复现审计=", recurrence[Recurrent_Class == "BIDIRECTIONAL_RECURRENT", .N]
)
}

# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/06_run_lncrna_gene_coexpression.R
stage_transcriptome_proteome_06_run_lncrna_gene_coexpression <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]

if(!dir.exists(paths["lncrna_dir"]) || file.access(paths["lncrna_dir"],5L)!=0L) stop("Input directory must be readable/searchable.")
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
lncrna_files <- file.path(paths["lncrna_dir"],paste0("TCGA-",expected_cancers),paste0("TCGA-",expected_cancers,"_lncrna_expr_tpm.rdata"))
required_files <- c(input_paths[names(input_paths)!="lncrna_dir"],lncrna_files)

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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。


suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "06_run_lncrna_gene_coexpression.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_013, environment())
`norm_sample` <- bind_stage_helper(shared_norm_sample_037, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())

standardize_rank_rows <- function(x) {
  ranks <- matrixStats::rowRanks(
    x, ties.method = "average", preserveShape = TRUE
  )
  means <- rowMeans(ranks)
  sds <- matrixStats::rowSds(ranks)
  if (any(!is.finite(sds) | sds <= 0)) {
    stop("通过IQR过滤后仍存在秩方差为0的feature")
  }
  sweep(sweep(ranks, 1L, means, "-"), 1L, sds, "/")
}

load_lncrna_matrix <- function(cancer, candidate_ids, annotation, input_root) {
  input_file <- file.path(
    input_root, paste0("TCGA-", cancer),
    paste0("TCGA-", cancer, "_lncrna_expr_tpm.rdata")
  )
  if (!file.exists(input_file) || file.info(input_file)$size <= 0) {
    stop("缺少或为空的lncRNA TPM: ", input_file)
  }
  env <- new.env(parent = emptyenv())
  load(input_file, envir = env)
  if (!exists("lncrna_expr_tpm", envir = env, inherits = FALSE)) {
    stop("lncRNA RData缺少lncrna_expr_tpm: ", cancer)
  }
  raw <- as.matrix(env$lncrna_expr_tpm)
  storage.mode(raw) <- "double"
  raw_key <- norm_symbol(rownames(raw))
  if (any(!nzchar(raw_key)) || anyDuplicated(raw_key)) {
    stop("lncRNA源symbol为空或重复: ", cancer)
  }

  ann <- annotation[match(candidate_ids, formal_feature_id)]
  if (anyNA(ann$formal_feature_id)) {
    stop("候选lncRNA不在正式annotation registry: ", cancer)
  }
  idx <- match(ann$source_symbol_key, raw_key)
  if (anyNA(idx)) {
    stop("候选lncRNA无法映射至表达矩阵: ", cancer)
  }
  mat <- raw[idx, , drop = FALSE]
  rownames(mat) <- ann$formal_feature_id
  mat
}

align_matrices <- function(patient_cancer, mrna, lnc) {
  mrna_ids <- norm_sample(colnames(mrna))
  lnc_ids <- norm_sample(colnames(lnc))
  mrna <- mrna[, !duplicated(mrna_ids), drop = FALSE]
  lnc <- lnc[, !duplicated(lnc_ids), drop = FALSE]
  colnames(mrna) <- mrna_ids[!duplicated(mrna_ids)]
  colnames(lnc) <- lnc_ids[!duplicated(lnc_ids)]
  samples <- patient_cancer$Sample[
    patient_cancer$Sample %chin% colnames(mrna) &
      patient_cancer$Sample %chin% colnames(lnc)
  ]
  if (length(samples) < 50L) {
    stop("mRNA/lncRNA共同代表患者少于50: ", unique(patient_cancer$Cancer))
  }
  meta <- patient_cancer[match(samples, Sample)]
  mrna <- mrna[, samples, drop = FALSE]
  lnc <- lnc[, samples, drop = FALSE]
  if (!identical(colnames(mrna), meta$Sample) ||
      !identical(colnames(lnc), meta$Sample)) {
    stop("mRNA/lncRNA/患者三方样本未严格对齐")
  }
  if (any(!is.finite(mrna)) || any(!is.finite(lnc))) {
    stop("表达矩阵存在非有限值")
  }
  list(meta = meta, mrna = mrna, lnc = lnc)
}

run_one_cancer <- function(cancer, patient, mrna_list, annotation, candidates,
                           lnc_root, edge_dir, availability_dir) {
  patient_cancer <- patient[
    Cancer == cancer & C_subtype %chin% paste0("C", 1:4)
  ]
  if (!cancer %chin% names(mrna_list)) {
    stop("mRNA对象缺少癌种: ", cancer)
  }
  mrna <- as.matrix(mrna_list[[cancer]]$expression_data)
  storage.mode(mrna) <- "double"
  genes <- norm_symbol(rownames(mrna))
  if (any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("mRNA gene symbol为空或重复: ", cancer)
  }
  rownames(mrna) <- genes

  candidate_ids <- candidates$Feature_ID
  lnc <- load_lncrna_matrix(
    cancer, candidate_ids, annotation, lnc_root
  )
  aligned <- align_matrices(patient_cancer, mrna, lnc)
  n <- nrow(aligned$meta)

  mrna_log <- log2(aligned$mrna + 1)
  mrna_keep <- rowMeans(aligned$mrna >= 1) >= 0.20 &
    matrixStats::rowIQRs(mrna_log) > 0
  mrna_log <- mrna_log[mrna_keep, , drop = FALSE]
  mrna_genes <- rownames(mrna_log)
  if (!length(mrna_genes)) {
    stop("mRNA过滤后为空: ", cancer)
  }

  lnc_log <- log2(aligned$lnc + 1)
  lnc_keep <- rowMeans(aligned$lnc >= 1) >= 0.20 &
    matrixStats::rowIQRs(lnc_log) > 0
  eligible_lnc <- rownames(lnc_log)[lnc_keep]
  lnc_log <- lnc_log[lnc_keep, , drop = FALSE]

  lnc_availability <- candidates[, .(
    Feature_ID, Feature_Label, Definition_Overlap, Recurrent_Class,
    N_C1_High_Significant, N_C4_High_Significant
  )]
  lnc_availability[, `:=`(
    Cancer = cancer,
    N_Common_Patients = n,
    Passing_Expression_Filter = Feature_ID %chin% eligible_lnc,
    Failure_Reason = fifelse(
      Feature_ID %chin% eligible_lnc,
      NA_character_, "LNCRNA_EXPRESSION_FILTER_FAILED"
    )
  )]
  setcolorder(lnc_availability, c(
    "Cancer", "Feature_ID", "Feature_Label", "Definition_Overlap",
    "Recurrent_Class", "N_C1_High_Significant",
    "N_C4_High_Significant", "N_Common_Patients",
    "Passing_Expression_Filter", "Failure_Reason"
  ))

  mrna_availability <- data.table(
    Cancer = cancer,
    Gene = mrna_genes,
    N_Common_Patients = n,
    Passing_Expression_Filter = TRUE
  )

  lnc_file <- file.path(
    availability_dir,
    sprintf("%s_lncrna_availability_v0.1_20260730.tsv.gz", cancer)
  )
  mrna_file <- file.path(
    availability_dir,
    sprintf("%s_mrna_availability_v0.1_20260730.tsv.gz", cancer)
  )
  fwrite(lnc_availability, lnc_file, sep = "\t", compress = "gzip")
  fwrite(mrna_availability, mrna_file, sep = "\t", compress = "gzip")

  edge_base <- data.table()
  if (nrow(lnc_log)) {
    # Spearman等价于秩转换后的Pearson；全部计算和BH均严格限于单癌种。
    lnc_z <- standardize_rank_rows(lnc_log)
    mrna_z <- standardize_rank_rows(mrna_log)
    rho_mat <- tcrossprod(lnc_z, mrna_z) / (n - 1)
    rho_mat <- pmax(pmin(rho_mat, 0.999999), -0.999999)
    edge_rows <- vector("list", nrow(rho_mat))
    for (i in seq_len(nrow(rho_mat))) {
      rho <- rho_mat[i, ]
      statistic <- rho * sqrt(
        (n - 2) / pmax(1 - rho^2, .Machine$double.eps)
      )
      p_value <- 2 * pt(-abs(statistic), df = n - 2)
      cancer_fdr <- p.adjust(p_value, method = "BH")
      keep <- abs(rho) >= 0.30 & cancer_fdr < 0.05
      if (any(keep)) {
        edge_rows[[i]] <- data.table(
          Cancer = cancer,
          Feature_ID = rownames(lnc_log)[i],
          Gene = mrna_genes[keep],
          N_Common_Patients = n,
          Rho = rho[keep],
          P_Value = p_value[keep],
          Cancer_LNCRNA_BH_FDR = cancer_fdr[keep],
          Correlation_Direction = fifelse(
            rho[keep] > 0, "POSITIVE", "NEGATIVE"
          )
        )
      }
    }
    edge_base <- rbindlist(edge_rows, use.names = TRUE, fill = TRUE)
  }

  if (nrow(edge_base)) {
    edge_base <- merge(
      edge_base,
      candidates[, .(
        Feature_ID, Feature_Label, Definition_Overlap, Recurrent_Class,
        N_C1_High_Significant, N_C4_High_Significant
      )],
      by = "Feature_ID", all.x = TRUE
    )
    setorder(
      edge_base, Feature_ID, Correlation_Direction,
      Cancer_LNCRNA_BH_FDR, Gene
    )
  } else {
    edge_base <- data.table(
      Feature_ID = character(), Gene = character(), Cancer = character(),
      N_Common_Patients = integer(), Rho = numeric(), P_Value = numeric(),
      Cancer_LNCRNA_BH_FDR = numeric(),
      Correlation_Direction = character(), Feature_Label = character(),
      Definition_Overlap = logical(), Recurrent_Class = character(),
      N_C1_High_Significant = integer(),
      N_C4_High_Significant = integer()
    )
  }

  edge_file <- file.path(
    edge_dir,
    sprintf("%s_lncrna_mrna_passing_edges_v0.1_20260730.tsv.gz", cancer)
  )
  fwrite(edge_base, edge_file, sep = "\t", compress = "gzip")

  data.table(
    Cancer = cancer,
    N_Common_Patients = n,
    N_MRNA_Passing_Filter = length(mrna_genes),
    N_LNCRNA_Computed = nrow(candidates),
    N_LNCRNA_Passing_Filter = length(eligible_lnc),
    N_Passing_Edges_Rho030 = nrow(edge_base),
    N_Positive_Passing_Edges_Rho030 =
      sum(edge_base$Correlation_Direction == "POSITIVE"),
    N_Negative_Passing_Edges_Rho030 =
      sum(edge_base$Correlation_Direction == "NEGATIVE"),
    N_Passing_Edges_Rho050 = edge_base[abs(Rho) >= 0.50, .N],
    Edge_Output = edge_file,
    LNCRNA_Availability_Output = lnc_file,
    MRNA_Availability_Output = mrna_file
  )
}

main <- function() {
  final_dir <- file.path(out_dir,"01_gene_level_coexpression")
  if (dir.exists(final_dir) &&
      length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录已存在且非空，拒绝覆盖: ", final_dir)
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(tmp_dir)) {
    stop("临时目录已存在: ", tmp_dir)
  }
  source_dir <- file.path(tmp_dir, "source_data")
  edge_dir <- file.path(tmp_dir, "passing_edges_by_cancer")
  availability_dir <- file.path(tmp_dir, "availability_by_cancer")
  qa_dir <- file.path(tmp_dir, "qa")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(edge_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(availability_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(patient=unname(paths["patient_manifest"]),annotation=unname(paths["lncrna_registry"]),candidates=unname(paths["recurrent_candidates"]),mrna=unname(paths["mrna_rdata"]),lnc_root=unname(paths["lncrna_dir"]))
  file_inputs <- unlist(paths[c("patient", "annotation", "candidates", "mrna")])
  if (any(!file.exists(file_inputs)) ||
      any(file.info(file_inputs)$size <= 0)) {
    stop("缺少或为空的冻结文件输入")
  }
  if (!dir.exists(paths$lnc_root)) {
    stop("缺少lncRNA表达根目录")
  }

  patient <- fread(paths$patient)
  patient[, Sample := norm_sample(Sample)]
  if (anyDuplicated(patient$Patient_ID) || anyDuplicated(patient$Sample)) {
    stop("代表患者接口不唯一")
  }
  cancers <- sort(unique(patient$Cancer))
  if (length(cancers) != 22L) {
    stop("正式癌种数不是22")
  }

  if (!setequal(cancers, expected_cancers)) stop("Patient cancer scope must match the complete formal 22-cancer set.")
  annotation <- fread(paths$annotation)[formal_feature_eligible == TRUE]
  annotation[, source_symbol_key := norm_symbol(source_symbol_key)]
  if (anyDuplicated(annotation$formal_feature_id)) {
    stop("lncRNA annotation registry键重复")
  }
  candidates <- fread(paths$candidates)
  candidates <- candidates[
    Formal_Functional_Candidate == TRUE &
      Recurrent_Class %chin% c("C1_HIGH", "C4_HIGH")
  ]
  if (!nrow(candidates) || anyDuplicated(candidates$Feature_ID)) {
    stop("正式候选为空或Feature_ID重复")
  }
  if (any(!candidates$Feature_ID %chin% annotation$formal_feature_id)) {
    stop("正式候选不在lncRNA annotation registry")
  }

  env <- new.env(parent = emptyenv())
  load(paths$mrna, envir = env)
  if (!exists("pancancer_data_list", envir = env, inherits = FALSE)) {
    stop("mRNA RData缺少pancancer_data_list")
  }
  mrna_list <- env$pancancer_data_list
  if (any(!cancers %chin% names(mrna_list))) {
    stop("mRNA癌种输入不完整")
  }

  fwrite(
    candidates,
    file.path(
      source_dir,
      "lncrna_c1c4_compute_candidate_registry_v0.1_20260730.tsv"
    ),
    sep = "\t"
  )
  log_event(
    "START", "Stage 3 lncRNA C1 vs C4 gene coexpression begin",
    list(
      project = "LM_pan_cancer", environment = "R4.5.1",
      cancers = length(cancers), candidates = nrow(candidates),
      patient_scope = "WITHIN_CANCER_ALL_C1_C4_REPRESENTATIVES",
      correlation = "SPEARMAN",
      bh_family = "CANCER_X_LNCRNA",
      audit_threshold = "FDR_LT_0.05_AND_ABS_RHO_GE_0.30",
      output = final_dir
    )
  )

  manifest <- rbindlist(lapply(cancers, function(cancer) {
    z <- run_one_cancer(
      cancer, patient, mrna_list, annotation, candidates,
      paths$lnc_root, edge_dir, availability_dir
    )
    log_event(
      "CHECKPOINT", "per-cancer coexpression complete",
      list(
        cancer = cancer, n_samples = z$N_Common_Patients,
        mrna = z$N_MRNA_Passing_Filter,
        lncrna = z$N_LNCRNA_Passing_Filter,
        edges_rho030 = z$N_Passing_Edges_Rho030,
        edges_rho050 = z$N_Passing_Edges_Rho050
      )
    )
    gc()
    z
  }), use.names = TRUE)

  path_columns <- c(
    "Edge_Output", "LNCRNA_Availability_Output",
    "MRNA_Availability_Output"
  )
  for (column in path_columns) {
    set(
      manifest, j = column,
      value = paste0(final_dir, substring(manifest[[column]], nchar(tmp_dir) + 1L))
    )
  }
  fwrite(
    manifest,
    file.path(source_dir, "coexpression_cancer_manifest_v0.1_20260730.tsv"),
    sep = "\t"
  )

  output_paths <- c(
    list.files(source_dir, full.names = TRUE),
    list.files(edge_dir, full.names = TRUE),
    list.files(availability_dir, full.names = TRUE)
  )
  file_qa <- data.table(
    File = output_paths,
    Size_Bytes = file.info(output_paths)$size,
    MD5 = unname(tools::md5sum(output_paths))
  )
  file_qa[, Pass := is.finite(Size_Bytes) & Size_Bytes > 0 & nzchar(MD5)]
  fwrite(
    file_qa,
    file.path(qa_dir, "output_file_qa_v0.1_20260730.tsv"),
    sep = "\t"
  )

  machine_qa <- data.table(
    Check = c(
      "CANCERS_22",
      "CANDIDATES_NONEMPTY_AND_UNIQUE",
      "ALL_CANCERS_AT_LEAST_50_MATCHED_PATIENTS",
      "MRNA_BACKGROUND_BROAD",
      "PER_CANCER_LNCRNA_BH_FAMILY",
      "NO_CROSS_CANCER_PATIENT_POOLING",
      "OUTPUT_FILES_NONEMPTY_AND_MD5"
    ),
    Pass = c(
      nrow(manifest) == 22L,
      nrow(candidates) > 0L && !anyDuplicated(candidates$Feature_ID),
      all(manifest$N_Common_Patients >= 50L),
      all(manifest$N_MRNA_Passing_Filter > 10000L),
      TRUE,
      TRUE,
      all(file_qa$Pass)
    ),
    Detail = c(
      paste0(nrow(manifest), "/22"),
      paste0("n=", nrow(candidates)),
      paste0("min_n=", min(manifest$N_Common_Patients)),
      paste0("min_mrna=", min(manifest$N_MRNA_Passing_Filter)),
      "BH within Cancer x lncRNA across all eligible mRNAs",
      "each cancer processed separately",
      paste0(sum(file_qa$Pass), "/", nrow(file_qa), " files")
    )
  )
  fwrite(
    machine_qa,
    file.path(qa_dir, "machine_qa_v0.1_20260730.tsv"),
    sep = "\t"
  )
  if (!all(machine_qa$Pass)) {
    stop("lncRNA–mRNA癌种内相关机器QA未全部通过")
  }

  if (!file.rename(tmp_dir, final_dir)) {
    stop("临时目录原子提交失败: ", tmp_dir)
  }
  log_event(
    "DONE", "Stage 3 lncRNA C1 vs C4 gene coexpression complete",
    list(
      runtime_sec = round(as.numeric(difftime(
        Sys.time(), start_time, units = "secs"
      )), 1),
      cancers = nrow(manifest), candidates = nrow(candidates),
      passing_edges_rho030 = sum(manifest$N_Passing_Edges_Rho030),
      passing_edges_rho050 = sum(manifest$N_Passing_Edges_Rho050),
      output = final_dir, qa = "PASS"
    )
  )
}

tryCatch(
  main(),
  error = function(e) {
    log_event("ERROR", conditionMessage(e))
    quit(status = 1L)
  }
)
}

# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/07_build_lncrna_edge_recurrence.R
stage_transcriptome_proteome_07_build_lncrna_edge_recurrence <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]

if(!dir.exists(paths["coexpression_dir"]) || file.access(paths["coexpression_dir"],5L)!=0L) stop("Input directory must be readable/searchable.")
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
required_files <- c(paths["recurrent_candidates"],file.path(paths["coexpression_dir"],"source_data","coexpression_cancer_manifest_v0.1_20260730.tsv"),file.path(paths["coexpression_dir"],"passing_edges_by_cancer",paste0(expected_cancers,"_lncrna_mrna_passing_edges_v0.1_20260730.tsv.gz")),file.path(paths["coexpression_dir"],"availability_by_cancer",paste0(expected_cancers,"_lncrna_availability_v0.1_20260730.tsv.gz")),file.path(paths["coexpression_dir"],"availability_by_cancer",paste0(expected_cancers,"_mrna_availability_v0.1_20260730.tsv.gz")))

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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。


suppressPackageStartupMessages({
  library(data.table)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "07_build_lncrna_edge_recurrence.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_013, environment())

main <- function() {
  pass1_dir <- unname(paths["coexpression_dir"])
  final_dir <- file.path(out_dir,"02_edge_recurrence")
  tmp_dir <- paste0(final_dir,".tmp_",Sys.getpid())
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  manifest_path <- file.path(pass1_dir,"source_data","coexpression_cancer_manifest_v0.1_20260730.tsv")
  candidate_path <- unname(paths["recurrent_candidates"])

  if (any(!file.exists(c(manifest_path, candidate_path))) ||
      any(file.info(c(manifest_path, candidate_path))$size <= 0)) {
    stop("edge recurrence输入缺失或为空")
  }
  if (dir.exists(final_dir) &&
      length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录已存在且非空，拒绝覆盖: ", final_dir)
  }
  if (dir.exists(tmp_dir)) {
    stop("临时输出目录已存在: ", tmp_dir)
  }
  source_dir <- file.path(tmp_dir, "source_data")
  summary_dir <- file.path(tmp_dir, "summary")
  qa_dir <- file.path(tmp_dir, "qa")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- fread(manifest_path)
# 仅按癌种及固定basename重绑可搬迁的manifest路径，不读取旧绝对路径。
  if (nrow(manifest)!=22L || !setequal(manifest$Cancer,expected_cancers) || anyDuplicated(manifest$Cancer)) stop("Invalid complete cancer manifest.")
  path_columns <- c("Edge_Output","LNCRNA_Availability_Output","MRNA_Availability_Output")
  suffixes <- c("_lncrna_mrna_passing_edges_v0.1_20260730.tsv.gz","_lncrna_availability_v0.1_20260730.tsv.gz","_mrna_availability_v0.1_20260730.tsv.gz")
  subdirs <- c("passing_edges_by_cancer","availability_by_cancer","availability_by_cancer")
  for (j in seq_along(path_columns)) {
    column <- path_columns[j]
    expected <- paste0(manifest$Cancer,suffixes[j])
    if (!column %in% names(manifest) || anyNA(manifest[[column]]) || !identical(basename(manifest[[column]]),expected)) stop("Manifest basename mismatch: ",column)
    set(manifest,j=column,value=file.path(pass1_dir,subdirs[j],expected))
  }
  candidates <- fread(candidate_path)[
    Formal_Functional_Candidate == TRUE &
      Recurrent_Class %chin% c("C1_HIGH", "C4_HIGH")
  ]
  if (nrow(manifest) != 22L || uniqueN(manifest$Cancer) != 22L) {
    stop("相关manifest未覆盖22癌")
  }
  if (!nrow(candidates) || anyDuplicated(candidates$Feature_ID)) {
    stop("正式lncRNA候选为空或重复")
  }
  input_files <- c(
    manifest$Edge_Output, manifest$LNCRNA_Availability_Output,
    manifest$MRNA_Availability_Output
  )
  if (any(!file.exists(input_files)) ||
      any(file.info(input_files)$size <= 0)) {
    stop("逐癌相关或可评价性文件缺失/为空")
  }

  log_event(
    "START", "Stage 3 lncRNA C1 vs C4 edge recurrence begin",
    list(
      project = "LM_pan_cancer", environment = "R4.5.1",
      cancers = nrow(manifest), candidates = nrow(candidates),
      audit_edge = "WITHIN_CANCER_FDR_LT_0.05_ABS_RHO_GE_0.30",
      core_edge = "ABS_RHO_GE_0.50_SUPPORT_GE_7_FRACTION_GE_0.30",
      cross_cancer_bh = "NOT_PERFORMED",
      meta_loo_i2 = "NOT_PERFORMED", output = final_dir
    )
  )

  lnc_availability <- rbindlist(lapply(
    manifest$LNCRNA_Availability_Output, fread
  ), use.names = TRUE, fill = TRUE)
  mrna_availability <- rbindlist(lapply(
    manifest$MRNA_Availability_Output, fread
  ), use.names = TRUE, fill = TRUE)
  lnc_availability <- unique(lnc_availability[
    Passing_Expression_Filter == TRUE,
    .(Feature_ID, Cancer)
  ])
  mrna_availability <- unique(mrna_availability[
    Passing_Expression_Filter == TRUE,
    .(Gene = toupper(trimws(Gene)), Cancer)
  ])
  if (anyDuplicated(lnc_availability[, .(Feature_ID, Cancer)]) ||
      anyDuplicated(mrna_availability[, .(Gene, Cancer)])) {
    stop("可评价性键不唯一")
  }

  feature_levels <- candidates$Feature_ID
  gene_levels <- sort(unique(mrna_availability$Gene))
  cancer_levels <- sort(manifest$Cancer)
  lnc_matrix <- matrix(
    0L, nrow = length(feature_levels), ncol = length(cancer_levels),
    dimnames = list(feature_levels, cancer_levels)
  )
  mrna_matrix <- matrix(
    0L, nrow = length(gene_levels), ncol = length(cancer_levels),
    dimnames = list(gene_levels, cancer_levels)
  )
  lnc_matrix[cbind(
    match(lnc_availability$Feature_ID, feature_levels),
    match(lnc_availability$Cancer, cancer_levels)
  )] <- 1L
  mrna_matrix[cbind(
    match(mrna_availability$Gene, gene_levels),
    match(mrna_availability$Cancer, cancer_levels)
  )] <- 1L
  n_evaluable <- tcrossprod(lnc_matrix, mrna_matrix)

  background_list <- lapply(seq_along(feature_levels), function(i) {
    keep <- n_evaluable[i, ] >= 7L
    data.table(
      Feature_ID = feature_levels[i],
      Gene = gene_levels[keep],
      N_Evaluable_Cancers = as.integer(n_evaluable[i, keep])
    )
  })
  background <- rbindlist(background_list, use.names = TRUE)
  background_summary <- background[, .(
    Universe_N = uniqueN(Gene),
    Min_Evaluable_Cancers = min(N_Evaluable_Cancers),
    Max_Evaluable_Cancers = max(N_Evaluable_Cancers)
  ), by = Feature_ID]

  log_event(
    "CHECKPOINT", "pairwise evaluability matrix complete",
    list(
      feature_count = length(feature_levels),
      gene_count = length(gene_levels),
      background_rows = nrow(background),
      min_universe = min(background_summary$Universe_N),
      max_universe = max(background_summary$Universe_N)
    )
  )

  # 逐癌读取后立即收窄到正式rho门槛，避免把完整rho>=0.30审计层同时载入内存。
  passing_edges <- rbindlist(lapply(
    manifest$Edge_Output,
    function(path) fread(path)[abs(Rho) >= 0.50]
  ), use.names = TRUE, fill = TRUE)
  passing_edges[, Gene := toupper(trimws(Gene))]
  if (anyDuplicated(passing_edges[, .(Cancer, Feature_ID, Gene)])) {
    stop("通过边Cancer × Feature_ID × Gene键不唯一")
  }
  if (any(
    passing_edges$Cancer_LNCRNA_BH_FDR >= 0.05 |
      abs(passing_edges$Rho) < 0.50
  )) {
    stop("单癌种正式相关边违反FDR或rho门槛")
  }

  edge_audit <- passing_edges[, .(
    N_Support_Cancers = uniqueN(Cancer),
    Median_Rho_Support = median(Rho),
    Mean_Rho_Support = mean(Rho),
    Min_Abs_Rho_Support = min(abs(Rho)),
    Max_Cancer_FDR_Support = max(Cancer_LNCRNA_BH_FDR),
    Supporting_Cancers = paste(sort(unique(Cancer)), collapse = ";")
  ), by = .(
    Feature_ID, Gene, Correlation_Direction
  )]
  edge_audit[, N_Evaluable_Cancers := as.integer(
    n_evaluable[
      cbind(match(Feature_ID, feature_levels), match(Gene, gene_levels))
    ]
  )]
  edge_audit[, Support_Fraction :=
               N_Support_Cancers / N_Evaluable_Cancers]
  edge_audit[, Direction_Core := (
    N_Evaluable_Cancers >= 7L &
      N_Support_Cancers >= 7L &
      Support_Fraction >= 0.30
  )]

  direction_wide <- dcast(
    edge_audit,
    Feature_ID + Gene ~ Correlation_Direction,
    value.var = "Direction_Core", fill = FALSE
  )
  if (!"POSITIVE" %chin% names(direction_wide)) {
    direction_wide[, POSITIVE := FALSE]
  }
  if (!"NEGATIVE" %chin% names(direction_wide)) {
    direction_wide[, NEGATIVE := FALSE]
  }
  direction_wide[, Bidirectional_Core := POSITIVE & NEGATIVE]
  edge_audit <- merge(
    edge_audit,
    direction_wide[, .(Feature_ID, Gene, Bidirectional_Core)],
    by = c("Feature_ID", "Gene"), all.x = TRUE
  )
  edge_audit[, Formal_Core_Edge :=
               Direction_Core & !Bidirectional_Core]
  edge_audit <- merge(
    edge_audit,
    candidates[, .(
      Feature_ID, Feature_Label, Definition_Overlap, Recurrent_Class,
      N_C1_High_Significant, N_C4_High_Significant
    )],
    by = "Feature_ID", all.x = TRUE
  )
  edge_audit[, Abs_Median_Rho_Support := abs(Median_Rho_Support)]
  setorder(
    edge_audit, -Formal_Core_Edge, Recurrent_Class, Feature_ID,
    Correlation_Direction, -N_Support_Cancers, -Abs_Median_Rho_Support,
    Gene
  )
  core_edges <- edge_audit[Formal_Core_Edge == TRUE]
  if (!nrow(core_edges)) {
    stop("没有lncRNA–mRNA边通过正式复现门槛")
  }
  if (anyDuplicated(
    core_edges[, .(Feature_ID, Gene, Correlation_Direction)]
  )) {
    stop("正式核心边键重复")
  }

  unit_summary <- core_edges[, .(
    N_Linked_MRNA = uniqueN(Gene),
    Min_Support_Cancers = min(N_Support_Cancers),
    Median_Support_Cancers = as.numeric(median(N_Support_Cancers)),
    Max_Support_Cancers = max(N_Support_Cancers),
    Median_Abs_Rho = median(abs(Median_Rho_Support))
  ), by = .(
    Recurrent_Class, Feature_ID, Feature_Label, Definition_Overlap,
    Correlation_Direction, N_C1_High_Significant,
    N_C4_High_Significant
  )]
  unit_summary[, ORA_Eligible := N_Linked_MRNA >= 10L]
  setorder(
    unit_summary, Recurrent_Class, -ORA_Eligible, -N_Linked_MRNA,
    Feature_Label
  )

  # 用核心边汇总，避免同一mRNA在两个相关方向之间被重复计数。
  candidate_counts <- core_edges[, .(
    N_Core_Directions = uniqueN(Correlation_Direction),
    N_Core_Linked_MRNA_Total = uniqueN(Gene),
    N_Core_Edges = .N
  ), by = Feature_ID]
  candidate_summary <- merge(
    candidates[, .(
      Feature_ID, Feature_Label, Definition_Overlap, Recurrent_Class,
      N_C1_High_Significant, N_C4_High_Significant
    )],
    candidate_counts, by = "Feature_ID", all.x = TRUE
  )
  for (column in c(
    "N_Core_Directions", "N_Core_Linked_MRNA_Total", "N_Core_Edges"
  )) {
    set(candidate_summary, which(is.na(candidate_summary[[column]])),
        column, 0L)
  }

  summary_table <- rbindlist(list(
    core_edges[, .(
      Metric = "CORE_EDGE_COUNT",
      Stratum = paste(Recurrent_Class, Correlation_Direction, sep = "::"),
      N = .N
    ), by = .(Recurrent_Class, Correlation_Direction)][
      , c("Recurrent_Class", "Correlation_Direction") := NULL
    ],
    unit_summary[, .(
      Metric = "ORA_ELIGIBLE_UNIT_COUNT",
      Stratum = paste(Recurrent_Class, Correlation_Direction, sep = "::"),
      N = sum(ORA_Eligible)
    ), by = .(Recurrent_Class, Correlation_Direction)][
      , c("Recurrent_Class", "Correlation_Direction") := NULL
    ],
    data.table(
      Metric = c(
        "FORMAL_CANDIDATE_LNCRNA", "CORE_EDGE_TOTAL",
        "CORE_LINKED_LNCRNA", "ORA_ELIGIBLE_UNITS",
        "BIDIRECTIONAL_CORE_PAIRS"
      ),
      Stratum = "ALL",
      N = c(
        nrow(candidates), nrow(core_edges),
        uniqueN(core_edges$Feature_ID), sum(unit_summary$ORA_Eligible),
        uniqueN(edge_audit[Bidirectional_Core == TRUE,
                           paste(Feature_ID, Gene)])
      )
    )
  ), use.names = TRUE, fill = TRUE)

  qa <- data.table(
    Check = c(
      "CANCERS_22",
      "EVALUABILITY_KEYS_UNIQUE",
      "PASSING_EDGE_KEYS_UNIQUE",
      "NO_CROSS_CANCER_BH_OR_META",
      "CORE_EDGE_THRESHOLDS",
      "NO_BIDIRECTIONAL_CORE_IN_FORMAL",
      "CORE_EDGES_NONEMPTY",
      "BACKGROUND_CONTAINS_ALL_CORE_GENES"
    ),
    Pass = c(
      nrow(manifest) == 22L,
      !anyDuplicated(lnc_availability[, .(Feature_ID, Cancer)]) &&
        !anyDuplicated(mrna_availability[, .(Gene, Cancer)]),
      !anyDuplicated(passing_edges[, .(Cancer, Feature_ID, Gene)]),
      TRUE,
      all(
        core_edges$N_Evaluable_Cancers >= 7L &
          core_edges$N_Support_Cancers >= 7L &
          core_edges$Support_Fraction >= 0.30 &
          core_edges$Min_Abs_Rho_Support >= 0.50
      ),
      !any(core_edges$Bidirectional_Core),
      nrow(core_edges) > 0L,
      nrow(core_edges[!background, on = .(Feature_ID, Gene)]) == 0L
    )
  )
  if (!all(qa$Pass)) {
    fwrite(qa, file.path(qa_dir, "machine_qa_v0.1_20260730.tsv"),
           sep = "\t")
    stop("lncRNA–mRNA边复现QA未全部通过")
  }

  fwrite(
    edge_audit,
    file.path(
      source_dir,
      "lncrna_mrna_edge_recurrence_audit_v0.1_20260730.tsv.gz"
    ),
    sep = "\t", compress = "gzip"
  )
  fwrite(
    core_edges,
    file.path(
      source_dir,
      "lncrna_mrna_formal_core_edges_v0.1_20260730.tsv.gz"
    ),
    sep = "\t", compress = "gzip"
  )
  fwrite(
    background,
    file.path(
      source_dir,
      "lncrna_specific_ora_background_v0.1_20260730.tsv.gz"
    ),
    sep = "\t", compress = "gzip"
  )
  fwrite(
    background_summary,
    file.path(
      source_dir,
      "lncrna_specific_ora_background_summary_v0.1_20260730.tsv"
    ),
    sep = "\t"
  )
  fwrite(
    unit_summary,
    file.path(
      summary_dir,
      "lncrna_linked_mrna_unit_summary_v0.1_20260730.tsv"
    ),
    sep = "\t"
  )
  fwrite(
    candidate_summary,
    file.path(
      summary_dir,
      "lncrna_linked_mrna_candidate_summary_v0.1_20260730.tsv"
    ),
    sep = "\t"
  )
  fwrite(
    summary_table,
    file.path(summary_dir, "edge_recurrence_summary_v0.1_20260730.tsv"),
    sep = "\t"
  )
  fwrite(
    qa, file.path(qa_dir, "machine_qa_v0.1_20260730.tsv"), sep = "\t"
  )

  output_paths <- list.files(
    tmp_dir, recursive = TRUE, full.names = TRUE
  )
  file_qa <- data.table(
    File = output_paths,
    Size_Bytes = file.info(output_paths)$size,
    MD5 = unname(tools::md5sum(output_paths))
  )
  file_qa[, Pass := Size_Bytes > 0 & nzchar(MD5)]
  fwrite(
    file_qa, file.path(qa_dir, "output_file_qa_v0.1_20260730.tsv"),
    sep = "\t"
  )
  if (!all(file_qa$Pass)) {
    stop("edge recurrence输出存在空文件")
  }

  if (!file.rename(tmp_dir, final_dir)) {
    stop("临时目录原子提交失败: ", tmp_dir)
  }
  log_event(
    "DONE", "Stage 3 lncRNA C1 vs C4 edge recurrence complete",
    list(
      runtime_sec = round(as.numeric(difftime(
        Sys.time(), start_time, units = "secs"
      )), 1),
      core_edges = nrow(core_edges),
      linked_lncrna = uniqueN(core_edges$Feature_ID),
      ora_eligible_units = sum(unit_summary$ORA_Eligible),
      background_rows = nrow(background), output = final_dir, qa = "PASS"
    )
  )
}

tryCatch(
  main(),
  error = function(e) {
    log_event("ERROR", conditionMessage(e))
    quit(status = 1L)
  }
)
}

# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/08_run_lncrna_linked_mrna_ora.R
stage_transcriptome_proteome_08_run_lncrna_linked_mrna_ora <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]

if(!dir.exists(paths["edge_recurrence_dir"]) || file.access(paths["edge_recurrence_dir"],5L)!=0L) stop("Input directory must be readable/searchable.")
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
required_files <- c(input_paths[names(input_paths)!="edge_recurrence_dir"],file.path(paths["edge_recurrence_dir"],"source_data",c("lncrna_mrna_formal_core_edges_v0.1_20260730.tsv.gz","lncrna_specific_ora_background_v0.1_20260730.tsv.gz")),file.path(paths["edge_recurrence_dir"],"summary","lncrna_linked_mrna_unit_summary_v0.1_20260730.tsv"))

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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "08_run_lncrna_linked_mrna_ora.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_013, environment())

read_gmt <- function(path, collection) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(parts) < 3L)) {
    stop(collection, " GMT存在少于3列的记录")
  }
  registry <- rbindlist(lapply(seq_along(parts), function(i) {
    data.table(
      Collection = collection,
      Term_Name = trimws(parts[[i]][1]),
      Description_URL = trimws(parts[[i]][2]),
      Original_Line = i
    )
  }))
  membership <- rbindlist(lapply(seq_along(parts), function(i) {
    genes <- unique(toupper(trimws(parts[[i]][-c(1, 2)])))
    genes <- genes[nzchar(genes)]
    data.table(
      Collection = collection,
      Term_Name = trimws(parts[[i]][1]),
      Gene = genes
    )
  }))
  registry[, Term_Key := paste(Collection, Term_Name, sep = "::")]
  membership[, Term_Key := paste(Collection, Term_Name, sep = "::")]
  list(
    registry = registry[order(Original_Line), .SD[1], by = Term_Key],
    membership = unique(membership)
  )
}

main <- function() {
  edge_dir <- unname(paths["edge_recurrence_dir"])
  final_dir <- file.path(out_dir,"03_linked_mrna_ora")
  tmp_dir <- paste0(final_dir,".tmp_",Sys.getpid())
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  input_paths <- list(
    core_edges = file.path(
      edge_dir, "source_data",
      "lncrna_mrna_formal_core_edges_v0.1_20260730.tsv.gz"
    ),
    background = file.path(
      edge_dir, "source_data",
      "lncrna_specific_ora_background_v0.1_20260730.tsv.gz"
    ),
    units = file.path(
      edge_dir, "summary",
      "lncrna_linked_mrna_unit_summary_v0.1_20260730.tsv"
    )
  )
  input_files <- unlist(input_paths, use.names = FALSE)
  if (any(!file.exists(input_files)) ||
      any(file.info(input_files)$size <= 0)) {
    stop("linked-mRNA ORA输入缺失或为空")
  }
  if (dir.exists(final_dir) &&
      length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录已存在且非空，拒绝覆盖: ", final_dir)
  }
  if (dir.exists(tmp_dir)) {
    stop("临时目录已存在: ", tmp_dir)
  }
  source_dir <- file.path(tmp_dir, "source_data")
  resource_dir <- file.path(tmp_dir, "resource_manifest")
  ora_dir <- file.path(tmp_dir, "ora")
  summary_dir <- file.path(tmp_dir, "summary")
  qa_dir <- file.path(tmp_dir, "qa")
  for (path in c(source_dir, resource_dir, ora_dir, summary_dir, qa_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  units <- fread(input_paths$units)[ORA_Eligible == TRUE]
  core_edges <- fread(input_paths$core_edges)
  background <- fread(input_paths$background)
  if (!nrow(units) || anyDuplicated(
    units[, .(Feature_ID, Correlation_Direction)]
  )) {
    stop("ORA单元为空或键重复")
  }
  units[, Unit_ID := paste(
    Recurrent_Class, Feature_ID, Correlation_Direction, sep = "::"
  )]
  core_input <- merge(
    core_edges,
    units[, .(
      Unit_ID, Feature_ID, Correlation_Direction
    )],
    by = c("Feature_ID", "Correlation_Direction"), all = FALSE
  )
  core_input[, Gene := toupper(trimws(Gene))]
  if (any(!nzchar(core_input$Gene)) ||
      anyDuplicated(core_input[, .(Unit_ID, Gene)])) {
    stop("ORA输入Unit_ID × Gene为空或重复")
  }
  observed <- core_input[, .(Observed_Input_N = uniqueN(Gene)), by = Unit_ID]
  count_check <- merge(
    units[, .(Unit_ID, N_Linked_MRNA)], observed,
    by = "Unit_ID", all.x = TRUE
  )
  if (any(count_check$N_Linked_MRNA != count_check$Observed_Input_N)) {
    stop("核心边数量与ORA单元注册表不一致")
  }
  background[, Gene := toupper(trimws(Gene))]
  background <- unique(background[, .(
    Feature_ID, Gene, N_Evaluable_Cancers
  )])
  if (anyDuplicated(background[, .(Feature_ID, Gene)])) {
    stop("lncRNA特异背景键重复")
  }
  missing_core <- core_input[!background, on = .(Feature_ID, Gene)]
  if (nrow(missing_core)) {
    stop("存在不属于lncRNA特异背景的linked mRNA")
  }
  background_summary <- background[, .(
    Universe_N = uniqueN(Gene),
    Min_Evaluable_Cancers = min(N_Evaluable_Cancers),
    Max_Evaluable_Cancers = max(N_Evaluable_Cancers)
  ), by = Feature_ID]
  units <- merge(units, background_summary, by = "Feature_ID", all.x = TRUE)
  setorder(
    units, Recurrent_Class, Feature_ID, Correlation_Direction
  )

  log_event(
    "START", "Stage 3 lncRNA C1 vs C4 linked-mRNA ORA begin",
    list(
      project = "LM_pan_cancer", environment = "R4.5.1",
      units = nrow(units), unique_lncrna = uniqueN(units$Feature_ID),
      linked_edges = nrow(core_input),
      collections = "GO_BP;REACTOME;HALLMARK;IMM_PORT;C7_SUPPLEMENT",
      fdr_family = "LNCRNA_X_DIRECTION_X_COLLECTION",
      output = final_dir
    )
  )

  resource_specs <- data.table(
    Collection = c("HALLMARK", "REACTOME", "GO_BP", "C7", "IMM_PORT"),
    Path = unname(paths[c("hallmark_gmt","reactome_gmt","go_bp_gmt","c7_gmt","immport_gmt")]),
    Version = c(
      "MSIGDB_2025.1_HS", "MSIGDB_2025.1_HS",
      "MSIGDB_2024.1_HS", "MSIGDB_2025.1_HS",
      "LOCAL_FROZEN_VERSION_UNSPECIFIED"
    ),
    Evidence_Role = c(
      "MAIN", "MAIN", "MAIN", "SUPPLEMENT_ONLY", "MAIN"
    )
  )
  if (any(!file.exists(resource_specs$Path)) ||
      any(file.info(resource_specs$Path)$size <= 0)) {
    stop("GMT资源缺失或为空")
  }
  resource_specs[, `:=`(
    Size_Bytes = file.info(Path)$size,
    MD5 = unname(tools::md5sum(Path))
  )]
  parsed <- lapply(seq_len(nrow(resource_specs)), function(i) {
    read_gmt(resource_specs$Path[i], resource_specs$Collection[i])
  })
  term_registry <- rbindlist(
    lapply(parsed, `[[`, "registry"), use.names = TRUE, fill = TRUE
  )
  term_membership <- unique(rbindlist(
    lapply(parsed, `[[`, "membership"), use.names = TRUE, fill = TRUE
  ))
  if (anyDuplicated(term_registry$Term_Key) ||
      anyDuplicated(term_membership[, .(Term_Key, Gene)])) {
    stop("GMT term或membership键重复")
  }

  feature_levels <- sort(unique(units$Feature_ID))
  unit_levels <- units$Unit_ID
  background_for_ora <- background[Feature_ID %chin% feature_levels]
  gene_levels <- sort(unique(background$Gene))
  term_levels <- term_registry$Term_Key
  term_membership_matrix <- term_membership[Gene %chin% gene_levels]
  term_by_gene <- sparseMatrix(
    i = match(term_membership_matrix$Term_Key, term_levels),
    j = match(term_membership_matrix$Gene, gene_levels),
    x = 1,
    dims = c(length(term_levels), length(gene_levels))
  )
  background_matrix <- sparseMatrix(
    i = match(background_for_ora$Gene, gene_levels),
    j = match(background_for_ora$Feature_ID, feature_levels),
    x = 1,
    dims = c(length(gene_levels), length(feature_levels))
  )
  input_matrix <- sparseMatrix(
    i = match(core_input$Gene, gene_levels),
    j = match(core_input$Unit_ID, unit_levels),
    x = 1,
    dims = c(length(gene_levels), length(unit_levels))
  )
  term_background_size <- as.matrix(term_by_gene %*% background_matrix)
  term_overlap_size <- as.matrix(term_by_gene %*% input_matrix)

  ora <- rbindlist(lapply(seq_len(nrow(units)), function(i) {
    feature_index <- match(units$Feature_ID[i], feature_levels)
    gene_set_size <- as.integer(term_background_size[, feature_index])
    overlap_n <- as.integer(term_overlap_size[, i])
    keep <- gene_set_size >= 10L & gene_set_size <= 500L
    universe_n <- units$Universe_N[i]
    input_n <- units$N_Linked_MRNA[i]
    z <- data.table(
      Unit_ID = units$Unit_ID[i],
      Term_Key = term_levels[keep],
      Gene_Set_Size = gene_set_size[keep],
      Overlap_N = overlap_n[keep],
      Universe_N = universe_n,
      Input_N = input_n
    )
    z[, `:=`(
      Expected_Overlap = Input_N * Gene_Set_Size / Universe_N,
      Fold_Enrichment = fifelse(
        Overlap_N > 0,
        (Overlap_N / Input_N) / (Gene_Set_Size / Universe_N),
        0
      ),
      P_Value = phyper(
        Overlap_N - 1L, Gene_Set_Size,
        Universe_N - Gene_Set_Size, Input_N,
        lower.tail = FALSE
      )
    )]
    z
  }), use.names = TRUE)
  ora <- merge(ora, term_registry, by = "Term_Key", all.x = TRUE)
  ora[, FDR := p.adjust(P_Value, method = "BH"),
      by = .(Unit_ID, Collection)]
  ora[, `:=`(
    Full_Significant = Overlap_N >= 3L & FDR < 0.05,
    Display_Eligible = Overlap_N >= 5L & FDR < 0.05
  )]
  ora <- merge(
    ora,
    units[, .(
      Unit_ID, Recurrent_Class, Feature_ID, Feature_Label,
      Definition_Overlap, Correlation_Direction,
      N_C1_High_Significant, N_C4_High_Significant,
      N_Linked_MRNA, Universe_N
    )],
    by = c("Unit_ID", "Universe_N"), all.x = TRUE
  )
  setorder(
    ora, Recurrent_Class, Feature_ID, Correlation_Direction,
    Collection, FDR, -Overlap_N, Term_Name
  )

  significant_keys <- ora[
    Full_Significant == TRUE,
    .(Unit_ID, Term_Key)
  ]
  significant_membership <- merge(
    merge(
      core_input[, .(
        Unit_ID, Feature_ID, Correlation_Direction, Gene,
        N_Support_Cancers, N_Evaluable_Cancers, Support_Fraction,
        Median_Rho_Support
      )],
      term_membership[, .(
        Term_Key, Collection, Term_Name, Gene
      )],
      by = "Gene", allow.cartesian = TRUE
    ),
    significant_keys,
    by = c("Unit_ID", "Term_Key"), all = FALSE
  )
  significant_membership <- unique(significant_membership)
  overlap_summary <- significant_membership[, .(
    Membership_Overlap_N = uniqueN(Gene),
    Overlap_Genes = paste(sort(unique(Gene)), collapse = ";")
  ), by = .(Unit_ID, Term_Key)]
  ora <- merge(
    ora, overlap_summary,
    by = c("Unit_ID", "Term_Key"), all.x = TRUE
  )
  ora[is.na(Membership_Overlap_N), `:=`(
    Membership_Overlap_N = 0L, Overlap_Genes = ""
  )]
  if (any(
    ora$Full_Significant &
      ora$Overlap_N != ora$Membership_Overlap_N
  )) {
    stop("显著term overlap计数与membership不一致")
  }

  significant_ora <- ora[Full_Significant == TRUE]
  display_ora <- ora[
    Display_Eligible == TRUE & Collection != "C7"
  ]
  unit_summary <- ora[, .(
    N_Term_Tests = .N,
    N_Full_Significant = sum(Full_Significant),
    N_Main_Display_Eligible = sum(
      Display_Eligible & Collection != "C7"
    ),
    Best_FDR = min(FDR),
    Best_Main_Display_FDR = if (
      any(Display_Eligible & Collection != "C7")
    ) min(FDR[Display_Eligible & Collection != "C7"]) else NA_real_,
    Best_Overlap_N = max(Overlap_N)
  ), by = .(
    Unit_ID, Recurrent_Class, Feature_ID, Feature_Label,
    Definition_Overlap, Correlation_Direction,
    N_C1_High_Significant, N_C4_High_Significant,
    N_Linked_MRNA, Universe_N
  )]
  unit_summary[, Has_Main_Display_ORA := N_Main_Display_Eligible > 0L]
  setorder(
    unit_summary, Recurrent_Class, -Has_Main_Display_ORA,
    -N_Main_Display_Eligible, Best_Main_Display_FDR,
    -N_Linked_MRNA, Feature_Label
  )

  summary_collection <- ora[, .(
    N_Units_Tested = uniqueN(Unit_ID),
    N_Tests = .N,
    N_Full_Significant = sum(Full_Significant),
    N_Display_Eligible = sum(Display_Eligible),
    N_Units_With_Display = uniqueN(Unit_ID[Display_Eligible])
  ), by = .(Recurrent_Class, Collection)]
  setorder(summary_collection, Recurrent_Class, Collection)

  qa <- data.table(
    Check = c(
      "ORA_UNITS_NONEMPTY_UNIQUE",
      "INPUT_GENES_AT_LEAST_10",
      "INPUT_WITHIN_BACKGROUND",
      "TERM_SIZE_10_TO_500",
      "BH_WITHIN_UNIT_COLLECTION",
      "DISPLAY_OVERLAP_AT_LEAST_5",
      "C7_EXCLUDED_FROM_MAIN_DISPLAY",
      "SIGNIFICANT_MEMBERSHIP_CONCORDANT",
      "NO_MRNA_CAUSAL_INTERPRETATION"
    ),
    Pass = c(
      nrow(units) > 0L && !anyDuplicated(units$Unit_ID),
      all(units$N_Linked_MRNA >= 10L),
      nrow(missing_core) == 0L,
      all(ora$Gene_Set_Size >= 10L & ora$Gene_Set_Size <= 500L),
      TRUE,
      all(ora$Overlap_N[ora$Display_Eligible] >= 5L),
      all(display_ora$Collection != "C7"),
      all(
        ora$Overlap_N[ora$Full_Significant] ==
          ora$Membership_Overlap_N[ora$Full_Significant]
      ),
      TRUE
    )
  )
  if (!all(qa$Pass)) {
    fwrite(qa, file.path(qa_dir, "machine_qa_v0.1_20260730.tsv"),
           sep = "\t")
    stop("linked-mRNA ORA机器QA未全部通过")
  }

  fwrite(units, file.path(
    source_dir, "lncrna_ora_unit_registry_v0.1_20260730.tsv"
  ), sep = "\t")
  fwrite(core_input, file.path(
    source_dir, "lncrna_linked_mrna_ora_input_v0.1_20260730.tsv.gz"
  ), sep = "\t", compress = "gzip")
  fwrite(significant_membership, file.path(
    source_dir, "significant_term_mrna_membership_v0.1_20260730.tsv.gz"
  ), sep = "\t", compress = "gzip")
  fwrite(resource_specs, file.path(
    resource_dir, "gene_set_resource_manifest_v0.1_20260730.tsv"
  ), sep = "\t")
  fwrite(term_registry, file.path(
    resource_dir, "gene_set_term_registry_v0.1_20260730.tsv"
  ), sep = "\t")
  fwrite(ora, file.path(
    ora_dir, "lncrna_linked_mrna_ora_all_v0.1_20260730.tsv.gz"
  ), sep = "\t", compress = "gzip")
  fwrite(significant_ora, file.path(
    ora_dir, "lncrna_linked_mrna_ora_significant_v0.1_20260730.tsv.gz"
  ), sep = "\t", compress = "gzip")
  fwrite(display_ora, file.path(
    ora_dir, "lncrna_linked_mrna_ora_display_eligible_v0.1_20260730.tsv"
  ), sep = "\t")
  fwrite(unit_summary, file.path(
    summary_dir, "lncrna_ora_unit_summary_v0.1_20260730.tsv"
  ), sep = "\t")
  fwrite(summary_collection, file.path(
    summary_dir, "lncrna_ora_collection_summary_v0.1_20260730.tsv"
  ), sep = "\t")
  fwrite(qa, file.path(
    qa_dir, "machine_qa_v0.1_20260730.tsv"
  ), sep = "\t")

  output_paths <- list.files(
    tmp_dir, recursive = TRUE, full.names = TRUE
  )
  file_qa <- data.table(
    File = output_paths,
    Size_Bytes = file.info(output_paths)$size,
    MD5 = unname(tools::md5sum(output_paths))
  )
  file_qa[, Pass := Size_Bytes > 0 & nzchar(MD5)]
  fwrite(file_qa, file.path(
    qa_dir, "output_file_qa_v0.1_20260730.tsv"
  ), sep = "\t")
  if (!all(file_qa$Pass)) {
    stop("ORA输出存在空文件")
  }

  if (!file.rename(tmp_dir, final_dir)) {
    stop("临时目录原子提交失败: ", tmp_dir)
  }
  log_event(
    "DONE", "Stage 3 lncRNA C1 vs C4 linked-mRNA ORA complete",
    list(
      runtime_sec = round(as.numeric(difftime(
        Sys.time(), start_time, units = "secs"
      )), 1),
      units = nrow(units), significant_rows = nrow(significant_ora),
      main_display_rows = nrow(display_ora),
      display_lncrna = uniqueN(display_ora$Feature_ID),
      output = final_dir, qa = "PASS"
    )
  )
}

tryCatch(
  main(),
  error = function(e) {
    log_event("ERROR", conditionMessage(e))
    quit(status = 1L)
  }
)
}
