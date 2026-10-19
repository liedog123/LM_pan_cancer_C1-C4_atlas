# Purpose: mrna programmes analysis stages.
# Configuration: config/molecular.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Fit mRNA state comparisons and derive recurrent genes and functional programmes.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/01_run_mrna_models.R
stage_transcriptome_proteome_01_run_mrna_models <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]

required_files <- input_paths[names(input_paths) != "model_dir"]

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

script_name <- "01_run_mrna_models.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_014, environment())
`norm_sample` <- bind_stage_helper(shared_norm_sample_037, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())


state_levels <- c("C1", "C2", "C3", "C4")

design_specs <- list(
  C1_vs_REST = list(
    type = "one_vs_rest_pooled_patients",
    target = "C1", reference = "REST_POOLED",
    evidence_role = "STATE_PORTRAIT_MAIN"
  ),
  C2_vs_REST = list(
    type = "one_vs_rest_pooled_patients",
    target = "C2", reference = "REST_POOLED",
    evidence_role = "STATE_PORTRAIT_MAIN"
  ),
  C3_vs_REST = list(
    type = "one_vs_rest_pooled_patients",
    target = "C3", reference = "REST_POOLED",
    evidence_role = "STATE_PORTRAIT_MAIN"
  ),
  C4_vs_REST = list(
    type = "one_vs_rest_pooled_patients",
    target = "C4", reference = "REST_POOLED",
    evidence_role = "STATE_PORTRAIT_MAIN"
  ),
  C1_vs_C4 = list(
    type = "pair",
    target = "C1", reference = "C4",
    evidence_role = "JOINT_ENDPOINT_SUPPORT"
  )
)

design_registry <- rbindlist(lapply(names(design_specs), function(id) {
  spec <- design_specs[[id]]
  data.table(
    Design_ID = id,
    Design_Type = spec$type,
    Model_Variable = "C_subtype",
    Target = spec$target,
    Reference = spec$reference,
    Candidate_Eligible = TRUE,
    Evidence_Role = spec$evidence_role,
    BH_Family = "WITHIN_CANCER_X_DESIGN",
    Purity_Adjustment = FALSE
  )
}))
`align_matrix_to_patient` <- bind_stage_helper(shared_align_matrix_to_patient_017, environment())

state_feature_counts <- function(y, meta) {
  valid <- is.finite(y)
  counts <- lapply(state_levels, function(state) {
    idx <- which(meta$C_subtype == state)
    if (!length(idx)) rep(0L, nrow(y)) else rowSums(valid[, idx, drop = FALSE])
  })
  names(counts) <- state_levels
  counts
}

run_design <- function(expr, meta, feature, cancer, design_id, spec) {
  if (spec$type == "one_vs_rest_pooled_patients") {
    keep <- meta$C_subtype %chin% state_levels
  } else {
    keep <- meta$C_subtype %chin% c(spec$reference, spec$target)
  }
  y <- expr[, keep, drop = FALSE]
  meta0 <- meta[keep]
  if (!ncol(y)) stop(cancer, " ", design_id, "无设计样本")

  n_state <- state_feature_counts(y, meta0)
  n_used <- rowSums(is.finite(y))
  if (spec$type == "one_vs_rest_pooled_patients") {
    n_target <- n_state[[spec$target]]
    n_reference <- n_used - n_target
    eligible <- Reduce(`&`, lapply(n_state, function(x) x >= 10L)) &
      n_used >= 60L
    grp <- factor(
      ifelse(meta0$C_subtype == spec$target, "TARGET", "REST"),
      levels = c("REST", "TARGET")
    )
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- c("REST", "TARGET")
    contrast <- c(REST = -1, TARGET = 1)
    failure_reason <- fifelse(
      eligible, NA_character_, "FOUR_STATE_SAMPLE_THRESHOLD_FAILED"
    )
  } else {
    n_target <- n_state[[spec$target]]
    n_reference <- n_state[[spec$reference]]
    eligible <- n_target >= 10L & n_reference >= 10L &
      (n_target + n_reference) >= 30L
    grp <- factor(
      meta0$C_subtype,
      levels = c(spec$reference, spec$target)
    )
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- c(spec$reference, spec$target)
    contrast <- setNames(c(-1, 1), c(spec$reference, spec$target))
    failure_reason <- fifelse(
      eligible, NA_character_, "PAIR_SAMPLE_THRESHOLD_FAILED"
    )
  }

  if (any(table(grp) == 0L)) {
    eligible[] <- FALSE
    failure_reason[] <- "TARGET_OR_REFERENCE_ABSENT"
  }
  fit <- lmFit(y, design)
  fit2 <- contrasts.fit(fit, contrasts = matrix(contrast, ncol = 1L))
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  estimate <- as.numeric(fit2$coefficients[, 1])
  se <- as.numeric(fit2$stdev.unscaled[, 1] * sqrt(fit2$s2.post))
  statistic <- as.numeric(fit2$t[, 1])
  p_value <- as.numeric(fit2$p.value[, 1])
  b_value <- as.numeric(fit2$lods[, 1])
  estimate[!eligible] <- NA_real_
  se[!eligible] <- NA_real_
  statistic[!eligible] <- NA_real_
  p_value[!eligible] <- NA_real_
  b_value[!eligible] <- NA_real_

  result <- cbind(
    data.table(Cancer = cancer, Layer = "MRNA"),
    feature,
    data.table(
      Design_ID = design_id,
      Design_Type = spec$type,
      Analysis_Set = "REPRESENTATIVE_MAIN",
      Model_Type = "UNADJUSTED",
      Target = spec$target,
      Reference = spec$reference,
      Candidate_Eligible = TRUE,
      Evidence_Role = spec$evidence_role,
      N_Input = ncol(y),
      N_Used = n_used,
      N_Target = n_target,
      N_Reference = n_reference,
      N_C1 = n_state$C1,
      N_C2 = n_state$C2,
      N_C3 = n_state$C3,
      N_C4 = n_state$C4,
      Estimate = estimate,
      SE = se,
      Statistic = statistic,
      P_Value = p_value,
      FDR = NA_real_,
      AveExpr = rowMeans(y, na.rm = TRUE),
      B = b_value,
      Status = fifelse(eligible, "ESTIMABLE", "NOT_ESTIMABLE_BY_CONTRACT"),
      Failure_Reason = failure_reason
    )
  )
  result[is.finite(P_Value), FDR := p.adjust(P_Value, method = "BH")]
  result
}

run_cancer <- function(cancer, mrna_list, patient, definition_symbols,
                       tmp_individual_dir, final_individual_dir) {
  patient_cancer <- patient[Cancer == cancer]
  if (is.null(mrna_list[[cancer]]$expression_data)) {
    stop("缺少mRNA表达矩阵: ", cancer)
  }
  raw <- as.matrix(mrna_list[[cancer]]$expression_data)
  storage.mode(raw) <- "double"
  if (any(raw < 0, na.rm = TRUE)) stop("mRNA出现负TPM: ", cancer)
  if (anyDuplicated(rownames(raw))) stop("mRNA gene symbol重复: ", cancer)

  aligned <- align_matrix_to_patient(raw, patient_cancer)
  keep_state <- aligned$meta$C_subtype %chin% state_levels
  aligned$meta <- aligned$meta[keep_state]
  aligned$mat <- aligned$mat[, keep_state, drop = FALSE]
  if (!ncol(aligned$mat)) stop("无canonical四状态样本: ", cancer)

  log_expr_all <- log2(aligned$mat + 1)
  keep_gene <- rowMeans(aligned$mat >= 1, na.rm = TRUE) >= 0.20 &
    apply(log_expr_all, 1L, IQR, na.rm = TRUE) > 0
  expr <- log_expr_all[keep_gene, , drop = FALSE]
  feature <- data.table(
    Feature_ID = rownames(expr),
    Feature_Label = rownames(expr),
    Feature_Class = "GENE",
    Definition_Overlap = norm_symbol(rownames(expr)) %in% definition_symbols
  )

  result <- rbindlist(lapply(names(design_specs), function(design_id) {
    run_design(
      expr, aligned$meta, feature, cancer, design_id,
      design_specs[[design_id]]
    )
  }), use.names = TRUE)
  setcolorder(result, c(
    "Cancer", "Layer", "Feature_ID", "Feature_Label", "Feature_Class",
    "Definition_Overlap", "Design_ID", "Design_Type", "Analysis_Set",
    "Model_Type", "Target", "Reference", "Candidate_Eligible",
    "Evidence_Role", "N_Input", "N_Used", "N_Target", "N_Reference",
    "N_C1", "N_C2", "N_C3", "N_C4", "Estimate", "SE", "Statistic",
    "P_Value", "FDR", "AveExpr", "B", "Status", "Failure_Reason"
  ))

  basename_out <- sprintf(
    "%s_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz", cancer
  )
  tmp_file <- file.path(tmp_individual_dir, basename_out)
  final_file <- file.path(final_individual_dir, basename_out)
  fwrite(
    result, tmp_file, sep = "\t", quote = TRUE, na = "NA",
    compress = "gzip"
  )
  if (!file.exists(tmp_file) || file.info(tmp_file)$size <= 0) {
    stop("逐癌输出为空: ", tmp_file)
  }

  coverage <- result[, .(
    N_Features = .N,
    N_Estimable = sum(Status == "ESTIMABLE"),
    N_FDR_005 = sum(FDR < 0.05, na.rm = TRUE),
    N_FDR_005_Effect_GT_1 = sum(
      FDR < 0.05 & is.finite(Estimate) & abs(Estimate) > 1,
      na.rm = TRUE
    )
  ), by = .(
    Cancer, Layer, Design_ID, Design_Type, Candidate_Eligible, Evidence_Role
  )]
  sample_counts <- table(factor(aligned$meta$C_subtype, levels = state_levels))
  manifest <- data.table(
    Cancer = cancer,
    Layer = "MRNA",
    N_Source_Features = nrow(raw),
    N_Filtered_Features = nrow(expr),
    N_Matched_Patients = ncol(expr),
    C1 = as.integer(sample_counts["C1"]),
    C2 = as.integer(sample_counts["C2"]),
    C3 = as.integer(sample_counts["C3"]),
    C4 = as.integer(sample_counts["C4"]),
    Four_State_Evaluable = all(sample_counts >= 10L) &&
      sum(sample_counts) >= 60L,
    Output = final_file,
    Temporary_Output = tmp_file
  )
  sample_manifest <- aligned$meta[, .(
    Cancer, Patient_ID, Sample, C_subtype, Acid_z, IG_z
  )]
  definition_audit <- data.table(
    Cancer = cancer,
    N_Filtered_Features = nrow(feature),
    N_Definition_Features_Retained = sum(feature$Definition_Overlap)
  )
  list(
    coverage = coverage,
    manifest = manifest,
    sample_manifest = sample_manifest,
    definition_audit = definition_audit
  )
}

main <- function() {
  args <- list(mode="formal", cancers="ALL")
  final_dir <- file.path(out_dir, "01_within_cancer_models")
  if (dir.exists(final_dir) &&
      length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录已存在且非空，拒绝覆盖: ", final_dir)
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  if (dir.exists(tmp_dir)) stop("临时目录已存在: ", tmp_dir)
  tmp_individual_dir <- file.path(tmp_dir, "individual")
  final_individual_dir <- file.path(final_dir, "individual")
  dir.create(tmp_individual_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "qa"), recursive = TRUE, showWarnings = FALSE)

  paths <- list(patient=unname(paths["patient_manifest"]),
                mrna=unname(paths["mrna_rdata"]),
                definition=unname(paths["definition_registry"]),
                gate=unname(paths["readiness_decision"]))
  if (any(!file.exists(unlist(paths)))) stop("缺少冻结输入")
  gate <- fread(paths$gate)
  if (nrow(gate) != 1L ||
      gate$decision != "PASS_STAGE3_READINESS_AUTHORIZE_FORMAL_ANALYSIS") {
    stop("Stage 3 readiness未授权")
  }

  patient <- fread(paths$patient)
  patient[, Sample := norm_sample(Sample)]
  if (anyDuplicated(patient$Patient_ID) || anyDuplicated(patient$Sample)) {
    stop("代表患者接口不唯一")
  }
  all_cancers <- sort(unique(patient$Cancer))
  cancers <- if (args$mode == "formal" && identical(args$cancers, "ALL")) {
    all_cancers
  } else {
    args$cancers
  }
  if (any(!cancers %in% all_cancers)) {
    stop("未知癌种: ", paste(setdiff(cancers, all_cancers), collapse = ","))
  }

  definition <- fread(paths$definition)
  definition_symbols <- norm_symbol(definition$gene_symbol)
  mrna_env <- new.env(parent = emptyenv())
  load(paths$mrna, envir = mrna_env)
  if (!exists("pancancer_data_list", envir = mrna_env, inherits = FALSE)) {
    stop("mRNA输入缺少pancancer_data_list")
  }
  mrna_list <- mrna_env$pancancer_data_list

  log_event("START", "mRNA five-design within-cancer models begin", list(
    project = "LM_pan_cancer",
    mode = args$mode,
    environment = "R4.5.1",
    cancers = paste(cancers, collapse = ","),
    designs = paste(names(design_specs), collapse = ","),
    output = final_dir,
    excluded_designs = "ATOMIC_EDGES,CONTINUOUS_AXES,OMNIBUS",
    patient_pooling_across_cancer = "PROHIBITED",
    purity_adjustment = "NONE"
  ))

  fwrite(
    design_registry,
    file.path(tmp_dir, "mrna_ovr_c1c4_design_registry_v0.2_20260729.tsv"),
    sep = "\t", quote = TRUE, na = "NA"
  )
  coverage_rows <- list()
  manifest_rows <- list()
  sample_rows <- list()
  definition_rows <- list()
  for (cancer in cancers) {
    run <- run_cancer(
      cancer, mrna_list, patient, definition_symbols,
      tmp_individual_dir, final_individual_dir
    )
    coverage_rows[[cancer]] <- run$coverage
    manifest_rows[[cancer]] <- run$manifest
    sample_rows[[cancer]] <- run$sample_manifest
    definition_rows[[cancer]] <- run$definition_audit
    log_event("CHECKPOINT", "cancer five-design models complete", list(
      cancer = cancer,
      n_features = run$manifest$N_Filtered_Features,
      n_patients = run$manifest$N_Matched_Patients,
      c1 = run$manifest$C1,
      c2 = run$manifest$C2,
      c3 = run$manifest$C3,
      c4 = run$manifest$C4
    ))
    gc()
  }

  coverage <- rbindlist(coverage_rows)
  manifest <- rbindlist(manifest_rows)
  sample_manifest <- unique(rbindlist(sample_rows))
  definition_audit <- rbindlist(definition_rows)
  fwrite(
    coverage,
    file.path(tmp_dir, "mrna_ovr_c1c4_model_coverage_v0.2_20260729.tsv"),
    sep = "\t", quote = TRUE, na = "NA"
  )
  manifest_export <- manifest[
    , setdiff(names(manifest), "Temporary_Output"), with = FALSE
  ]
  fwrite(
    manifest_export,
    file.path(tmp_dir, "mrna_ovr_c1c4_feature_manifest_v0.2_20260729.tsv"),
    sep = "\t", quote = TRUE, na = "NA"
  )
  fwrite(
    sample_manifest,
    file.path(tmp_dir, "mrna_ovr_c1c4_sample_manifest_v0.2_20260729.tsv.gz"),
    sep = "\t", quote = TRUE, na = "NA", compress = "gzip"
  )
  fwrite(
    definition_audit,
    file.path(
      tmp_dir, "mrna_ovr_c1c4_definition_retention_audit_v0.2_20260729.tsv"
    ),
    sep = "\t", quote = TRUE, na = "NA"
  )

  forbidden_patterns <- c(
    "^AG[123]_IG[12]_vs_", "^CONTINUOUS_", "OMNIBUS"
  )
  qa <- data.table(
    Check = c(
      "EXACTLY_FIVE_DESIGNS",
      "FOUR_POOLED_PATIENT_OVR",
      "ONE_C1_VS_C4",
      "NO_ATOMIC_CONTINUOUS_OMNIBUS",
      "NO_EQUAL_STATE_MEAN_OVR",
      "NO_PURITY_MODEL_OR_FIELD",
      "ONE_NONEMPTY_FILE_PER_CANCER",
      "ESTIMABLE_SE_FINITE",
      "FDR_RANGE_VALID",
      "WITHIN_CANCER_BH_FAMILY_FIXED",
      "DEFINITION_GENES_NOT_EXCLUDED"
    ),
    Pass = c(
      uniqueN(coverage$Design_ID) == 5L,
      sum(design_registry$Design_Type == "one_vs_rest_pooled_patients") == 4L,
      identical(design_registry[Design_ID == "C1_vs_C4", .N], 1L),
      !any(vapply(
        forbidden_patterns,
        function(pattern) any(grepl(pattern, design_registry$Design_ID)),
        logical(1)
      )),
      !any(grepl("EQUAL_MEAN", design_registry$Reference, fixed = TRUE)),
      !any(design_registry$Purity_Adjustment) &&
        !any(c("purity_z", "TumorPurity") %in% names(sample_manifest)),
      nrow(manifest) == length(cancers) &&
        all(file.exists(manifest$Temporary_Output)) &&
        all(file.info(manifest$Temporary_Output)$size > 0),
      all(vapply(manifest$Temporary_Output, function(path) {
        x <- fread(path, select = c("Status", "Estimate", "SE"))
        x <- x[Status == "ESTIMABLE" & is.finite(Estimate)]
        !nrow(x) || all(is.finite(x$SE))
      }, logical(1))),
      all(is.na(coverage$N_FDR_005) | coverage$N_FDR_005 >= 0L),
      all(design_registry$BH_Family == "WITHIN_CANCER_X_DESIGN"),
      all(definition_audit$N_Definition_Features_Retained >= 0L)
    )
  )
  fwrite(
    qa,
    file.path(
      tmp_dir, "qa",
      "mrna_ovr_c1c4_model_machine_qa_v0.2_20260729.tsv"
    ),
    sep = "\t", quote = TRUE
  )
  if (!all(qa$Pass)) stop("mRNA五设计逐癌模型机器QA失败")

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) stop("无法提交mRNA五设计模型目录")
  log_event("DONE", "mRNA five-design within-cancer models complete", list(
    runtime_sec = round(
      as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1
    ),
    mode = args$mode,
    n_cancers = length(cancers),
    n_designs = nrow(design_registry),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)),
    output = final_dir,
    evidence_role = if (args$mode == "formal") "MAINLINE" else "SMOKE_QC"
  ))
}

tryCatch(
  main(),
  error = function(e) {
    log_event("ERROR", conditionMessage(e), list(
      runtime_sec = round(
        as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1
      )
    ))
    quit(status = 1L)
  }
)
}

# Purpose: Fit mRNA state comparisons and derive recurrent genes and functional programmes.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/02_run_mrna_recurrence_ora.R
stage_transcriptome_proteome_02_run_mrna_recurrence_ora <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("model_dir must be readable/searchable.")

required_files <- input_paths[names(input_paths) != "model_dir"]
required_files <- c(required_files,file.path(paths["model_dir"],c("mrna_ovr_c1c4_design_registry_v0.2_20260729.tsv","mrna_ovr_c1c4_feature_manifest_v0.2_20260729.tsv","individual/BLCA_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/BRCA_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/CESC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/COAD_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/ESCA_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/GBM_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/HNSC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/KIRC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/KIRP_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/LGG_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/LIHC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/LUAD_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/LUSC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/OV_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/PCPG_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/PRAD_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/READ_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/SARC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/SKCM_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/STAD_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/THCA_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz","individual/UCEC_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz")))

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

suppressPackageStartupMessages(library(data.table))

options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "02_run_mrna_recurrence_ora.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_014, environment())


write_tsv <- function(x, path, compress = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(
    x, path, sep = "\t", quote = TRUE, na = "NA",
    compress = if (compress) "gzip" else "none"
  )
}

sha256_file <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(value) || !grepl("^[0-9a-fA-F]{64}", value[1])) {
    stop("无法计算SHA-256: ", path)
  }
  sub("\\s+.*$", "", value[1])
}

read_gmt <- function(path, collection) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(parts) < 3L)) {
    stop(collection, " GMT存在少于3列的记录")
  }
  registry <- rbindlist(lapply(seq_along(parts), function(i) {
    x <- parts[[i]]
    data.table(
      Collection = collection,
      Term_Name = trimws(x[1]),
      Description_URL = trimws(x[2]),
      Original_Line = i
    )
  }))
  membership <- rbindlist(lapply(seq_along(parts), function(i) {
    x <- parts[[i]]
    genes <- unique(toupper(trimws(x[-c(1, 2)])))
    genes <- genes[nzchar(genes)]
    data.table(
      Collection = collection,
      Term_Name = trimws(x[1]),
      Gene = genes
    )
  }))
  registry[, Term_Key := paste(Collection, Term_Name, sep = "::")]
  membership[, Term_Key := paste(Collection, Term_Name, sep = "::")]
  registry <- registry[order(Original_Line), .SD[1], by = Term_Key]
  membership <- unique(membership[, .(
    Term_Key, Collection, Term_Name, Gene
  )])
  list(registry = registry, membership = membership)
}

clean_display_label <- function(term, collection) {
  value <- term
  if (collection == "GO_BP") value <- sub("^GOBP_", "", value)
  if (collection == "REACTOME") value <- sub("^REACTOME_", "", value)
  if (collection == "HALLMARK") value <- sub("^HALLMARK_", "", value)
  if (collection != "IMM_PORT") {
    value <- tools::toTitleCase(tolower(gsub("_", " ", value, fixed = TRUE)))
  }
  trimws(value)
}

median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

collapse_cancers <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  paste(x, collapse = ";")
}

compute_ora <- function(design_id, direction, input_genes, universe,
                        registry, membership) {
  if (!length(input_genes)) return(data.table())
  bg_membership <- membership[Gene %chin% universe]
  term_size <- bg_membership[, .(Gene_Set_Size = uniqueN(Gene)), by = Term_Key]
  tested <- merge(registry, term_size, by = "Term_Key", all.x = TRUE)
  tested[is.na(Gene_Set_Size), Gene_Set_Size := 0L]
  tested <- tested[Gene_Set_Size >= 10L & Gene_Set_Size <= 500L]
  overlap <- bg_membership[
    Term_Key %chin% tested$Term_Key & Gene %chin% input_genes,
    .(
      Overlap_N = uniqueN(Gene),
      Overlap_Genes = paste(sort(unique(Gene)), collapse = ";")
    ),
    by = Term_Key
  ]
  ans <- merge(tested, overlap, by = "Term_Key", all.x = TRUE)
  ans[is.na(Overlap_N), `:=`(Overlap_N = 0L, Overlap_Genes = "")]
  universe_n <- length(universe)
  input_n <- length(input_genes)
  ans[, `:=`(
    Design_ID = design_id,
    Direction = direction,
    Universe_N = universe_n,
    Input_N = input_n,
    Expected_Overlap = input_n * Gene_Set_Size / universe_n,
    Fold_Enrichment = (Overlap_N / input_n) /
      (Gene_Set_Size / universe_n)
  )]
  ans[, P_Value := phyper(
    Overlap_N - 1L,
    Gene_Set_Size,
    Universe_N - Gene_Set_Size,
    Input_N,
    lower.tail = FALSE
  )]
  ans[, Odds_Ratio_Haldane := {
    a <- Overlap_N
    b <- Input_N - Overlap_N
    c <- Gene_Set_Size - Overlap_N
    d <- Universe_N - a - b - c
    ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  }]
  ans[, FDR := p.adjust(P_Value, method = "BH"), by = Collection]
  ans[, Full_Significant := FDR < 0.05 & Overlap_N >= 3L]
  ans[, Display_Eligible := Full_Significant &
        Gene_Set_Size <= 300L & Overlap_N >= 5L]
  ans[, Display_Label := mapply(
    clean_display_label, Term_Name, Collection, USE.NAMES = FALSE
  )]
  setorder(ans, Design_ID, Direction, Collection, FDR, -Overlap_N,
           Gene_Set_Size, Term_Name)
  ans
}

main <- function() {
  args <- list(mode="formal", cancers="ALL")
  final_dir <- file.path(out_dir, "02_recurrence_ora")
  model_root <- unname(paths["model_dir"])
  if (dir.exists(final_dir) &&
      length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录已存在且非空，拒绝覆盖: ", final_dir)
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  if (dir.exists(tmp_dir)) stop("临时目录已存在: ", tmp_dir)
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    design_registry = file.path(
      model_root, "mrna_ovr_c1c4_design_registry_v0.2_20260729.tsv"
    ),
    feature_manifest = file.path(
      model_root, "mrna_ovr_c1c4_feature_manifest_v0.2_20260729.tsv"
    )
  )
  if (any(!file.exists(unlist(paths)))) stop("缺少V0.2逐癌模型输出")
  design_registry <- fread(paths$design_registry)
  manifest <- fread(paths$feature_manifest)
  expected_names <- paste0(c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC"), "_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz")
  if(nrow(manifest)!=length(expected_names) || anyDuplicated(basename(manifest$Output)) ||
     !setequal(basename(manifest$Output),expected_names)) stop("Expected complete 22-cancer model manifest.")
  manifest[, Output := file.path(model_root,"individual",basename(Output))]

  if (!nrow(manifest) || any(!file.exists(manifest$Output))) {
    stop("V0.2逐癌模型manifest为空或文件缺失")
  }
  expected_designs <- c(
    "C1_vs_REST", "C2_vs_REST", "C3_vs_REST", "C4_vs_REST", "C1_vs_C4"
  )
  candidate_designs <- design_registry[
    Candidate_Eligible == TRUE, Design_ID
  ]
  if (!setequal(candidate_designs, expected_designs) ||
      length(candidate_designs) != 5L) {
    stop("V0.2方向性设计不是固定五设计")
  }

  keep_columns <- c(
    "Cancer", "Layer", "Feature_ID", "Feature_Label", "Feature_Class",
    "Definition_Overlap", "Design_ID", "Design_Type", "Target",
    "Reference", "Candidate_Eligible", "Evidence_Role", "N_Input",
    "N_Used", "N_Target", "N_Reference", "Estimate", "SE", "Statistic",
    "P_Value", "FDR", "Status"
  )
  model <- rbindlist(lapply(manifest$Output, function(path) {
    fread(path, select = keep_columns)
  }), use.names = TRUE)
  model <- model[Design_ID %chin% expected_designs]
  model[, Gene := toupper(trimws(Feature_ID))]
  if (any(!nzchar(model$Gene))) stop("存在空gene symbol")
  if (model[, anyDuplicated(paste(Cancer, Design_ID, Gene, sep = "::"))]) {
    stop("逐癌结果存在重复Cancer × Design × Gene")
  }
  model[, `:=`(
    Sig_Up = Status == "ESTIMABLE" & is.finite(FDR) & FDR < 0.05 &
      is.finite(Estimate) & Estimate > 1.0,
    Sig_Down = Status == "ESTIMABLE" & is.finite(FDR) & FDR < 0.05 &
      is.finite(Estimate) & Estimate < -1.0
  )]

  log_event("START", "mRNA five-design recurrence and ORA begin", list(
    project = "LM_pan_cancer",
    mode = args$mode,
    environment = "R4.5.1",
    cancers = uniqueN(model$Cancer),
    designs = paste(expected_designs, collapse = ","),
    within_cancer_fdr = 0.05,
    abs_effect = 1.0,
    recurrence_fraction = 0.30,
    cross_cancer_bh = "PROHIBITED",
    output = final_dir
  ))

  design_evaluable <- model[
    Status == "ESTIMABLE",
    .(N_Design_Evaluable_Cancers = uniqueN(Cancer)),
    by = Design_ID
  ]
  recurrence <- model[, .(
    Feature_Label = Feature_Label[1],
    Feature_Class = Feature_Class[1],
    Definition_Overlap = any(Definition_Overlap),
    Design_Type = Design_Type[1],
    Target = Target[1],
    Reference = Reference[1],
    Evidence_Role = Evidence_Role[1],
    N_Evaluable_Cancers = uniqueN(Cancer[Status == "ESTIMABLE"]),
    N_Sig_Up = uniqueN(Cancer[Sig_Up]),
    N_Sig_Down = uniqueN(Cancer[Sig_Down]),
    Sig_Up_Cancers = collapse_cancers(Cancer[Sig_Up]),
    Sig_Down_Cancers = collapse_cancers(Cancer[Sig_Down]),
    Median_Estimate_Up = median_or_na(Estimate[Sig_Up]),
    Median_Abs_Estimate_Up = median_or_na(abs(Estimate[Sig_Up])),
    Median_Estimate_Down = median_or_na(Estimate[Sig_Down]),
    Median_Abs_Estimate_Down = median_or_na(abs(Estimate[Sig_Down]))
  ), by = .(Design_ID, Gene)]
  recurrence <- merge(
    recurrence, design_evaluable, by = "Design_ID", all.x = TRUE
  )
  recurrence[
    is.na(N_Design_Evaluable_Cancers),
    N_Design_Evaluable_Cancers := 0L
  ]
  recurrence[, `:=`(
    Min_Recurrence_Cancers = ceiling(0.30 * N_Evaluable_Cancers),
    Recurrence_Up = N_Sig_Up / N_Evaluable_Cancers,
    Recurrence_Down = N_Sig_Down / N_Evaluable_Cancers
  )]
  recurrence[, `:=`(
    Recurrent_Up = N_Evaluable_Cancers > 0L &
      N_Sig_Up >= Min_Recurrence_Cancers,
    Recurrent_Down = N_Evaluable_Cancers > 0L &
      N_Sig_Down >= Min_Recurrence_Cancers,
    Common_Measurable = N_Design_Evaluable_Cancers > 0L &
      N_Evaluable_Cancers == N_Design_Evaluable_Cancers
  )]
  recurrence[, Bidirectional_Recurrent := Recurrent_Up & Recurrent_Down]
  recurrence[, Recurrence_Class := fifelse(
    Bidirectional_Recurrent, "BIDIRECTIONAL_RECURRENT",
    fifelse(
      Recurrent_Up, "RECURRENT_UP",
      fifelse(Recurrent_Down, "RECURRENT_DOWN", "NOT_RECURRENT")
    )
  )]
  setorder(recurrence, Design_ID, -N_Sig_Up, -N_Sig_Down, Gene)

  common_background <- recurrence[
    Common_Measurable == TRUE,
    .(
      Design_ID, Gene, Feature_Label, Definition_Overlap,
      N_Design_Evaluable_Cancers, N_Evaluable_Cancers
    )
  ]
  directional <- rbindlist(list(
    recurrence[, .(
      Design_ID, Gene, Feature_Label, Feature_Class, Definition_Overlap,
      Design_Type, Target, Reference, Evidence_Role,
      Direction = "UP",
      N_Design_Evaluable_Cancers, N_Evaluable_Cancers,
      Min_Recurrence_Cancers,
      N_Sig = N_Sig_Up,
      Recurrence = Recurrence_Up,
      Sig_Cancers = Sig_Up_Cancers,
      Median_Estimate = Median_Estimate_Up,
      Median_Abs_Estimate = Median_Abs_Estimate_Up,
      Recurrent = Recurrent_Up,
      Common_Measurable,
      Bidirectional_Recurrent
    )],
    recurrence[, .(
      Design_ID, Gene, Feature_Label, Feature_Class, Definition_Overlap,
      Design_Type, Target, Reference, Evidence_Role,
      Direction = "DOWN",
      N_Design_Evaluable_Cancers, N_Evaluable_Cancers,
      Min_Recurrence_Cancers,
      N_Sig = N_Sig_Down,
      Recurrence = Recurrence_Down,
      Sig_Cancers = Sig_Down_Cancers,
      Median_Estimate = Median_Estimate_Down,
      Median_Abs_Estimate = Median_Abs_Estimate_Down,
      Recurrent = Recurrent_Down,
      Common_Measurable,
      Bidirectional_Recurrent
    )]
  ), use.names = TRUE)
  directional[, ORA_Eligible := Recurrent & Common_Measurable &
                !Bidirectional_Recurrent]
  recurrent_directional <- directional[ORA_Eligible == TRUE]
  setorder(
    recurrent_directional, Design_ID, Direction, -N_Sig, -Recurrence,
    -Median_Abs_Estimate, -Median_Estimate, Gene
  )
  recurrent_directional[, Evidence_Rank := seq_len(.N),
                        by = .(Design_ID, Direction)]

  write_tsv(
    recurrence,
    file.path(
      tmp_dir, "recurrence",
      "mrna_ovr_c1c4_recurrence_all_v0.2_20260729.tsv.gz"
    ),
    compress = TRUE
  )
  write_tsv(
    recurrent_directional,
    file.path(
      tmp_dir, "recurrence",
      "mrna_ovr_c1c4_recurrent_directional_genes_v0.2_20260729.tsv"
    )
  )
  write_tsv(
    common_background,
    file.path(
      tmp_dir, "recurrence",
      "mrna_ovr_c1c4_common_background_v0.2_20260729.tsv.gz"
    ),
    compress = TRUE
  )
  evidence <- merge(
    recurrent_directional[, .(Design_ID, Gene, Direction, Evidence_Rank)],
    model,
    by = c("Design_ID", "Gene"),
    allow.cartesian = TRUE
  )
  evidence[, Significant_In_Direction := fifelse(
    Direction == "UP", Sig_Up, Sig_Down
  )]
  setorder(evidence, Design_ID, Direction, Evidence_Rank, Cancer)
  write_tsv(
    evidence,
    file.path(
      tmp_dir, "source_data",
      "mrna_ovr_c1c4_recurrent_gene_cancer_evidence_v0.2_20260729.tsv.gz"
    ),
    compress = TRUE
  )

  resource_specs <- data.table(
    Collection = c("HALLMARK", "REACTOME", "GO_BP", "C7", "IMM_PORT"),
    Path = c(
      unname(get("paths",envir=.GlobalEnv)["hallmark_gmt"]),
      unname(get("paths",envir=.GlobalEnv)["reactome_gmt"]),
      unname(get("paths",envir=.GlobalEnv)["go_bp_gmt"]),
      unname(get("paths",envir=.GlobalEnv)["c7_gmt"]),
      unname(get("paths",envir=.GlobalEnv)["immport_gmt"])
    ),
    Version = c(
      "MSIGDB_2025.1_HS", "MSIGDB_2025.1_HS", "MSIGDB_2024.1_HS",
      "MSIGDB_2025.1_HS", "LOCAL_FROZEN_VERSION_UNSPECIFIED"
    ),
    Evidence_Role = c(
      "SUPPORT", "SUPPORT", "MAIN_NAMING", "SUPPLEMENT_ONLY",
      "MAIN_IMMUNE_NAMING"
    )
  )
  if (any(!file.exists(resource_specs$Path))) stop("存在缺失GMT资源")
  resource_specs[, `:=`(
    File_Size = file.info(Path)$size,
    Modified_Time = format(file.info(Path)$mtime, "%Y-%m-%d %H:%M:%S"),
    SHA256 = vapply(Path, sha256_file, character(1))
  )]
  parsed <- lapply(seq_len(nrow(resource_specs)), function(i) {
    read_gmt(resource_specs$Path[i], resource_specs$Collection[i])
  })
  registry <- rbindlist(lapply(parsed, `[[`, "registry"), fill = TRUE)
  membership <- unique(rbindlist(
    lapply(parsed, `[[`, "membership"), fill = TRUE
  ))
  resource_specs <- merge(
    resource_specs,
    registry[, .(Term_Count = .N), by = Collection],
    by = "Collection", all.x = TRUE
  )
  write_tsv(
    resource_specs,
    file.path(
      tmp_dir, "resource_manifest",
      "gene_set_resource_manifest_v0.2_20260729.tsv"
    )
  )

  ora_jobs <- CJ(
    Design_ID = expected_designs,
    Direction = c("UP", "DOWN"),
    unique = TRUE
  )
  ora_list <- vector("list", nrow(ora_jobs))
  ora_job_summary <- vector("list", nrow(ora_jobs))
  for (i in seq_len(nrow(ora_jobs))) {
    design_id <- ora_jobs$Design_ID[i]
    direction <- ora_jobs$Direction[i]
    universe <- sort(unique(
      common_background[Design_ID == design_id, Gene]
    ))
    input_genes <- sort(unique(recurrent_directional[
      Design_ID == design_id & Direction == direction, Gene
    ]))
    ora_list[[i]] <- compute_ora(
      design_id, direction, input_genes, universe, registry, membership
    )
    ora_job_summary[[i]] <- data.table(
      Design_ID = design_id,
      Direction = direction,
      Universe_N = length(universe),
      Input_N = length(input_genes),
      Tested_Terms = nrow(ora_list[[i]]),
      Significant_Terms = if (nrow(ora_list[[i]])) {
        sum(ora_list[[i]]$Full_Significant)
      } else {
        0L
      },
      Display_Eligible_Terms = if (nrow(ora_list[[i]])) {
        sum(ora_list[[i]]$Display_Eligible)
      } else {
        0L
      }
    )
    log_event("CHECKPOINT", "design-direction ORA complete", list(
      design = design_id,
      direction = direction,
      universe_n = length(universe),
      input_n = length(input_genes),
      significant_terms = ora_job_summary[[i]]$Significant_Terms
    ))
  }
  ora <- rbindlist(ora_list, use.names = TRUE, fill = TRUE)
  ora_summary <- rbindlist(ora_job_summary)
  if (!nrow(ora)) stop("所有方向均无ORA输入或结果")
  significant_ora <- ora[Full_Significant == TRUE]
  write_tsv(
    ora,
    file.path(
      tmp_dir, "ora",
      "mrna_ovr_c1c4_ora_all_v0.2_20260729.tsv.gz"
    ),
    compress = TRUE
  )
  write_tsv(
    significant_ora,
    file.path(
      tmp_dir, "ora",
      "mrna_ovr_c1c4_ora_significant_v0.2_20260729.tsv"
    )
  )

  significant_membership <- merge(
    significant_ora[, .(
      Design_ID, Direction, Term_Key, Collection, Term_Name, FDR,
      Full_Significant, Display_Eligible
    )],
    membership,
    by = c("Term_Key", "Collection", "Term_Name"),
    allow.cartesian = TRUE
  )
  significant_membership <- merge(
    significant_membership,
    recurrent_directional[, .(
      Design_ID, Direction, Gene, Evidence_Rank, Definition_Overlap
    )],
    by = c("Design_ID", "Direction", "Gene")
  )
  setorder(
    significant_membership, Design_ID, Direction, Collection, FDR,
    Term_Name, Evidence_Rank, Gene
  )
  write_tsv(
    significant_membership,
    file.path(
      tmp_dir, "source_data",
      "mrna_ovr_c1c4_significant_term_gene_membership_v0.2_20260729.tsv.gz"
    ),
    compress = TRUE
  )

  recurrence_summary <- directional[, .(
    N_Evaluable_Genes = sum(N_Evaluable_Cancers > 0L),
    N_Common_Measurable_Genes = sum(Common_Measurable),
    N_Recurrent_Genes = sum(Recurrent),
    N_Bidirectional_Recurrent_Genes = sum(Bidirectional_Recurrent),
    N_ORA_Input_Genes = sum(ORA_Eligible),
    N_Definition_ORA_Input_Genes = sum(ORA_Eligible & Definition_Overlap)
  ), by = .(Design_ID, Direction)]
  summary_table <- merge(
    recurrence_summary, ora_summary,
    by = c("Design_ID", "Direction"), all = TRUE
  )
  write_tsv(
    summary_table,
    file.path(
      tmp_dir, "summary",
      "mrna_ovr_c1c4_recurrence_ora_summary_v0.2_20260729.tsv"
    )
  )
  top_terms <- significant_ora[
    order(Design_ID, Direction, Collection, FDR, -Overlap_N, Gene_Set_Size),
    head(.SD, 20L),
    by = .(Design_ID, Direction, Collection)
  ]
  write_tsv(
    top_terms,
    file.path(
      tmp_dir, "summary",
      "mrna_ovr_c1c4_ora_top_terms_v0.2_20260729.tsv"
    )
  )

  recurrence_count_exact <- sum(recurrence$N_Sig_Up) == sum(model$Sig_Up) &&
    sum(recurrence$N_Sig_Down) == sum(model$Sig_Down)
  forbidden_columns <- c(
    "Meta_FDR", "I2", "Tau2", "LOO_Direction_Stable",
    "Direction_Consistency", "Max_Random_Weight", "Pooled_Estimate"
  )
  forbidden_designs <- grepl(
    "^AG[123]_IG[12]_vs_|^CONTINUOUS_|OMNIBUS",
    recurrence$Design_ID
  )
  qa <- data.table(
    Check = c(
      "EXACTLY_FIVE_DIRECTIONAL_DESIGNS",
      "TEN_DESIGN_DIRECTION_JOBS",
      "NO_ATOMIC_CONTINUOUS_OMNIBUS",
      "RECURRENCE_COUNTS_EXACT",
      "NO_CROSS_CANCER_BH_OR_META_FIELDS",
      "COMMON_BACKGROUND_REQUIRES_ALL_EVALUABLE_CANCERS",
      "BIDIRECTIONAL_GENES_EXCLUDED_FROM_ORA",
      "FIVE_GENESET_COLLECTIONS_PRESENT",
      "ORA_FDR_RANGE_VALID",
      "DEFINITION_FLAG_PRESERVED",
      "OUTPUT_TABLES_NONEMPTY"
    ),
    Pass = c(
      uniqueN(recurrence$Design_ID) == 5L &&
        setequal(unique(recurrence$Design_ID), expected_designs),
      nrow(ora_jobs) == 10L,
      !any(forbidden_designs),
      recurrence_count_exact,
      !any(forbidden_columns %in% names(recurrence)) &&
        !any(forbidden_columns %in% names(ora)),
      all(common_background$N_Evaluable_Cancers ==
            common_background$N_Design_Evaluable_Cancers),
      !any(recurrent_directional$Bidirectional_Recurrent),
      uniqueN(registry$Collection) == 5L,
      all(is.finite(ora$FDR) & ora$FDR >= 0 & ora$FDR <= 1),
      "Definition_Overlap" %in% names(recurrent_directional) &&
        "Definition_Overlap" %in% names(evidence),
      nrow(recurrence) > 0L && nrow(recurrent_directional) > 0L &&
        nrow(common_background) > 0L && nrow(evidence) > 0L &&
        nrow(ora) > 0L
    )
  )
  write_tsv(
    qa,
    file.path(
      tmp_dir, "qa",
      "mrna_ovr_c1c4_recurrence_ora_machine_qa_v0.2_20260729.tsv"
    )
  )
  if (!all(qa$Pass)) stop("mRNA五设计复现与ORA机器QA失败")

  report <- c(
    "# Stage 3 mRNA四状态one-vs-rest与C1 vs C4复现结果 V0.2",
    "",
    "日期：2026-07-29",
    "",
    paste0(
      "状态：`", if (args$mode == "formal") "FORMAL" else "SMOKE",
      " / FIVE_DESIGNS_ONLY / MACHINE_QA_PASS / MANUAL_REVIEW_PENDING`"
    ),
    "",
    "## 方法边界",
    "",
    paste0(
      "仅运行C1/C2/C3/C4 vs REST和C1 vs C4五个设计。所有差异均在",
      uniqueN(model$Cancer),
      "个单癌种内部运行并在Cancer × Design内做BH；不含原子边、continuous轴、omnibus、跨癌BH、meta或purity校正。"
    ),
    "",
    "方向性候选要求癌种内FDR<0.05、绝对效应>1.0并在至少30%可评价癌种中同方向显著。双向复现基因保留审计但不进入ORA。",
    "",
    "## 数量接口",
    "",
    paste(
      apply(summary_table, 1L, function(x) {
        sprintf(
          "%s %s：共同背景%s个，ORA输入%s个，显著term %s条，初始展示term %s条。",
          x[["Design_ID"]], x[["Direction"]],
          x[["N_Common_Measurable_Genes"]],
          x[["N_ORA_Input_Genes"]],
          x[["Significant_Terms"]],
          x[["Display_Eligible_Terms"]]
        )
      }),
      collapse = "\n\n"
    ),
    "",
    "## 下一步",
    "",
    "Current formal selection uses C1 vs C4 only. Four OVR designs are retained historical support, not formal feature/display selection."
  )
  writeLines(
    report,
    file.path(
      tmp_dir,
      "stage3_mrna_ovr_c1c4_recurrence_ora_report_v0.2_20260729.md"
    )
  )

  output_files <- list.files(tmp_dir, recursive = TRUE, full.names = TRUE)
  if (!length(output_files) || any(file.info(output_files)$size <= 0)) {
    stop("存在缺失或空输出文件")
  }
  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) {
    stop("无法提交mRNA五设计复现与ORA目录")
  }
  log_event("DONE", "mRNA five-design recurrence and ORA complete", list(
    runtime_sec = round(
      as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1
    ),
    mode = args$mode,
    n_cancers = uniqueN(model$Cancer),
    n_designs = uniqueN(recurrence$Design_ID),
    recurrent_directional_genes = nrow(recurrent_directional),
    significant_terms = nrow(significant_ora),
    display_eligible_terms = sum(ora$Display_Eligible),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)),
    output = final_dir,
    manual_review = "PENDING"
  ))
}

tryCatch(
  main(),
  error = function(e) {
    log_event("ERROR", conditionMessage(e), list(
      runtime_sec = round(
        as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1
      )
    ))
    quit(status = 1L)
  }
)
}

# Purpose: Fit mRNA state comparisons and derive recurrent genes and functional programmes.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/16_prepare_mrna_representatives.R
stage_transcriptome_proteome_16_prepare_mrna_representatives <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

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


options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))
# 固定人工代表项来自source37；不重新挑选term或重算ORA/BH。
cluster_registry <- data.table(
  Direction = c(rep("UP", 8), rep("DOWN", 2)),
  State = c(rep("C1", 8), rep("C4", 2)),
  Cluster_ID = c(
    sprintf("C1_FC%02d", 1:8),
    sprintf("C4_FC%02d", 1:2)
  ),
  Cluster_Order = 1:10,
  Cluster_Label = c(
    "Cell cycle and chromosome organization",
    "Hypoxia and glycolytic metabolic adaptation",
    "ECM, integrin and cell adhesion",
    "EMT, motility and tissue morphogenesis",
    "Growth-factor, RTK and mTOR/TGF-beta signaling",
    "Angiogenesis, vascular remodeling and hemostasis",
    "Inflammatory, cytokine and IFN-associated remodeling",
    "Protein secretion, Golgi-ER trafficking and stress response",
    "Aerobic respiration and oxidative energy metabolism",
    "Electron transport, oxidative phosphorylation and ATP synthesis"
  ),
  Parent_Domain = c(
    rep(NA_character_, 8),
    rep("Mitochondrial respiratory metabolism", 2)
  ),
  Cluster_Color = c(
    "#CC79A7", "#E69F00", "#009E73", "#D55E00",
    "#7B61A8", "#56B4E9", "#B8860B", "#7A7A7A",
    "#0072B2", "#4C78A8"
  ),
  Semantic_Boundary = c(
    "Cell-cycle-associated bulk expression; not a proliferative dependency.",
    "Hypoxia/glycolysis-associated enrichment; not proof that lactate or low pH drives C1.",
    "Bulk ECM/adhesion association may include stromal composition; not direct matrix causality.",
    "Bulk EMT/motility annotation; not proof of tumor-cell EMT or invasion.",
    "Growth-factor pathway-associated expression; not evidence of pathway addiction.",
    "Vascular/hemostatic association; not proof of angiogenic dependence.",
    "Inflammatory/IFN-associated expression; not direct immune-cell activation or causal signaling.",
    "Secretory/trafficking/stress-associated expression; not proof of secretory dependency.",
    "C4-high mitochondrial respiration annotation supported by a small recurrent gene set.",
    "C4-high electron-transport annotation; related to, not independent from, C4_FC01."
  )
)

reg_term <- function(
    direction, cluster_id, term_key, display_label, term_order,
    semantic_boundary = NA_character_) {
  data.table(
    Direction = direction,
    Cluster_ID = cluster_id,
    Term_Key = term_key,
    Display_Term_Label = display_label,
    Term_Order = term_order,
    Manual_Review_Status = "FROZEN_SAFE_LABEL",
    Semantic_Boundary = semantic_boundary
  )
}

# 显式人工注册表：只选已核验的非C7正式显著term。
term_registry <- rbindlist(list(
  reg_term("UP", "C1_FC01", "HALLMARK::HALLMARK_G2M_CHECKPOINT",
           "G2/M checkpoint-associated expression", 1),
  reg_term("UP", "C1_FC01", "HALLMARK::HALLMARK_MITOTIC_SPINDLE",
           "Mitotic spindle organization", 2),
  reg_term("UP", "C1_FC01", "GO_BP::GOBP_SISTER_CHROMATID_SEGREGATION",
           "Sister chromatid segregation", 3),
  reg_term("UP", "C1_FC02", "HALLMARK::HALLMARK_HYPOXIA",
           "Hypoxia-associated expression", 1),
  reg_term("UP", "C1_FC02", "HALLMARK::HALLMARK_GLYCOLYSIS",
           "Glycolytic metabolism", 2),
  reg_term("UP", "C1_FC02", "GO_BP::GOBP_RESPONSE_TO_OXYGEN_LEVELS",
           "Response to oxygen levels", 3),
  reg_term("UP", "C1_FC02", "GO_BP::GOBP_PYRUVATE_METABOLIC_PROCESS",
           "Pyruvate metabolism", 4),
  reg_term("UP", "C1_FC03", "REACTOME::REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
           "Extracellular matrix organization", 1),
  reg_term("UP", "C1_FC03", "REACTOME::REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
           "Integrin cell-surface interactions", 2),
  reg_term("UP", "C1_FC03", "REACTOME::REACTOME_COLLAGEN_FORMATION",
           "Collagen formation", 3),
  reg_term("UP", "C1_FC04", "HALLMARK::HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "EMT-associated expression", 1),
  reg_term("UP", "C1_FC04", "GO_BP::GOBP_WOUND_HEALING",
           "Wound-healing-associated program", 2),
  reg_term("UP", "C1_FC04", "GO_BP::GOBP_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
           "Regulation of epithelial cell migration", 3),
  reg_term("UP", "C1_FC05", "HALLMARK::HALLMARK_MTORC1_SIGNALING",
           "mTORC1-signaling-associated expression", 1),
  reg_term("UP", "C1_FC05", "HALLMARK::HALLMARK_TGF_BETA_SIGNALING",
           "TGF-beta-signaling-associated expression", 2),
  reg_term("UP", "C1_FC05", "GO_BP::GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
           "Response to TGF-beta", 3),
  reg_term(
    "UP", "C1_FC05",
    "GO_BP::GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS",
    "Regulation of growth-factor response", 4
  ),
  reg_term("UP", "C1_FC06", "GO_BP::GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT",
           "Regulation of vasculature development", 1),
  reg_term("UP", "C1_FC06", "HALLMARK::HALLMARK_ANGIOGENESIS",
           "Angiogenesis-associated expression", 2),
  reg_term(
    "UP", "C1_FC06",
    "REACTOME::REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
    "Platelet activation and aggregation annotation", 3
  ),
  reg_term("UP", "C1_FC07", "HALLMARK::HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "TNF-NF-kappaB-associated expression", 1),
  reg_term("UP", "C1_FC07", "HALLMARK::HALLMARK_INFLAMMATORY_RESPONSE",
           "Inflammatory-response-associated expression", 2),
  reg_term("UP", "C1_FC07", "REACTOME::REACTOME_INTERFERON_SIGNALING",
           "Interferon-signaling-associated expression", 3),
  reg_term("UP", "C1_FC07", "GO_BP::GOBP_MYELOID_LEUKOCYTE_DIFFERENTIATION",
           "Myeloid leukocyte differentiation", 4),
  reg_term("UP", "C1_FC08", "HALLMARK::HALLMARK_PROTEIN_SECRETION",
           "Protein secretion", 1),
  reg_term(
    "UP", "C1_FC08",
    "REACTOME::REACTOME_INTRA_GOLGI_AND_RETROGRADE_GOLGI_TO_ER_TRAFFIC",
    "Intra-Golgi and Golgi-to-ER trafficking", 2
  ),
  reg_term("UP", "C1_FC08", "HALLMARK::HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
           "Unfolded-protein-response-associated expression", 3),
  reg_term(
    "UP", "C1_FC08",
    "GO_BP::GOBP_ENDOPLASMIC_RETICULUM_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT",
    "ER-to-Golgi vesicle transport", 4
  ),
  reg_term("DOWN", "C4_FC01", "GO_BP::GOBP_AEROBIC_RESPIRATION",
           "Aerobic respiration", 1),
  reg_term("DOWN", "C4_FC01", "GO_BP::GOBP_CELLULAR_RESPIRATION",
           "Cellular respiration", 2),
  reg_term(
    "DOWN", "C4_FC02",
    "GO_BP::GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
    "ATP synthesis-coupled electron transport", 1
  ),
  reg_term("DOWN", "C4_FC02", "GO_BP::GOBP_OXIDATIVE_PHOSPHORYLATION",
           "Oxidative phosphorylation", 2),
  reg_term(
    "DOWN", "C4_FC02",
    "REACTOME::REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "Respiratory electron transport", 3
  )
))

first_token <- function(x) {
  vapply(strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE), function(z) {
    if (!length(z) || !nzchar(z[[1]])) NA_character_ else z[[1]]
  }, character(1))
}


prepare_v1_display <- function(ora, genes, membership, upstream_qa) {
  required_ora <- c(
    "Term_Key", "Collection", "Term_Name", "Gene_Set_Size", "Overlap_N",
    "Design_ID", "Direction", "Fold_Enrichment", "FDR",
    "Full_Significant", "Display_Eligible"
  )
  required_gene <- c(
    "Design_ID", "Gene", "Direction", "N_Evaluable_Cancers", "N_Sig",
    "Recurrence", "Median_Estimate", "Median_Abs_Estimate",
    "Definition_Overlap", "ORA_Eligible", "Evidence_Rank"
  )
  required_membership <- c(
    "Design_ID", "Direction", "Gene", "Term_Key", "Collection",
    "FDR", "Full_Significant", "Display_Eligible"
  )
  if (!all(required_ora %in% names(ora))) stop("ORA输入列不完整")
  if (!all(required_gene %in% names(genes))) stop("复现基因输入列不完整")
  if (!all(required_membership %in% names(membership))) {
    stop("term-gene membership输入列不完整")
  }
  if ("Pass" %in% names(upstream_qa) && any(!upstream_qa$Pass)) {
    stop("上游复现/ORA机器QA存在未通过项")
  }

  main_terms <- copy(ora[
    Design_ID == "C1_vs_C4" &
      Full_Significant &
      Display_Eligible &
      Gene_Set_Size <= 300L &
      Overlap_N >= 5L &
      Collection != "C7"
  ])
  setorder(main_terms, Direction, FDR, -Overlap_N, -Fold_Enrichment, Term_Key)
  if (main_terms[Direction == "UP", .N] != 990L) {
    stop("C1-high非C7展示候选数量偏离冻结输入: ", main_terms[Direction == "UP", .N])
  }
  if (main_terms[Direction == "DOWN", .N] != 8L) {
    stop("C4-high非C7展示候选数量偏离冻结输入: ", main_terms[Direction == "DOWN", .N])
  }
  selected_terms <- merge(
    term_registry,
    cluster_registry[
      ,
      .(
        Direction, State, Cluster_ID, Cluster_Order, Cluster_Label,
        Parent_Domain, Cluster_Color,
        Cluster_Semantic_Boundary = Semantic_Boundary
      )
    ],
    by = c("Direction", "Cluster_ID"),
    all.x = TRUE,
    sort = FALSE
  )
  selected_terms <- merge(
    selected_terms,
    main_terms,
    by = c("Direction", "Term_Key"),
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_Manual", "")
  )
  if (nrow(selected_terms) != nrow(term_registry) ||
      any(is.na(selected_terms$Collection))) {
    stop("人工代表term未全部匹配正式候选")
  }
  if (any(selected_terms$Collection == "C7") ||
      any(!selected_terms$Full_Significant) ||
      any(!selected_terms$Display_Eligible)) {
    stop("人工代表term存在非正式主展示记录")
  }
  setorder(selected_terms, Cluster_Order, Term_Order)
  selected_terms[, Term_Node_ID := paste("TERM", Cluster_ID, Term_Key, sep = "::")]

  selected_membership <- unique(membership[
    Design_ID == "C1_vs_C4" &
      Term_Key %chin% selected_terms$Term_Key,
    .(Direction, Gene, Term_Key, Collection)
  ])
  selected_membership <- merge(
    selected_membership,
    selected_terms[
      ,
      .(
        Direction, Term_Key, Cluster_ID, Cluster_Order, Cluster_Label,
        Term_Order, Term_Node_ID, Display_Term_Label
      )
    ],
    by = c("Direction", "Term_Key"),
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(selected_membership$Cluster_ID))) {
    stop("代表term membership映射失败")
  }

  gene_evidence <- copy(genes[Design_ID == "C1_vs_C4" & ORA_Eligible])
  gene_support <- selected_membership[
    ,
    .(
      N_Selected_Terms = uniqueN(Term_Key),
      N_Selected_Clusters = uniqueN(Cluster_ID),
      Selected_Term_Keys = paste(sort(unique(Term_Key)), collapse = ";"),
      Selected_Cluster_IDs = paste(sort(unique(Cluster_ID)), collapse = ";")
    ),
    by = .(Direction, Gene)
  ]
  gene_candidates <- merge(
    gene_support,
    gene_evidence,
    by = c("Direction", "Gene"),
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(gene_candidates$N_Sig))) {
    stop("代表term成员基因未全部匹配复现证据")
  }
  setorderv(
    gene_candidates,
    c(
      "N_Sig", "Recurrence", "Median_Abs_Estimate",
      "N_Selected_Terms", "Gene"
    ),
    c(-1L, -1L, -1L, -1L, 1L),
    na.last = TRUE
  )
  gene_candidates[, Display_Evidence_Rank := seq_len(.N)]

  # 每个term至少3个、每个cluster至少15个C1基因；C4所有可连接基因保留。
  term_seeds <- selected_membership[
    gene_candidates[
      ,
      .(Direction, Gene, Display_Evidence_Rank)
    ],
    on = .(Direction, Gene),
    nomatch = 0L
  ][
    order(Cluster_Order, Term_Order, Display_Evidence_Rank, Gene),
    head(.SD, min(3L, .N)),
    by = .(Direction, Term_Key)
  ][, .(Direction, Gene)]

  cluster_seeds <- selected_membership[
    gene_candidates[
      ,
      .(Direction, Gene, Display_Evidence_Rank)
    ],
    on = .(Direction, Gene),
    nomatch = 0L
  ][
    Direction == "UP",
    .SD[order(Display_Evidence_Rank, Gene)][
      !duplicated(Gene)
    ][seq_len(min(15L, .N))],
    by = Cluster_ID
  ][, .(Direction, Gene)]

  c4_connected <- gene_candidates[Direction == "DOWN", .(Direction, Gene)]
  seed_pairs <- unique(rbindlist(list(term_seeds, cluster_seeds, c4_connected)))
  if (nrow(seed_pairs) > 240L) {
    stop("最小覆盖种子超过240个gene node")
  }

  target_gene_n <- 220L
  seed_keys <- paste(seed_pairs$Direction, seed_pairs$Gene, sep = "\r")
  remaining <- gene_candidates[
    !(paste(Direction, Gene, sep = "\r") %chin% seed_keys)
  ][order(Display_Evidence_Rank)]
  fill_n <- max(0L, target_gene_n - nrow(seed_pairs))
  selected_pairs <- unique(rbindlist(list(
    seed_pairs,
    remaining[seq_len(min(fill_n, .N)), .(Direction, Gene)]
  )))
  if (nrow(selected_pairs) < 200L || nrow(selected_pairs) > 240L) {
    stop("最终展示gene node数量不在200–240: ", nrow(selected_pairs))
  }

  display_genes <- merge(
    selected_pairs,
    gene_candidates,
    by = c("Direction", "Gene"),
    all.x = TRUE,
    sort = FALSE
  )
  setorderv(
    display_genes,
    c(
      "Direction", "N_Sig", "Recurrence", "Median_Abs_Estimate",
      "N_Selected_Terms", "Gene"
    ),
    c(1L, -1L, -1L, -1L, -1L, 1L)
  )
  display_genes[, State := fifelse(Direction == "UP", "C1", "C4")]
  display_genes[, Gene_Node_ID := paste("GENE", Gene, sep = "::")]
  display_genes[, Selected_For_Cytoscape := TRUE]
  display_genes[, Display_Gene_Rank := seq_len(.N)]

  final_membership <- selected_membership[
    display_genes[
      ,
      .(Direction, Gene, Gene_Node_ID)
    ],
    on = .(Direction, Gene),
    nomatch = 0L
  ]
  final_membership <- unique(final_membership)
  final_membership[, Selected_For_Cytoscape := TRUE]

  directional_gene_freeze <- copy(gene_evidence)
  directional_gene_freeze[, State := fifelse(Direction == "UP", "C1", "C4")]
  directional_gene_freeze[, Connected_To_Representative_Term :=
    paste(Direction, Gene, sep = "\r") %chin%
      paste(gene_candidates$Direction, gene_candidates$Gene, sep = "\r")]
  directional_gene_freeze[, Selected_For_Cytoscape :=
    paste(Direction, Gene, sep = "\r") %chin%
      paste(display_genes$Direction, display_genes$Gene, sep = "\r")]
  setorder(directional_gene_freeze, Direction, Evidence_Rank, Gene)


  list(display = display_genes, directional = directional_gene_freeze,
       terms = selected_terms, membership = selected_membership,
       final_membership = final_membership)
}
# 忠实保留source39的C4合簇、Other和最终字段/排序变换。
refine_v2_display <- function(display_v1, directional) {
  c4_cluster_id <- "C4_FC01"
  other_cluster_id <- "OTHER"
  display <- copy(display_v1)
  display[Direction == "DOWN", Selected_Cluster_IDs := c4_cluster_id]
  display[, N_Selected_Clusters := uniqueN(na.omit(unlist(strsplit(
    paste(Selected_Cluster_IDs, collapse = ";"), ";", fixed = TRUE
  )))), by = Gene]
  display[, `:=`(
    Display_Assignment = "FUNCTIONAL_TERM",
    Functional_Term_Mapped = TRUE,
    Other_Reason = NA_character_
  )]

  other_genes <- directional[
    State == "C4" & Recurrent & ORA_Eligible & !Connected_To_Representative_Term
  ]
  if (nrow(other_genes) != 10L) {
    stop("预期10个C4未连接代表term基因，实际为: ", nrow(other_genes))
  }

  new_display <- other_genes[, .(
    Direction,
    Gene,
    N_Selected_Terms = 0L,
    N_Selected_Clusters = 0L,
    Selected_Term_Keys = NA_character_,
    Selected_Cluster_IDs = NA_character_,
    Design_ID,
    Feature_Label,
    Feature_Class,
    Definition_Overlap,
    Design_Type,
    Target,
    Reference,
    Evidence_Role,
    N_Design_Evaluable_Cancers,
    N_Evaluable_Cancers,
    Min_Recurrence_Cancers,
    N_Sig,
    Recurrence,
    Sig_Cancers,
    Median_Estimate,
    Median_Abs_Estimate,
    Recurrent,
    Common_Measurable,
    Bidirectional_Recurrent,
    ORA_Eligible,
    Evidence_Rank,
    Display_Evidence_Rank = Evidence_Rank,
    State,
    Gene_Node_ID = paste0("GENE::", Gene),
    Selected_For_Cytoscape = TRUE,
    Display_Gene_Rank = NA_integer_,
    Display_Assignment = "OTHER",
    Functional_Term_Mapped = FALSE,
    Other_Reason = "FORMAL_RECURRENT_GENE_NO_DISPLAYED_REPRESENTATIVE_TERM"
  )]

  display <- rbindlist(list(display, new_display), use.names = TRUE, fill = TRUE)
  setorder(display, State, -N_Sig, -Recurrence, -Median_Abs_Estimate, Gene)
  display[, Display_Gene_Rank := seq_len(.N), by = State]
  display[, Primary_Cluster_ID := fifelse(
    Display_Assignment == "OTHER",
    other_cluster_id,
    first_token(Selected_Cluster_IDs)
  )]


  display
}
main <- function() {
  ora <- fread(paths["ora_significant"])
  genes <- fread(paths["recurrent_directional_genes"])
  membership <- fread(paths["significant_term_gene_membership"])
  upstream_qa <- fread(paths["upstream_qa"])
  if (!("Pass" %in% names(upstream_qa)) || !nrow(upstream_qa) ||
      anyNA(upstream_qa$Pass) || !all(upstream_qa$Pass)) stop("Upstream QA must be nonempty and all passed.")
  v1 <- prepare_v1_display(ora, genes, membership, upstream_qa)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  final_dir <- file.path(out_dir,"01_mrna_display_preparation")
  stage_dir <- paste0(final_dir,".tmp_",Sys.getpid())
  if(file.exists(final_dir) || file.exists(stage_dir)) stop("Output stage already exists.")
  if(!dir.create(stage_dir,showWarnings=FALSE)) stop("Cannot create temporary stage.")
  # 保留37→39历史TSV序列化边界，仅写两个必要v0.1中间表。
  intermediate <- file.path(stage_dir,"intermediate_v01")
  dir.create(intermediate,showWarnings=FALSE)
  v1_display_path <- file.path(intermediate,"mrna_c1c4_extreme_display_genes_v0.1_20260730.tsv")
  v1_directional_path <- file.path(intermediate,"mrna_c1c4_extreme_directional_gene_freeze_v0.1_20260730.tsv")
  fwrite(v1$display,v1_display_path,sep="\t",quote=TRUE,na="NA")
  fwrite(v1$directional,v1_directional_path,sep="\t",quote=TRUE,na="NA")
  display_v1 <- fread(v1_display_path)
  directional <- fread(v1_directional_path)
  display <- refine_v2_display(display_v1, directional)
  terms <- v1$terms
  # 仅保留代表项输入相关QA，不声称运行网络/社区/legacy QA。
  qa <- rbindlist(list(
    data.table(Check_ID="UPSTREAM_QA_PASS", Pass=all(upstream_qa$Pass)),
    data.table(Check_ID="REPRESENTATIVE_TERMS_33", Pass=nrow(terms)==33L && uniqueN(terms$Term_Key)==33L),
    data.table(Check_ID="TERMS_FORMAL_NON_C7", Pass=all(terms$Full_Significant) && all(terms$Display_Eligible) && all(terms$Collection!="C7")),
    data.table(Check_ID="V01_DISPLAY_200_240", Pass=nrow(v1$display) %between% c(200L,240L)),
    data.table(Check_ID="DISPLAY_FUNCTIONAL_MEMBERSHIP", Pass=all(v1$display$Gene_Node_ID %chin% v1$final_membership$Gene_Node_ID)),
    data.table(Check_ID="TERMS_CONNECTED", Pass=all(terms$Term_Node_ID %chin% v1$final_membership$Term_Node_ID)),
    data.table(Check_ID="TERM_OVERLAP_MEMBERSHIP_MATCH",
      Pass=v1$membership[,.(Membership_N=uniqueN(Gene)),by=.(Direction,Term_Key)][
        terms[,.(Direction,Term_Key,Overlap_N)],on=.(Direction,Term_Key)][,all(Membership_N==Overlap_N)]),
    data.table(Check_ID="FORMAL_C1_3416", Pass=directional[State=="C1" & Recurrent & ORA_Eligible,uniqueN(Gene)]==3416L),
    data.table(Check_ID="FORMAL_C4_16", Pass=directional[State=="C4" & Recurrent & ORA_Eligible,uniqueN(Gene)]==16L),
    data.table(Check_ID="DISPLAY_230", Pass=uniqueN(display$Gene)==230L),
    data.table(Check_ID="DISPLAY_C1_214", Pass=display[State=="C1",uniqueN(Gene)]==214L),
    data.table(Check_ID="DISPLAY_C4_16", Pass=display[State=="C4",uniqueN(Gene)]==16L),
    data.table(Check_ID="C4_FUNCTIONAL6_OTHER10",
      Pass=display[State=="C4" & Display_Assignment=="FUNCTIONAL_TERM",uniqueN(Gene)]==6L &&
        display[State=="C4" & Display_Assignment=="OTHER",uniqueN(Gene)]==10L),
    data.table(Check_ID="OTHER_ONLY_UNMAPPED_C4", Pass=all(display[Display_Assignment=="OTHER",State=="C4" & !Functional_Term_Mapped])),
    data.table(Check_ID="ALL_C4_FORMAL_DISPLAYED", Pass=setequal(directional[State=="C4" & Recurrent & ORA_Eligible,Gene],display[State=="C4",Gene])),
    data.table(Check_ID="DEFINITION_MARKER_PRESERVED", Pass="Definition_Overlap" %in% names(display) && !anyNA(display$Definition_Overlap))
  ))
  if(anyNA(qa$Pass) || !all(qa$Pass)) stop("Representative-input QA failed: ",paste(qa[!Pass,Check_ID],collapse="; "))
  fwrite(display,file.path(stage_dir,"mrna_c1c4_extreme_display_genes_v0.2_20260730.tsv"),sep="	",quote=TRUE,na="NA")
  fwrite(qa,file.path(stage_dir,"mrna_representative_input_machine_qa.tsv"),sep="	",quote=TRUE,na="NA")
  fwrite(data.table(Input=names(input_paths),Path=unname(input_paths)),
         file.path(stage_dir,"mrna_representative_input_sources.tsv"),sep="	",quote=TRUE,na="NA")
  writeLines(c(
    "Source: 37_freeze_stage3_mrna_c1c4_extreme_display_v0.1_20260730.R; 39_refine_stage3_mrna_c1c4_other_legacy_cytoscape_interface_v0.2_20260730.R",
    "Fixed 33 reviewed terms; source gene selection; source C4 merge/Other. Upstream P/FDR and 22-cancer recurrence are inherited.",
    "Output is the full-schema display input for public 13. Jaccard/Louvain, Cytoscape/legacy tables and figures are not generated.",
    "Two necessary v0.1 intermediate TSVs preserve the source37 fwrite/source39 fread boundary.",
    "These run-specific QA checks concern representative input only; no complete historical network QA is claimed."
  ),file.path(stage_dir,"README.txt"))
  writeLines(capture.output(sessionInfo()),file.path(stage_dir,"sessionInfo.txt"))
  if(!file.rename(stage_dir,final_dir)) stop("Cannot commit output stage.")
  cat("COMPLETED: ",final_dir,"
",sep="")
}
tryCatch(main(),error=function(e) { message(conditionMessage(e)); quit(status=1L) })
}
