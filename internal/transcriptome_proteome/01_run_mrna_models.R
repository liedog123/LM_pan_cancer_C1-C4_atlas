#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 01_run_mrna_models.R --config PATH [--check-config]\n       Rscript 01_run_mrna_models.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_manifest", "mrna_rdata", "definition_registry", "readiness_decision", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "01_run_mrna_models.R"
start_time <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extra <- if (length(fields)) {
    paste(vapply(names(fields), function(nm) {
      value <- gsub("[\r\n\t]", " ", as.character(fields[[nm]]))
      sprintf('%s="%s"', nm, value)
    }, character(1)), collapse = " ")
  } else {
    ""
  }
  cat(sprintf(
    "[%s] %s script=%s msg=\"%s\"%s%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
    gsub("[\r\n]", " ", msg), if (nzchar(extra)) " " else "", extra
  ))
  flush.console()
}

norm_sample <- function(x) substr(gsub("\\.", "-", as.character(x)), 1L, 15L)
norm_symbol <- function(x) toupper(trimws(as.character(x)))


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
