# Purpose: boundary concordance analysis stages.
# Configuration: config/molecular.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Test the three adjacent-state comparisons using their original ranked-gene families.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/09_run_mrna_boundary_models.R
stage_transcriptome_proteome_09_run_mrna_boundary_models <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
required_files <- input_paths
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
script_name <- "09_run_mrna_boundary_models.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_011, environment())
`norm_sample` <- bind_stage_helper(shared_norm_sample_037, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())


# 三条边是AG/IG二维状态图中的局部原子边，不是full subtype OVR。
edge_specs <- list(
  E1 = list(target = "AG1_IG1", reference = "AG2_IG1",
            boundary = "FIXED_IG1_AG1_TO_AG2"),
  E2 = list(target = "AG2_IG1", reference = "AG2_IG2",
            boundary = "FIXED_AG2_IG1_TO_IG2"),
  E3 = list(target = "AG2_IG2", reference = "AG3_IG2",
            boundary = "FIXED_IG2_AG2_TO_AG3")
)

edge_registry <- rbindlist(lapply(names(edge_specs), function(id) {
  x <- edge_specs[[id]]
  data.table(
    Edge_ID = id, Design_Type = "ATOMIC_PAIRWISE_EDGE", Model_Variable = "AGIG",
    Target = x$target, Reference = x$reference, Boundary = x$boundary,
    BH_Family = "WITHIN_CANCER_X_EDGE", Min_Per_Group = 10L,
    Min_Total = 30L, Purity_Adjustment = FALSE,
    Figure_Eligibility = "PENDING_POST_FREEZE_GATE"
  )
}))
`align_matrix_to_patient` <- bind_stage_helper(shared_align_matrix_to_patient_017, environment())

run_edge <- function(expr, meta, feature, cancer, edge_id, spec) {
  keep <- meta$AGIG %chin% c(spec$reference, spec$target)
  y <- expr[, keep, drop = FALSE]
  meta0 <- meta[keep]
  if (!ncol(y)) stop(cancer, " ", edge_id, "无设计样本")

  valid <- is.finite(y)
  target_idx <- which(meta0$AGIG == spec$target)
  reference_idx <- which(meta0$AGIG == spec$reference)
  n_target <- if (length(target_idx)) rowSums(valid[, target_idx, drop = FALSE]) else rep(0L, nrow(y))
  n_reference <- if (length(reference_idx)) rowSums(valid[, reference_idx, drop = FALSE]) else rep(0L, nrow(y))
  n_used <- n_target + n_reference
  eligible <- n_target >= 10L & n_reference >= 10L & n_used >= 30L

  if (!length(target_idx) || !length(reference_idx)) {
    eligible[] <- FALSE
    failure_reason <- rep("TARGET_OR_REFERENCE_ABSENT", nrow(y))
  } else {
    failure_reason <- fifelse(eligible, NA_character_, "PAIR_SAMPLE_THRESHOLD_FAILED")
  }

  grp <- factor(meta0$AGIG, levels = c(spec$reference, spec$target))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- c(spec$reference, spec$target)
  contrast <- setNames(c(-1, 1), c(spec$reference, spec$target))
  fit <- lmFit(y, design)
  fit2 <- contrasts.fit(fit, contrasts = matrix(contrast, ncol = 1L))
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  estimate <- as.numeric(fit2$coefficients[, 1])
  se <- as.numeric(fit2$stdev.unscaled[, 1] * sqrt(fit2$s2.post))
  statistic <- as.numeric(fit2$t[, 1])
  p_value <- as.numeric(fit2$p.value[, 1])
  b_value <- as.numeric(fit2$lods[, 1])
  for (x in list(estimate, se, statistic, p_value, b_value)) invisible(x)
  estimate[!eligible] <- se[!eligible] <- statistic[!eligible] <- NA_real_
  p_value[!eligible] <- b_value[!eligible] <- NA_real_

  result <- cbind(
    data.table(Cancer = cancer, Layer = "MRNA"), feature,
    data.table(
      Edge_ID = edge_id, Design_Type = "ATOMIC_PAIRWISE_EDGE",
      Analysis_Set = "POST_FREEZE_FIG3_BOUNDARY_EXTENSION",
      Model_Type = "UNADJUSTED", Target = spec$target,
      Reference = spec$reference, Boundary = spec$boundary,
      N_Input = ncol(y), N_Used = n_used, N_Target = n_target,
      N_Reference = n_reference, Estimate = estimate, SE = se,
      Statistic = statistic, P_Value = p_value, FDR = NA_real_,
      AveExpr = rowMeans(y, na.rm = TRUE), B = b_value,
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
  if (is.null(mrna_list[[cancer]]$expression_data)) stop("缺少mRNA表达矩阵: ", cancer)
  raw <- as.matrix(mrna_list[[cancer]]$expression_data)
  storage.mode(raw) <- "double"
  if (any(raw < 0, na.rm = TRUE)) stop("mRNA出现负TPM: ", cancer)
  if (anyDuplicated(rownames(raw))) stop("mRNA gene symbol重复: ", cancer)

  aligned <- align_matrix_to_patient(raw, patient_cancer)
  keep_state <- aligned$meta$AGIG %chin% unique(c(
    vapply(edge_specs, `[[`, character(1), "target"),
    vapply(edge_specs, `[[`, character(1), "reference")
  ))
  aligned$meta <- aligned$meta[keep_state]
  aligned$mat <- aligned$mat[, keep_state, drop = FALSE]
  if (!ncol(aligned$mat)) stop("无边界原子组样本: ", cancer)

  log_expr_all <- log2(aligned$mat + 1)
  keep_gene <- rowMeans(aligned$mat >= 1, na.rm = TRUE) >= 0.20 &
    apply(log_expr_all, 1L, IQR, na.rm = TRUE) > 0
  expr <- log_expr_all[keep_gene, , drop = FALSE]
  feature <- data.table(
    Feature_ID = rownames(expr), Feature_Label = rownames(expr),
    Feature_Class = "GENE",
    Definition_Overlap = norm_symbol(rownames(expr)) %in% definition_symbols
  )

  result <- rbindlist(lapply(names(edge_specs), function(edge_id) {
    run_edge(expr, aligned$meta, feature, cancer, edge_id, edge_specs[[edge_id]])
  }), use.names = TRUE)
  basename_out <- sprintf("%s_mrna_boundary_models_v0.1_20260806.tsv.gz", cancer)
  tmp_file <- file.path(tmp_individual_dir, basename_out)
  final_file <- file.path(final_individual_dir, basename_out)
  fwrite(result, tmp_file, sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  if (!file.exists(tmp_file) || file.info(tmp_file)$size <= 0) stop("逐癌输出为空: ", tmp_file)

  coverage <- result[, .(
    N_Features = .N, N_Estimable = sum(Status == "ESTIMABLE"),
    N_FDR_005 = sum(FDR < 0.05, na.rm = TRUE),
    N_FDR_005_Effect_GT_1 = sum(FDR < 0.05 & abs(Estimate) > 1, na.rm = TRUE)
  ), by = .(Cancer, Layer, Edge_ID, Target, Reference, Boundary)]
  agig_counts <- table(aligned$meta$AGIG)
  edge_counts <- rbindlist(lapply(names(edge_specs), function(id) {
    spec <- edge_specs[[id]]
    data.table(
      Cancer = cancer, Edge_ID = id, Target = spec$target,
      Reference = spec$reference,
      N_Target_Samples = sum(aligned$meta$AGIG == spec$target),
      N_Reference_Samples = sum(aligned$meta$AGIG == spec$reference),
      Edge_Evaluable = sum(aligned$meta$AGIG == spec$target) >= 10L &&
        sum(aligned$meta$AGIG == spec$reference) >= 10L &&
        sum(aligned$meta$AGIG %chin% c(spec$target, spec$reference)) >= 30L
    )
  }))
  manifest <- data.table(
    Cancer = cancer, Layer = "MRNA", N_Source_Features = nrow(raw),
    N_Filtered_Features = nrow(expr), N_Matched_Patients = ncol(expr),
    Output = final_file, Temporary_Output = tmp_file
  )
  sample_manifest <- aligned$meta[, .(
    Cancer, Patient_ID, Sample, C_subtype, AGIG, AG, IG, Acid_z, IG_z
  )]
  list(coverage = coverage, edge_counts = edge_counts, manifest = manifest,
       sample_manifest = sample_manifest,
       definition_audit = data.table(
         Cancer = cancer, N_Filtered_Features = nrow(feature),
         N_Definition_Features_Retained = sum(feature$Definition_Overlap)
       ))
}

main <- function() {
  args <- list(mode = "formal", cancers = "ALL")
  final_dir <- file.path(out_dir,"01_within_cancer_models")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(final_dir) && length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    stop("输出目录已存在且非空，拒绝覆盖: ", final_dir)
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  tmp_individual_dir <- file.path(tmp_dir, "individual")
  final_individual_dir <- file.path(final_dir, "individual")
  dir.create(tmp_individual_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(tmp_dir, "qa"), recursive = TRUE, showWarnings = FALSE)

  paths <- as.list(setNames(unname(paths[c("patient_manifest","mrna_rdata","definition_registry","readiness_decision")]),c("patient","mrna","definition","gate")))
  if (any(!file.exists(unlist(paths)))) stop("缺少冻结输入")
  gate <- fread(paths$gate)
  if (nrow(gate) != 1L || gate$decision != "PASS_STAGE3_READINESS_AUTHORIZE_FORMAL_ANALYSIS") {
    stop("Stage 3 readiness未授权")
  }

  patient <- fread(paths$patient)
  required <- c("Patient_ID", "Sample", "Cancer", "C_subtype", "AGIG", "AG", "IG")
  if (any(!required %in% names(patient))) stop("代表患者接口缺少AG/IG字段")
  patient[, Sample := norm_sample(Sample)]
  if (anyDuplicated(patient$Patient_ID) || anyDuplicated(patient$Sample)) stop("代表患者接口不唯一")
  canonical_map <- unique(patient[, .(C_subtype, AGIG)])
  expected_map <- data.table(
    C_subtype = c("C1", "C2", "C2", "C3", "C3", "C4"),
    AGIG = c("AG1_IG1", "AG1_IG2", "AG2_IG2", "AG2_IG1", "AG3_IG1", "AG3_IG2")
  )
  if (!setequal(paste(canonical_map$C_subtype, canonical_map$AGIG),
                paste(expected_map$C_subtype, expected_map$AGIG))) {
    stop("canonical C subtype–AGIG映射不符合冻结接口")
  }

  all_cancers <- sort(unique(patient$Cancer))
  if(!setequal(all_cancers,expected_cancers) || length(all_cancers)!=22L) stop("Complete 22-cancer patient scope required.")
  cancers <- if (args$mode == "formal" && identical(args$cancers, "ALL")) all_cancers else args$cancers
  if (any(!cancers %in% all_cancers)) stop("未知癌种: ", paste(setdiff(cancers, all_cancers), collapse = ","))
  definition <- fread(paths$definition)
  definition_symbols <- norm_symbol(definition$gene_symbol)
  mrna_env <- new.env(parent = emptyenv())
  load(paths$mrna, envir = mrna_env)
  if (!exists("pancancer_data_list", envir = mrna_env, inherits = FALSE)) stop("mRNA输入缺少pancancer_data_list")
  mrna_list <- mrna_env$pancancer_data_list

  log_event("START", "Fig3 mRNA boundary within-cancer models begin", list(
    mode = args$mode, environment = "R4.5.1", cancers = paste(cancers, collapse = ","),
    edges = paste(names(edge_specs), collapse = ","), output = final_dir,
    forbidden = "OVR,META,REML,LOO,CROSS_CANCER_BH"))

  fwrite(edge_registry, file.path(tmp_dir, "edge_design_registry.tsv"), sep = "\t", quote = TRUE)
  runs <- lapply(cancers, function(cancer) {
    ans <- run_cancer(cancer, mrna_list, patient, definition_symbols,
                      tmp_individual_dir, final_individual_dir)
    log_event("CHECKPOINT", "cancer boundary models complete", list(
      cancer = cancer, n_features = ans$manifest$N_Filtered_Features,
      n_patients = ans$manifest$N_Matched_Patients,
      evaluable_edges = sum(ans$edge_counts$Edge_Evaluable)))
    gc()
    ans
  })
  coverage <- rbindlist(lapply(runs, `[[`, "coverage"))
  edge_counts <- rbindlist(lapply(runs, `[[`, "edge_counts"))
  manifest <- rbindlist(lapply(runs, `[[`, "manifest"))
  sample_manifest <- unique(rbindlist(lapply(runs, `[[`, "sample_manifest")))
  definition_audit <- rbindlist(lapply(runs, `[[`, "definition_audit"))

  fwrite(coverage, file.path(tmp_dir, "model_coverage.tsv"), sep = "\t", quote = TRUE, na = "NA")
  fwrite(edge_counts, file.path(tmp_dir, "sample_coverage.tsv"), sep = "\t", quote = TRUE)
  fwrite(manifest[, !"Temporary_Output"], file.path(tmp_dir, "input_manifest.tsv"), sep = "\t", quote = TRUE)
  fwrite(sample_manifest, file.path(tmp_dir, "sample_manifest.tsv.gz"), sep = "\t", quote = TRUE, compress = "gzip")
  fwrite(definition_audit, file.path(tmp_dir, "definition_retention_audit.tsv"), sep = "\t", quote = TRUE)

  forbidden_patterns <- c("vs_REST", "OVR", "CONTINUOUS", "OMNIBUS")
  qa <- data.table(
    Check = c("EXACTLY_THREE_ATOMIC_EDGES", "EXACT_AGIG_ENDPOINTS", "NO_OVR_CONTINUOUS_OMNIBUS",
              "NO_PURITY_MODEL_OR_FIELD", "ONE_NONEMPTY_FILE_PER_CANCER", "ESTIMABLE_SE_FINITE",
              "FDR_RANGE_VALID", "WITHIN_CANCER_EDGE_BH_FIXED", "DEFINITION_FLAG_PRESERVED"),
    Pass = c(
      nrow(edge_registry) == 3L && setequal(edge_registry$Edge_ID, c("E1", "E2", "E3")),
      identical(edge_registry$Target, c("AG1_IG1", "AG2_IG1", "AG2_IG2")) &&
        identical(edge_registry$Reference, c("AG2_IG1", "AG2_IG2", "AG3_IG2")),
      !any(vapply(forbidden_patterns, function(p) any(grepl(p, edge_registry$Edge_ID, ignore.case = TRUE)), logical(1))),
      !any(edge_registry$Purity_Adjustment) && !any(c("purity_z", "TumorPurity") %in% names(sample_manifest)),
      nrow(manifest) == length(cancers) && all(file.exists(manifest$Temporary_Output)) &&
        all(file.info(manifest$Temporary_Output)$size > 0),
      all(vapply(manifest$Temporary_Output, function(path) {
        x <- fread(path, select = c("Status", "Estimate", "SE"))
        x <- x[Status == "ESTIMABLE" & is.finite(Estimate)]
        !nrow(x) || all(is.finite(x$SE))
      }, logical(1))),
      all(vapply(manifest$Temporary_Output, function(path) {
        x <- fread(path, select = "FDR")$FDR
        all(is.na(x) | (x >= 0 & x <= 1))
      }, logical(1))),
      all(edge_registry$BH_Family == "WITHIN_CANCER_X_EDGE"),
      all(definition_audit$N_Definition_Features_Retained >= 0L)
    )
  )
  fwrite(qa, file.path(tmp_dir, "qa", "machine_qa.tsv"), sep = "\t", quote = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(tmp_dir, "sessionInfo.txt"))
  if (!all(qa$Pass)) stop("Fig3 mRNA边界逐癌模型机器QA失败")

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) stop("无法提交Fig3 mRNA边界模型目录")
  log_event("DONE", "Fig3 mRNA boundary models complete", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    n_cancers = length(cancers), evaluable_cancer_edges = sum(edge_counts$Edge_Evaluable),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)), output = final_dir))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)))
  quit(status = 1L)
})
}

# Purpose: Test the three adjacent-state comparisons using their original ranked-gene families.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/10_run_e2_full_rank_gsea.R
stage_transcriptome_proteome_10_run_e2_full_rank_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("model_dir must be readable/searchable.")
required_files <- c(input_paths[names(input_paths)!="model_dir"],file.path(paths["model_dir"],"input_manifest.tsv"),file.path(paths["model_dir"],"individual",paste0(expected_cancers,"_mrna_boundary_models_v0.1_20260806.tsv.gz")))
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
  library(fgsea)
})

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
script_name <- "10_run_e2_full_rank_gsea.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_011, environment())
`fail` <- bind_stage_helper(shared_fail_018, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())
`seed_from_keys` <- bind_stage_helper(shared_seed_from_keys_001, environment())
`read_gmt` <- bind_stage_helper(shared_read_gmt_023, environment())

run_fgsea <- function(ranks, registry, membership, cancer, sensitivity, collection) {
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  if (anyDuplicated(names(ranks))) fail(paste(cancer, sensitivity, "rank gene重复"))
  present <- membership[Gene %chin% names(ranks)]
  sizes <- present[, .(Mapped_Size = uniqueN(Gene)), by = Set_ID]
  eligible <- sizes[Mapped_Size >= 10L & Mapped_Size <= 500L]
  meta <- merge(registry, eligible, by = "Set_ID")
  pathways <- split(present[Set_ID %chin% meta$Set_ID, Gene],
                    present[Set_ID %chin% meta$Set_ID, Set_ID])
  pathways <- lapply(pathways, unique)
  if (!length(pathways)) return(data.table())

  set.seed(seed_from_keys(cancer, sensitivity, collection))
  fg <- suppressWarnings(as.data.table(fgseaMultilevel(
    pathways = pathways, stats = ranks, sampleSize = 101L,
    minSize = 10L, maxSize = 500L, eps = 0, scoreType = "std",
    nproc = 1L, nPermSimple = 10000L
  )))
  setnames(fg, "pathway", "Set_ID")
  fg[, GSEA_Method := "fgseaMultilevel"]

  failed <- fg[!is.finite(pval) | !is.finite(NES) | !is.finite(ES), Set_ID]
  if (length(failed)) {
    set.seed(seed_from_keys(cancer, sensitivity, collection, "fallback"))
    fb <- suppressWarnings(as.data.table(fgseaSimple(
      pathways = pathways[failed], stats = ranks, nperm = 100000L,
      minSize = 10L, maxSize = 500L, scoreType = "std", nproc = 1L
    )))
    setnames(fb, "pathway", "Set_ID")
    for (id in failed) {
      i <- which(fg$Set_ID == id)
      z <- fb[Set_ID == id]
      if (length(i) != 1L || nrow(z) != 1L) {
        fail(paste("GSEA fallback无法唯一匹配", cancer, sensitivity, id))
      }
      if (!is.finite(z$pval) || !is.finite(z$NES) || !is.finite(z$ES)) {
        if (!is.finite(fg$ES[i]) || fg$ES[i] == 0) {
          fail(paste("GSEA方向fallback缺少有限ES", cancer, sensitivity, id))
        }
        score_type <- ifelse(fg$ES[i] > 0, "pos", "neg")
        set.seed(seed_from_keys(cancer, sensitivity, collection, id, score_type))
        z <- suppressWarnings(as.data.table(fgseaSimple(
          pathways = pathways[id], stats = ranks, nperm = 100000L,
          minSize = 10L, maxSize = 500L, scoreType = score_type, nproc = 1L
        )))
        setnames(z, "pathway", "Set_ID")
        if (nrow(z) != 1L || !is.finite(z$pval) ||
            !is.finite(z$NES) || !is.finite(z$ES)) {
          fail(paste("GSEA方向fallback仍失败", cancer, sensitivity, id))
        }
        method_label <- paste0("fgseaSimple_nperm100000_", score_type,
                               "_directional_fallback")
      } else {
        method_label <- "fgseaSimple_nperm100000_fallback"
      }
      fg$pval[i] <- z$pval; fg$ES[i] <- z$ES; fg$NES[i] <- z$NES
      fg$size[i] <- z$size; fg$leadingEdge[[i]] <- z$leadingEdge[[1L]]
      fg$log2err[i] <- NA_real_; fg$GSEA_Method[i] <- method_label
    }
  }
  if (any(!is.finite(fg$pval)) || any(!is.finite(fg$NES)) || any(!is.finite(fg$ES))) {
    fail(paste("GSEA仍有非有限结果", cancer, sensitivity, collection))
  }
  fg[, Leading_Edge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  fg[, leadingEdge := NULL]
  fg[, FDR := p.adjust(pval, method = "BH")]
  fg[, Direction := fifelse(NES > 0, "C3_LOCAL_HIGH", "C2_LOCAL_HIGH")]
  fg <- merge(fg, meta, by = "Set_ID", all.x = TRUE)
  fg[, `:=`(Cancer = cancer, Edge_ID = "E2", Sensitivity = sensitivity,
            Collection = collection, Target = "AG2_IG1", Reference = "AG2_IG2")]
  setnames(fg, c("pval", "log2err", "size"), c("P_Value", "Log2_Error", "Set_Size_Used"))
  fg[]
}

main <- function() {
  final_dir <- file.path(out_dir,"results")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(final_dir) && length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    fail(paste("输出目录已存在且非空，拒绝覆盖", final_dir))
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(file.path(tmp_dir, "individual"), recursive = TRUE, showWarnings = FALSE)

  model_dir <- unname(paths["model_dir"])
  paths <- as.list(setNames(unname(paths[c("definition_registry","immune_registry","immune_membership","lm12_membership","hallmark_gmt")]),c("definition","immune_registry","immune_long","lm12_long","hallmark")))
  paths$model_manifest <- file.path(model_dir,"input_manifest.tsv")
  if (any(!file.exists(unlist(paths)))) fail("缺少E2 GSEA冻结输入")

  manifest <- fread(paths$model_manifest)
  # 只按完整癌种及固定basename重绑路径，不读取历史绝对路径。
  if(!all(c("Cancer","Output") %in% names(manifest)) || nrow(manifest)!=22L || anyNA(manifest$Cancer) || anyDuplicated(manifest$Cancer) || !setequal(manifest$Cancer,expected_cancers)) fail("Complete unique 22-cancer model manifest required.")
  expected_names <- paste0(manifest$Cancer,"_mrna_boundary_models_v0.1_20260806.tsv.gz")
  if(anyNA(manifest$Output) || !identical(basename(manifest$Output),expected_names)) fail("Model manifest basenames do not match Cancer.")
  manifest[, Output := file.path(model_dir,"individual",expected_names)]
  definition <- norm_symbol(fread(paths$definition)$gene_symbol)
  immune_registry <- fread(paths$immune_registry)[, .(
    Set_ID = set_id, Display_Name = display_name, Collection = "IMMUNE119_PRIMARY",
    Functional_Domain = functional_domain, Evidence_Role = evidence_role,
    Source_Gene_Count = source_gene_count
  )]
  immune_membership <- fread(paths$immune_long)[, .(Set_ID = set_id, Gene = norm_symbol(gene))]
  if (uniqueN(immune_registry$Set_ID) != 119L) fail("冻结免疫背景不是119条")

  lm12 <- fread(paths$lm12_long)
  lm12_registry <- lm12[, .(
    Display_Name = pathway[1],
    Collection = "LM12_DEFINITION_CONTROL", Functional_Domain = functional_group[1],
    Evidence_Role = "DEFINITION_POSITIVE_CONTROL",
    Source_Gene_Count = uniqueN(gene)
  ), by = .(Set_ID = feature_id)]
  lm12_membership <- unique(lm12[, .(Set_ID = feature_id, Gene = norm_symbol(gene))])
  if (uniqueN(lm12_registry$Set_ID) != 12L) fail("LM12定义对照不是12条")

  hallmark <- read_gmt(paths$hallmark, "HALLMARK_CONTEXT", "GENERAL_CONTEXT_ONLY")
  registry <- rbindlist(list(immune_registry, lm12_registry, hallmark$registry), use.names = TRUE)
  membership <- rbindlist(list(immune_membership, lm12_membership, hallmark$membership), use.names = TRUE)
  membership <- unique(membership[nzchar(Gene)])

  log_event("START", "E2 full-rank GSEA begin", list(
    environment = "R4.5.1", model = "reuse_existing_E2_moderated_t",
    cancers = nrow(manifest), collections = paste(unique(registry$Collection), collapse = ","),
    sensitivities = "FULL_RANK,NON_DEFINITION_RANK", output = final_dir,
    prohibited = "rerun_limma,purity,meta,LOO,cross_cancer_BH,C7"))

  all_results <- list()
  coverage <- list()
  out_manifest <- list()
  for (i in seq_len(nrow(manifest))) {
    cancer <- manifest$Cancer[i]
    x <- fread(manifest$Output[i])
    x <- x[Edge_ID == "E2" & Status == "ESTIMABLE" & is.finite(Statistic)]
    if (!nrow(x)) next
    x[, Gene := norm_symbol(Feature_ID)]
    if (anyDuplicated(x$Gene)) fail(paste(cancer, "E2 rank gene重复"))
    full_rank <- setNames(x$Statistic, x$Gene)
    rank_list <- list(
      FULL_RANK = full_rank,
      NON_DEFINITION_RANK = full_rank[!names(full_rank) %chin% definition]
    )
    cancer_rows <- list()
    for (sensitivity in names(rank_list)) {
      for (collection in unique(registry$Collection)) {
        reg0 <- registry[Collection == collection]
        mem0 <- membership[Set_ID %chin% reg0$Set_ID]
        z <- run_fgsea(rank_list[[sensitivity]], reg0, mem0, cancer,
                       sensitivity, collection)
        cancer_rows[[paste(sensitivity, collection)]] <- z
        coverage[[paste(cancer, sensitivity, collection)]] <- data.table(
          Cancer = cancer, Sensitivity = sensitivity, Collection = collection,
          Rank_N = length(rank_list[[sensitivity]]), N_Sets_Evaluable = nrow(z),
          N_FDR_005 = if (nrow(z)) sum(z$FDR < 0.05) else 0L
        )
      }
    }
    cancer_result <- rbindlist(cancer_rows, use.names = TRUE, fill = TRUE)
    out_path <- file.path(tmp_dir, "individual",
      paste0(cancer, "_e2_full_rank_gsea_v0.1_20260806.tsv.gz"))
    fwrite(cancer_result, out_path, sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
    if (!file.exists(out_path) || file.info(out_path)$size <= 0) fail(paste(cancer, "逐癌GSEA输出为空"))
    out_manifest[[cancer]] <- data.table(
      Cancer = cancer,
      Output = file.path(final_dir, "individual", basename(out_path)),
      N_Rows = nrow(cancer_result)
    )
    all_results[[cancer]] <- cancer_result
    log_event("CHECKPOINT", "cancer E2 GSEA complete", list(
      cancer = cancer, full_rank_n = length(full_rank),
      nondefinition_rank_n = length(rank_list$NON_DEFINITION_RANK),
      result_rows = nrow(cancer_result), significant_rows = sum(cancer_result$FDR < 0.05)))
  }
  gsea <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  if (!nrow(gsea)) fail("所有癌种E2 GSEA均为空")
  fwrite(gsea, file.path(tmp_dir, "e2_gsea_all_results.tsv.gz"), sep = "\t",
         quote = TRUE, na = "NA", compress = "gzip")
  fwrite(rbindlist(coverage), file.path(tmp_dir, "e2_gsea_coverage.tsv"), sep = "\t", quote = TRUE)

  # 跨癌只做等权方向与显著票数，不合并P值。
  vote <- gsea[, .(
    Display_Name = Display_Name[1], Functional_Domain = Functional_Domain[1],
    Evidence_Role = Evidence_Role[1], N_Evaluable = uniqueN(Cancer),
    N_C3_Local_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES > 0]),
    N_C2_Local_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES < 0]),
    C3_Local_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES > 0])), collapse = ";"),
    C2_Local_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES < 0])), collapse = ";"),
    N_C3_Local_High_NES = sum(NES > 0), N_C2_Local_High_NES = sum(NES < 0),
    Median_NES = median(NES), Median_Abs_NES = median(abs(NES)),
    Leading_Edge_Union = paste(sort(unique(unlist(strsplit(Leading_Edge, ";", fixed = TRUE)))), collapse = ";")
  ), by = .(Sensitivity, Collection, Set_ID)]
  vote[, `:=`(
    Min_Recurrence = pmax(3L, ceiling(0.30 * N_Evaluable)),
    Min_Strong_Support = pmax(4L, ceiling(0.40 * N_Evaluable)),
    Max_Opposite_Allowed = floor(0.10 * N_Evaluable),
    C3_NES_Direction_Fraction = N_C3_Local_High_NES / N_Evaluable,
    C2_NES_Direction_Fraction = N_C2_Local_High_NES / N_Evaluable
  )]
  vote[, `:=`(
    Recurrent_C3_Local_High = N_C3_Local_High_Sig >= Min_Recurrence,
    Recurrent_C2_Local_High = N_C2_Local_High_Sig >= Min_Recurrence,
    Strong_C3_Local_High = N_C3_Local_High_Sig >= Min_Strong_Support &
      N_C2_Local_High_Sig <= Max_Opposite_Allowed,
    Strong_C2_Local_High = N_C2_Local_High_Sig >= Min_Strong_Support &
      N_C3_Local_High_Sig <= Max_Opposite_Allowed
  )]
  vote[, Support_Class := fifelse(
    Strong_C3_Local_High, "STRONG_C3_LOCAL_HIGH",
    fifelse(Strong_C2_Local_High, "STRONG_C2_LOCAL_HIGH",
      fifelse(Recurrent_C3_Local_High, "RECURRENT_C3_LOCAL_HIGH",
        fifelse(Recurrent_C2_Local_High, "RECURRENT_C2_LOCAL_HIGH", "NOT_RECURRENT"))))]
  setorder(vote, Sensitivity, Collection, -N_C3_Local_High_Sig,
           -N_C2_Local_High_Sig, -Median_Abs_NES, Set_ID)
  fwrite(vote, file.path(tmp_dir, "e2_gsea_cross_cancer_votes.tsv"), sep = "\t", quote = TRUE)

  candidates <- vote[
    Sensitivity == "NON_DEFINITION_RANK" & Collection == "IMMUNE119_PRIMARY" &
      Support_Class %chin% c("STRONG_C3_LOCAL_HIGH", "STRONG_C2_LOCAL_HIGH")
  ]
  candidates[, `:=`(Manual_Redundancy_Review = "PENDING",
                    Fig3B_Eligibility = "PENDING_MANUAL_GATE")]
  fwrite(candidates, file.path(tmp_dir, "e2_gsea_nondefinition_immune_candidates.tsv"),
         sep = "\t", quote = TRUE)

  # FULL与去定义层按同一集合并排，便于识别定义驱动。
  comparison <- dcast(vote[Collection == "IMMUNE119_PRIMARY"],
    Collection + Set_ID + Display_Name + Functional_Domain ~ Sensitivity,
    value.var = c("N_C3_Local_High_Sig", "N_C2_Local_High_Sig", "Support_Class"))
  fwrite(comparison, file.path(tmp_dir, "e2_gsea_definition_sensitivity_comparison.tsv"),
         sep = "\t", quote = TRUE)

  qa <- data.table(
    Check = c("E2_ONLY", "EXPECTED_COLLECTIONS", "TWO_SENSITIVITIES",
              "IMMUNE119_EXACT", "LM12_EXACT", "NO_C7", "NO_META_FIELDS",
              "FDR_RANGE_VALID", "NONDEFINITION_RANK_EXCLUDES_DEFINITION",
              "NINETEEN_EVALUABLE_CANCERS", "OUTPUT_NONEMPTY"),
    Pass = c(
      all(gsea$Edge_ID == "E2"),
      setequal(unique(gsea$Collection), c("IMMUNE119_PRIMARY", "LM12_DEFINITION_CONTROL", "HALLMARK_CONTEXT")),
      setequal(unique(gsea$Sensitivity), c("FULL_RANK", "NON_DEFINITION_RANK")),
      uniqueN(immune_registry$Set_ID) == 119L,
      uniqueN(lm12_registry$Set_ID) == 12L,
      !any(grepl("C7", gsea$Collection)),
      !any(c("Meta_FDR", "I2", "Tau2", "LOO") %in% names(gsea)),
      all(is.finite(gsea$FDR) & gsea$FDR >= 0 & gsea$FDR <= 1),
      length(intersect(
        unique(unlist(strsplit(
          gsea[Sensitivity == "NON_DEFINITION_RANK", Leading_Edge],
          ";", fixed = TRUE
        ))),
        definition
      )) == 0L,
      uniqueN(gsea$Cancer) == 19L,
      nrow(gsea) > 0L && nrow(vote) > 0L
    )
  )
  fwrite(qa, file.path(tmp_dir, "machine_qa.tsv"), sep = "\t", quote = TRUE)
  if (!all(qa$Pass)) fail("E2 full-rank GSEA机器QA失败")

  summary <- vote[, .(
    N_Sets = .N,
    N_Recurrent_C3 = sum(Recurrent_C3_Local_High),
    N_Recurrent_C2 = sum(Recurrent_C2_Local_High),
    N_Strong_C3 = sum(Strong_C3_Local_High),
    N_Strong_C2 = sum(Strong_C2_Local_High)
  ), by = .(Sensitivity, Collection)]
  fwrite(summary, file.path(tmp_dir, "e2_gsea_summary.tsv"), sep = "\t", quote = TRUE)
  fwrite(rbindlist(out_manifest), file.path(tmp_dir, "output_manifest.tsv"), sep = "\t", quote = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(tmp_dir, "sessionInfo.txt"))
  report <- c(
    "# E2 full-rank GSEA执行报告 V0.1", "", "日期：2026-08-06", "",
    "状态：`MACHINE_QA_PASS / MANUAL_GATE_PENDING / STAGE3_UNCHANGED`", "",
    "## 方法边界", "",
    "只复用既有E2逐癌moderated t，正NES表示AG2_IG1/C3局部侧高，负NES表示AG2_IG2/C2局部侧高。每癌种独立GSEA；跨癌仅做等权方向与显著票数，不运行meta、LOO或跨癌BH。", "",
    "FULL_RANK为描述层；NON_DEFINITION_RANK同时从排名和集合中去除全部Acid/LM12 definition genes。只有去定义后的冻结免疫119层可以进入新增免疫模块人工审核。", "",
    "## 数量摘要", "",
    paste(apply(summary, 1L, function(x) sprintf(
      "%s × %s：%s个可评价集合；C3侧复现/强支持%s/%s，C2侧复现/强支持%s/%s。",
      x[["Sensitivity"]], x[["Collection"]], x[["N_Sets"]],
      x[["N_Recurrent_C3"]], x[["N_Strong_C3"]],
      x[["N_Recurrent_C2"]], x[["N_Strong_C2"]]
    )), collapse = "\n\n"), "",
    "## 下一步", "",
    "人工审核NON_DEFINITION_RANK免疫强支持term的语义冗余和功能域数量；未达到两个独立免疫域时按合同停止Fig3B晋级。"
  )
  writeLines(report, file.path(tmp_dir, "execution_report.md"))

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) fail("无法提交E2 full-rank GSEA结果目录")
  log_event("DONE", "E2 full-rank GSEA complete", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    cancers = uniqueN(gsea$Cancer), gsea_rows = nrow(gsea),
    nondefinition_immune_strong_candidates = nrow(candidates),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)), output = final_dir,
    manual_gate = "PENDING"))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)))
  quit(status = 1L)
})
}

# Purpose: Test the three adjacent-state comparisons using their original ranked-gene families.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/11_run_e1_e3_full_rank_gsea.R
stage_transcriptome_proteome_11_run_e1_e3_full_rank_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("model_dir must be readable/searchable.")
required_files <- c(input_paths[names(input_paths)!="model_dir"],file.path(paths["model_dir"],"input_manifest.tsv"),file.path(paths["model_dir"],"individual",paste0(expected_cancers,"_mrna_boundary_models_v0.1_20260806.tsv.gz")))
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
  library(fgsea)
})

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
script_name <- "11_run_e1_e3_full_rank_gsea.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_011, environment())
`fail` <- bind_stage_helper(shared_fail_018, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())
`seed_from_keys` <- bind_stage_helper(shared_seed_from_keys_001, environment())
`read_gmt` <- bind_stage_helper(shared_read_gmt_023, environment())

run_fgsea <- function(ranks, registry, membership, cancer, sensitivity, collection,
                      edge_meta) {
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  if (anyDuplicated(names(ranks))) fail(paste(cancer, sensitivity, "rank gene重复"))
  present <- membership[Gene %chin% names(ranks)]
  sizes <- present[, .(Mapped_Size = uniqueN(Gene)), by = Set_ID]
  eligible <- sizes[Mapped_Size >= 10L & Mapped_Size <= 500L]
  meta <- merge(registry, eligible, by = "Set_ID")
  pathways <- split(present[Set_ID %chin% meta$Set_ID, Gene],
                    present[Set_ID %chin% meta$Set_ID, Set_ID])
  pathways <- lapply(pathways, unique)
  if (!length(pathways)) return(data.table())

  set.seed(seed_from_keys(cancer, sensitivity, collection))
  fg <- suppressWarnings(as.data.table(fgseaMultilevel(
    pathways = pathways, stats = ranks, sampleSize = 101L,
    minSize = 10L, maxSize = 500L, eps = 0, scoreType = "std",
    nproc = 1L, nPermSimple = 10000L
  )))
  setnames(fg, "pathway", "Set_ID")
  fg[, GSEA_Method := "fgseaMultilevel"]

  failed <- fg[!is.finite(pval) | !is.finite(NES) | !is.finite(ES), Set_ID]
  if (length(failed)) {
    set.seed(seed_from_keys(cancer, sensitivity, collection, "fallback"))
    fb <- suppressWarnings(as.data.table(fgseaSimple(
      pathways = pathways[failed], stats = ranks, nperm = 100000L,
      minSize = 10L, maxSize = 500L, scoreType = "std", nproc = 1L
    )))
    setnames(fb, "pathway", "Set_ID")
    for (id in failed) {
      i <- which(fg$Set_ID == id)
      z <- fb[Set_ID == id]
      if (length(i) != 1L || nrow(z) != 1L) {
        fail(paste("GSEA fallback无法唯一匹配", cancer, sensitivity, id))
      }
      if (!is.finite(z$pval) || !is.finite(z$NES) || !is.finite(z$ES)) {
        if (!is.finite(fg$ES[i]) || fg$ES[i] == 0) {
          fail(paste("GSEA方向fallback缺少有限ES", cancer, sensitivity, id))
        }
        score_type <- ifelse(fg$ES[i] > 0, "pos", "neg")
        set.seed(seed_from_keys(cancer, sensitivity, collection, id, score_type))
        z <- suppressWarnings(as.data.table(fgseaSimple(
          pathways = pathways[id], stats = ranks, nperm = 100000L,
          minSize = 10L, maxSize = 500L, scoreType = score_type, nproc = 1L
        )))
        setnames(z, "pathway", "Set_ID")
        if (nrow(z) != 1L || !is.finite(z$pval) ||
            !is.finite(z$NES) || !is.finite(z$ES)) {
          fail(paste("GSEA方向fallback仍失败", cancer, sensitivity, id))
        }
        method_label <- paste0("fgseaSimple_nperm100000_", score_type,
                               "_directional_fallback")
      } else {
        method_label <- "fgseaSimple_nperm100000_fallback"
      }
      fg$pval[i] <- z$pval; fg$ES[i] <- z$ES; fg$NES[i] <- z$NES
      fg$size[i] <- z$size; fg$leadingEdge[[i]] <- z$leadingEdge[[1L]]
      fg$log2err[i] <- NA_real_; fg$GSEA_Method[i] <- method_label
    }
  }
  if (any(!is.finite(fg$pval)) || any(!is.finite(fg$NES)) || any(!is.finite(fg$ES))) {
    fail(paste("GSEA仍有非有限结果", cancer, sensitivity, collection))
  }
  fg[, Leading_Edge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  fg[, leadingEdge := NULL]
  fg[, FDR := p.adjust(pval, method = "BH")]
  fg[, Direction := fifelse(NES > 0, edge_meta$Positive_Direction,
                            edge_meta$Negative_Direction)]
  fg <- merge(fg, meta, by = "Set_ID", all.x = TRUE)
  fg[, `:=`(Cancer = cancer, Edge_ID = edge_meta$Edge_ID,
            Sensitivity = sensitivity, Collection = collection,
            Target = edge_meta$Target, Reference = edge_meta$Reference,
            Target_State = edge_meta$Target_State,
            Reference_State = edge_meta$Reference_State)]
  setnames(fg, c("pval", "log2err", "size"), c("P_Value", "Log2_Error", "Set_Size_Used"))
  fg[]
}

main <- function() {
  edge_id <- unname(cfg["edge_id"])
  edge_registry <- data.table(
    Edge_ID = c("E1", "E3"),
    Target = c("AG1_IG1", "AG2_IG2"),
    Reference = c("AG2_IG1", "AG3_IG2"),
    Target_State = c("C1", "C2"),
    Reference_State = c("C3", "C4"),
    Positive_Direction = c("C1_LOCAL_HIGH", "C2_LOCAL_HIGH"),
    Negative_Direction = c("C3_LOCAL_HIGH", "C4_LOCAL_HIGH"),
    Expected_Cancers = c(20L, 21L)
  )
  edge_meta <- edge_registry[Edge_ID == edge_id]
  final_dir <- file.path(out_dir,"results","edges",edge_id)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(final_dir) && length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    fail(paste("输出目录已存在且非空，拒绝覆盖", final_dir))
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(file.path(tmp_dir, "individual"), recursive = TRUE, showWarnings = FALSE)

  model_dir <- unname(paths["model_dir"])
  paths <- as.list(setNames(unname(paths[c("definition_registry","immune_registry","immune_membership","lm12_membership","hallmark_gmt")]),c("definition","immune_registry","immune_long","lm12_long","hallmark")))
  paths$model_manifest <- file.path(model_dir,"input_manifest.tsv")
  if (any(!file.exists(unlist(paths)))) fail(paste("缺少", edge_id, "GSEA冻结输入"))

  manifest <- fread(paths$model_manifest)
  # 只按完整癌种及固定basename重绑路径，不读取历史绝对路径。
  if(!all(c("Cancer","Output") %in% names(manifest)) || nrow(manifest)!=22L || anyNA(manifest$Cancer) || anyDuplicated(manifest$Cancer) || !setequal(manifest$Cancer,expected_cancers)) fail("Complete unique 22-cancer model manifest required.")
  expected_names <- paste0(manifest$Cancer,"_mrna_boundary_models_v0.1_20260806.tsv.gz")
  if(anyNA(manifest$Output) || !identical(basename(manifest$Output),expected_names)) fail("Model manifest basenames do not match Cancer.")
  manifest[, Output := file.path(model_dir,"individual",expected_names)]
  definition <- norm_symbol(fread(paths$definition)$gene_symbol)
  immune_registry <- fread(paths$immune_registry)[, .(
    Set_ID = set_id, Display_Name = display_name, Collection = "IMMUNE119_PRIMARY",
    Functional_Domain = functional_domain, Evidence_Role = evidence_role,
    Source_Gene_Count = source_gene_count
  )]
  immune_membership <- fread(paths$immune_long)[, .(Set_ID = set_id, Gene = norm_symbol(gene))]
  if (uniqueN(immune_registry$Set_ID) != 119L) fail("冻结免疫背景不是119条")

  lm12 <- fread(paths$lm12_long)
  lm12_registry <- lm12[, .(
    Display_Name = pathway[1],
    Collection = "LM12_DEFINITION_CONTROL", Functional_Domain = functional_group[1],
    Evidence_Role = "DEFINITION_POSITIVE_CONTROL",
    Source_Gene_Count = uniqueN(gene)
  ), by = .(Set_ID = feature_id)]
  lm12_membership <- unique(lm12[, .(Set_ID = feature_id, Gene = norm_symbol(gene))])
  if (uniqueN(lm12_registry$Set_ID) != 12L) fail("LM12定义对照不是12条")

  hallmark <- read_gmt(paths$hallmark, "HALLMARK_CONTEXT", "GENERAL_CONTEXT_ONLY")
  registry <- rbindlist(list(immune_registry, lm12_registry, hallmark$registry), use.names = TRUE)
  membership <- rbindlist(list(immune_membership, lm12_membership, hallmark$membership), use.names = TRUE)
  membership <- unique(membership[nzchar(Gene)])

  log_event("START", "single-edge full-rank GSEA begin", list(
    edge = edge_id, environment = "R4.5.1",
    model = paste0("reuse_existing_", edge_id, "_moderated_t"),
    cancers = nrow(manifest), collections = paste(unique(registry$Collection), collapse = ","),
    sensitivities = "FULL_RANK,NON_DEFINITION_RANK", output = final_dir,
    prohibited = "rerun_limma,purity,meta,LOO,cross_cancer_BH,C7"))

  all_results <- list()
  coverage <- list()
  out_manifest <- list()
  for (i in seq_len(nrow(manifest))) {
    cancer <- manifest$Cancer[i]
    x <- fread(manifest$Output[i])
    x <- x[Edge_ID == edge_id & Status == "ESTIMABLE" & is.finite(Statistic)]
    if (!nrow(x)) next
    x[, Gene := norm_symbol(Feature_ID)]
    if (anyDuplicated(x$Gene)) fail(paste(cancer, edge_id, "rank gene重复"))
    full_rank <- setNames(x$Statistic, x$Gene)
    rank_list <- list(
      FULL_RANK = full_rank,
      NON_DEFINITION_RANK = full_rank[!names(full_rank) %chin% definition]
    )
    cancer_rows <- list()
    for (sensitivity in names(rank_list)) {
      for (collection in unique(registry$Collection)) {
        reg0 <- registry[Collection == collection]
        mem0 <- membership[Set_ID %chin% reg0$Set_ID]
        z <- run_fgsea(rank_list[[sensitivity]], reg0, mem0, cancer,
                       sensitivity, collection, edge_meta)
        cancer_rows[[paste(sensitivity, collection)]] <- z
        coverage[[paste(cancer, sensitivity, collection)]] <- data.table(
          Cancer = cancer, Sensitivity = sensitivity, Collection = collection,
          Rank_N = length(rank_list[[sensitivity]]), N_Sets_Evaluable = nrow(z),
          N_FDR_005 = if (nrow(z)) sum(z$FDR < 0.05) else 0L
        )
      }
    }
    cancer_result <- rbindlist(cancer_rows, use.names = TRUE, fill = TRUE)
    out_path <- file.path(tmp_dir, "individual",
      paste0(cancer, "_", tolower(edge_id), "_full_rank_gsea_v0.1_20260806.tsv.gz"))
    fwrite(cancer_result, out_path, sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
    if (!file.exists(out_path) || file.info(out_path)$size <= 0) fail(paste(cancer, "逐癌GSEA输出为空"))
    out_manifest[[cancer]] <- data.table(
      Cancer = cancer,
      Output = file.path(final_dir, "individual", basename(out_path)),
      N_Rows = nrow(cancer_result)
    )
    all_results[[cancer]] <- cancer_result
    log_event("CHECKPOINT", "cancer edge GSEA complete", list(
      edge = edge_id, cancer = cancer, full_rank_n = length(full_rank),
      nondefinition_rank_n = length(rank_list$NON_DEFINITION_RANK),
      result_rows = nrow(cancer_result), significant_rows = sum(cancer_result$FDR < 0.05)))
  }
  gsea <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  if (!nrow(gsea)) fail(paste("所有癌种", edge_id, "GSEA均为空"))
  prefix <- tolower(edge_id)
  fwrite(gsea, file.path(tmp_dir, paste0(prefix, "_gsea_all_results.tsv.gz")), sep = "\t",
         quote = TRUE, na = "NA", compress = "gzip")
  fwrite(rbindlist(coverage), file.path(tmp_dir, paste0(prefix, "_gsea_coverage.tsv")), sep = "\t", quote = TRUE)

  # 跨癌只做等权方向与显著票数，不合并P值。
  vote <- gsea[, .(
    Display_Name = Display_Name[1], Functional_Domain = Functional_Domain[1],
    Evidence_Role = Evidence_Role[1], N_Evaluable = uniqueN(Cancer),
    N_Target_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES > 0]),
    N_Reference_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES < 0]),
    Target_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES > 0])), collapse = ";"),
    Reference_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES < 0])), collapse = ";"),
    N_Target_High_NES = sum(NES > 0), N_Reference_High_NES = sum(NES < 0),
    Median_NES = median(NES), Median_Abs_NES = median(abs(NES)),
    Leading_Edge_Union = paste(sort(unique(unlist(strsplit(Leading_Edge, ";", fixed = TRUE)))), collapse = ";")
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection, Set_ID)]
  vote[, `:=`(
    Min_Recurrence = pmax(3L, ceiling(0.30 * N_Evaluable)),
    Min_Strong_Support = pmax(4L, ceiling(0.40 * N_Evaluable)),
    Max_Opposite_Allowed = floor(0.10 * N_Evaluable),
    Target_NES_Direction_Fraction = N_Target_High_NES / N_Evaluable,
    Reference_NES_Direction_Fraction = N_Reference_High_NES / N_Evaluable
  )]
  vote[, `:=`(
    Recurrent_Target_High = N_Target_High_Sig >= Min_Recurrence,
    Recurrent_Reference_High = N_Reference_High_Sig >= Min_Recurrence,
    Strong_Target_High = N_Target_High_Sig >= Min_Strong_Support &
      N_Reference_High_Sig <= Max_Opposite_Allowed,
    Strong_Reference_High = N_Reference_High_Sig >= Min_Strong_Support &
      N_Target_High_Sig <= Max_Opposite_Allowed
  )]
  vote[, Support_Class := fifelse(
    Strong_Target_High, "STRONG_TARGET_HIGH",
    fifelse(Strong_Reference_High, "STRONG_REFERENCE_HIGH",
      fifelse(Recurrent_Target_High, "RECURRENT_TARGET_HIGH",
        fifelse(Recurrent_Reference_High, "RECURRENT_REFERENCE_HIGH", "NOT_RECURRENT"))))]
  setorder(vote, Sensitivity, Collection, -N_Target_High_Sig,
           -N_Reference_High_Sig, -Median_Abs_NES, Set_ID)
  fwrite(vote, file.path(tmp_dir, paste0(prefix, "_gsea_cross_cancer_votes.tsv")), sep = "\t", quote = TRUE)

  candidates <- vote[
    Sensitivity == "NON_DEFINITION_RANK" & Collection == "IMMUNE119_PRIMARY" &
      Support_Class %chin% c("STRONG_TARGET_HIGH", "STRONG_REFERENCE_HIGH")
  ]
  candidates[, `:=`(Manual_Redundancy_Review = "PENDING",
                    Fig3B_Eligibility = "PENDING_MANUAL_GATE")]
  fwrite(candidates, file.path(tmp_dir, paste0(prefix, "_gsea_nondefinition_immune_candidates.tsv")),
         sep = "\t", quote = TRUE)

  # FULL与去定义层按同一集合并排，便于识别定义驱动。
  comparison <- dcast(vote[Collection == "IMMUNE119_PRIMARY"],
    Collection + Set_ID + Display_Name + Functional_Domain ~ Sensitivity,
    value.var = c("N_Target_High_Sig", "N_Reference_High_Sig", "Support_Class"))
  fwrite(comparison, file.path(tmp_dir, paste0(prefix, "_gsea_definition_sensitivity_comparison.tsv")),
         sep = "\t", quote = TRUE)

  qa <- data.table(
    Check = c("REQUESTED_EDGE_ONLY", "EXPECTED_COLLECTIONS", "TWO_SENSITIVITIES",
              "IMMUNE119_EXACT", "LM12_EXACT", "NO_C7", "NO_META_FIELDS",
              "FDR_RANGE_VALID", "NONDEFINITION_RANK_EXCLUDES_DEFINITION",
              "EXPECTED_EVALUABLE_CANCERS", "OUTPUT_NONEMPTY"),
    Pass = c(
      all(gsea$Edge_ID == edge_id),
      setequal(unique(gsea$Collection), c("IMMUNE119_PRIMARY", "LM12_DEFINITION_CONTROL", "HALLMARK_CONTEXT")),
      setequal(unique(gsea$Sensitivity), c("FULL_RANK", "NON_DEFINITION_RANK")),
      uniqueN(immune_registry$Set_ID) == 119L,
      uniqueN(lm12_registry$Set_ID) == 12L,
      !any(grepl("C7", gsea$Collection)),
      !any(c("Meta_FDR", "I2", "Tau2", "LOO") %in% names(gsea)),
      all(is.finite(gsea$FDR) & gsea$FDR >= 0 & gsea$FDR <= 1),
      length(intersect(
        unique(unlist(strsplit(
          gsea[Sensitivity == "NON_DEFINITION_RANK", Leading_Edge],
          ";", fixed = TRUE
        ))),
        definition
      )) == 0L,
      uniqueN(gsea$Cancer) == edge_meta$Expected_Cancers,
      nrow(gsea) > 0L && nrow(vote) > 0L
    )
  )
  fwrite(qa, file.path(tmp_dir, "machine_qa.tsv"), sep = "\t", quote = TRUE)
  if (!all(qa$Pass)) fail(paste(edge_id, "full-rank GSEA机器QA失败"))

  summary <- vote[, .(
    N_Sets = .N,
    N_Recurrent_Target = sum(Recurrent_Target_High),
    N_Recurrent_Reference = sum(Recurrent_Reference_High),
    N_Strong_Target = sum(Strong_Target_High),
    N_Strong_Reference = sum(Strong_Reference_High)
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection)]
  fwrite(summary, file.path(tmp_dir, paste0(prefix, "_gsea_summary.tsv")), sep = "\t", quote = TRUE)
  fwrite(rbindlist(out_manifest), file.path(tmp_dir, "output_manifest.tsv"), sep = "\t", quote = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(tmp_dir, "sessionInfo.txt"))
  report <- c(
    paste0("# ", edge_id, " full-rank GSEA执行报告 V0.1"), "", "日期：2026-08-06", "",
    "状态：`MACHINE_QA_PASS / MANUAL_GATE_PENDING / STAGE3_UNCHANGED`", "",
    "## 方法边界", "",
    sprintf("只复用既有%s逐癌moderated t；正NES表示%s/%s局部侧高，负NES表示%s/%s局部侧高。每癌种独立GSEA；跨癌仅做等权方向与显著票数，不运行meta、LOO或跨癌BH。", edge_id, edge_meta$Target, edge_meta$Target_State, edge_meta$Reference, edge_meta$Reference_State), "",
    "FULL_RANK为描述层；NON_DEFINITION_RANK同时从排名和集合中去除全部Acid/LM12 definition genes。只有去定义后的冻结免疫119层可以进入新增免疫模块人工审核。", "",
    "## 数量摘要", "",
    paste(apply(summary, 1L, function(x) sprintf(
      "%s × %s：%s个可评价集合；target侧复现/强支持%s/%s，reference侧复现/强支持%s/%s。",
      x[["Sensitivity"]], x[["Collection"]], x[["N_Sets"]],
      x[["N_Recurrent_Target"]], x[["N_Strong_Target"]],
      x[["N_Recurrent_Reference"]], x[["N_Strong_Reference"]]
    )), collapse = "\n\n"), "",
    "## 下一步", "",
    "与既有E2结果统一汇总后，人工审核三条边NON_DEFINITION_RANK强支持term的语义冗余、definition敏感性和Fig3B展示资格。"
  )
  writeLines(report, file.path(tmp_dir, "execution_report.md"))

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) fail(paste("无法提交", edge_id, "full-rank GSEA结果目录"))
  log_event("DONE", "single-edge full-rank GSEA complete", list(
    edge = edge_id,
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    cancers = uniqueN(gsea$Cancer), gsea_rows = nrow(gsea),
    nondefinition_immune_strong_candidates = nrow(candidates),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)), output = final_dir,
    manual_gate = "PENDING"))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)))
  quit(status = 1L)
})
}

# Purpose: Test the three adjacent-state comparisons using their original ranked-gene families.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/12_combine_three_edge_gsea.R
stage_transcriptome_proteome_12_combine_three_edge_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
if(any(!dir.exists(paths[c("e1_dir","e2_dir","e3_dir")])) || any(file.access(paths[c("e1_dir","e2_dir","e3_dir")],5L)!=0L)) stop("Edge directories must be readable/searchable.")
required_files <- c(paths["definition_registry"],unlist(lapply(1:3,function(j) file.path(paths[paste0("e",j,"_dir")],c(paste0("e",j,"_gsea_all_results.tsv.gz"),"machine_qa.tsv"))),use.names=FALSE))
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
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

script_name <- "12_combine_three_edge_gsea.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_011, environment())
`fail` <- bind_stage_helper(shared_fail_018, environment())
`norm_symbol` <- bind_stage_helper(shared_norm_symbol_039, environment())

main <- function() {
  final_dir <- file.path(out_dir,"results","combined")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(final_dir) && length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    fail(paste("三边汇总结果目录已存在且非空，拒绝覆盖", final_dir))
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  edge_registry <- data.table(
    Edge_ID = c("E1", "E2", "E3"),
    Target = c("AG1_IG1", "AG2_IG1", "AG2_IG2"),
    Reference = c("AG2_IG1", "AG2_IG2", "AG3_IG2"),
    Target_State = c("C1", "C3", "C2"),
    Reference_State = c("C3", "C2", "C4"),
    Expected_Cancers = c(20L, 19L, 21L)
  )
  edge_dirs <- paths[c("e1_dir","e2_dir","e3_dir")]
  definition_path <- unname(paths["definition_registry"])
  paths <- data.table(Edge_ID=c("E1","E2","E3"),Source=file.path(edge_dirs,paste0("e",1:3,"_gsea_all_results.tsv.gz")),Reuse_Status=c("NEW_E1","REUSED_FROZEN_E2","NEW_E3"))
  if (any(!file.exists(paths$Source)) || any(file.info(paths$Source)$size <= 0)) {
    fail("E1/E2/E3逐癌GSEA汇总输入缺失或为空")
  }

  log_event("START", "three-edge GSEA combine begin", list(
    environment = "R4.5.1", edges = "E1,E2,E3",
    e2 = "reuse_frozen_no_rerun", output = final_dir))

  gsea <- rbindlist(lapply(seq_len(nrow(paths)), function(i) {
    z <- fread(paths$Source[i])
    if (!all(z$Edge_ID == paths$Edge_ID[i])) fail(paste(paths$Edge_ID[i], "来源边标识不一致"))
    drop_cols <- intersect(c("Target_State", "Reference_State"), names(z))
    if (length(drop_cols)) z[, (drop_cols) := NULL]
    merge(z, edge_registry[, .(Edge_ID, Target_State, Reference_State)],
          by = "Edge_ID", all.x = TRUE)
  }), use.names = TRUE, fill = TRUE)

  # 从三条边的逐癌结果统一重算跨癌票数，避免摘要格式差异。
  vote <- gsea[, .(
    Display_Name = Display_Name[1], Functional_Domain = Functional_Domain[1],
    Evidence_Role = Evidence_Role[1], N_Evaluable = uniqueN(Cancer),
    N_Target_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES > 0]),
    N_Reference_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES < 0]),
    Target_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES > 0])), collapse = ";"),
    Reference_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES < 0])), collapse = ";"),
    N_Target_High_NES = sum(NES > 0), N_Reference_High_NES = sum(NES < 0),
    Median_NES = median(NES), Median_Abs_NES = median(abs(NES)),
    Leading_Edge_Union = paste(sort(unique(unlist(strsplit(Leading_Edge, ";", fixed = TRUE)))), collapse = ";")
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection, Set_ID)]
  vote[, `:=`(
    Min_Recurrence = pmax(3L, ceiling(0.30 * N_Evaluable)),
    Min_Strong_Support = pmax(4L, ceiling(0.40 * N_Evaluable)),
    Max_Opposite_Allowed = floor(0.10 * N_Evaluable),
    Target_NES_Direction_Fraction = N_Target_High_NES / N_Evaluable,
    Reference_NES_Direction_Fraction = N_Reference_High_NES / N_Evaluable
  )]
  vote[, `:=`(
    Recurrent_Target_High = N_Target_High_Sig >= Min_Recurrence,
    Recurrent_Reference_High = N_Reference_High_Sig >= Min_Recurrence,
    Strong_Target_High = N_Target_High_Sig >= Min_Strong_Support &
      N_Reference_High_Sig <= Max_Opposite_Allowed,
    Strong_Reference_High = N_Reference_High_Sig >= Min_Strong_Support &
      N_Target_High_Sig <= Max_Opposite_Allowed
  )]
  vote[, Support_Class := fifelse(
    Strong_Target_High, "STRONG_TARGET_HIGH",
    fifelse(Strong_Reference_High, "STRONG_REFERENCE_HIGH",
      fifelse(Recurrent_Target_High, "RECURRENT_TARGET_HIGH",
        fifelse(Recurrent_Reference_High, "RECURRENT_REFERENCE_HIGH", "NOT_RECURRENT"))))]
  setorder(vote, Edge_ID, Sensitivity, Collection, -N_Target_High_Sig,
           -N_Reference_High_Sig, -Median_Abs_NES, Set_ID)

  candidates <- vote[
    Sensitivity == "NON_DEFINITION_RANK" &
      Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT") &
      Support_Class %chin% c("STRONG_TARGET_HIGH", "STRONG_REFERENCE_HIGH")
  ]
  candidates[, `:=`(
    Evidence_Tier = fifelse(Collection == "IMMUNE119_PRIMARY",
                            "PRIMARY_IMMUNE", "HALLMARK_CONTEXT"),
    Manual_Redundancy_Review = "PENDING",
    Fig3B_Eligibility = "PENDING_MANUAL_GATE"
  )]

  comparison <- dcast(vote[Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT")],
    Edge_ID + Target_State + Reference_State + Collection + Set_ID +
      Display_Name + Functional_Domain ~ Sensitivity,
    value.var = c("N_Target_High_Sig", "N_Reference_High_Sig", "Support_Class"))

  # 长表直接支持未来三比较×两方向复现矩阵；本轮不绘图。
  matrix_source <- rbind(
    vote[Sensitivity == "NON_DEFINITION_RANK" &
           Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT"), .(
      Edge_ID, Comparison = paste0(Target_State, " vs ", Reference_State),
      State = Target_State, Opposite_State = Reference_State, Collection,
      Set_ID, Display_Name, Functional_Domain, N_Evaluable,
      N_Directional_Significant = N_Target_High_Sig,
      N_Opposite_Significant = N_Reference_High_Sig,
      Direction_Fraction = Target_NES_Direction_Fraction,
      Strong = Strong_Target_High, Recurrent = Recurrent_Target_High,
      Side = "TARGET")],
    vote[Sensitivity == "NON_DEFINITION_RANK" &
           Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT"), .(
      Edge_ID, Comparison = paste0(Target_State, " vs ", Reference_State),
      State = Reference_State, Opposite_State = Target_State, Collection,
      Set_ID, Display_Name, Functional_Domain, N_Evaluable,
      N_Directional_Significant = N_Reference_High_Sig,
      N_Opposite_Significant = N_Target_High_Sig,
      Direction_Fraction = Reference_NES_Direction_Fraction,
      Strong = Strong_Reference_High, Recurrent = Recurrent_Reference_High,
      Side = "REFERENCE")]
  )
  setorder(matrix_source, Edge_ID, Collection, Set_ID, Side)

  summary <- vote[, .(
    N_Sets = .N,
    N_Recurrent_Target = sum(Recurrent_Target_High),
    N_Recurrent_Reference = sum(Recurrent_Reference_High),
    N_Strong_Target = sum(Strong_Target_High),
    N_Strong_Reference = sum(Strong_Reference_High)
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection)]

  definition <- norm_symbol(fread(definition_path)$gene_symbol)
  nondef_leading <- unique(unlist(strsplit(
    gsea[Sensitivity == "NON_DEFINITION_RANK", Leading_Edge], ";", fixed = TRUE)))
  cancer_counts <- unique(gsea[, .(N_Cancers = uniqueN(Cancer)), by = Edge_ID])
  cancer_counts <- merge(edge_registry[, .(Edge_ID, Expected_Cancers)], cancer_counts,
                         by = "Edge_ID", all.x = TRUE)
  qa <- data.table(
    Check = c("THREE_EDGES_EXACT", "EXPECTED_CANCER_COUNTS", "EXPECTED_COLLECTIONS",
              "TWO_SENSITIVITIES", "E2_REUSED_NOT_RERUN", "NO_C7",
              "NO_META_FIELDS", "FDR_RANGE_VALID",
              "NONDEFINITION_LEADING_EDGE_EXCLUDES_DEFINITION",
              "ALL_EDGE_QA_PASS", "OUTPUT_NONEMPTY"),
    Pass = c(
      setequal(unique(gsea$Edge_ID), c("E1", "E2", "E3")),
      all(cancer_counts$N_Cancers == cancer_counts$Expected_Cancers),
      setequal(unique(gsea$Collection), c("IMMUNE119_PRIMARY",
        "LM12_DEFINITION_CONTROL", "HALLMARK_CONTEXT")),
      setequal(unique(gsea$Sensitivity), c("FULL_RANK", "NON_DEFINITION_RANK")),
      paths[Edge_ID == "E2", Reuse_Status] == "REUSED_FROZEN_E2",
      !any(grepl("C7", gsea$Collection)),
      !any(c("Meta_FDR", "I2", "Tau2", "LOO") %in% names(gsea)),
      all(is.finite(gsea$FDR) & gsea$FDR >= 0 & gsea$FDR <= 1),
      length(intersect(nondef_leading, definition)) == 0L,
      all(vapply(c("E1", "E3"), function(e) {
        q <- file.path(edge_dirs[paste0(tolower(e),"_dir")], "machine_qa.tsv")
        file.exists(q) && all(fread(q)$Pass)
      }, logical(1))) && all(fread(file.path(edge_dirs["e2_dir"], "machine_qa.tsv"))$Pass),
      nrow(gsea) > 0L && nrow(vote) > 0L && nrow(matrix_source) > 0L
    )
  )
  if (!all(qa$Pass)) fail("三边GSEA统一汇总机器QA失败")

  paths[, `:=`(Size_Bytes = file.info(Source)$size,
               MD5 = unname(tools::md5sum(Source)))]
  fwrite(gsea, file.path(tmp_dir, "three_edge_gsea_all_results.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  fwrite(vote, file.path(tmp_dir, "three_edge_gsea_cross_cancer_votes.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(candidates, file.path(tmp_dir, "three_edge_nondefinition_strong_candidates.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(comparison, file.path(tmp_dir, "three_edge_definition_sensitivity_comparison.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(matrix_source, file.path(tmp_dir, "fig3b_recurrence_matrix_source.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(summary, file.path(tmp_dir, "three_edge_gsea_summary.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(paths, file.path(tmp_dir, "source_reuse_manifest.tsv"), sep = "\t", quote = TRUE)
  fwrite(qa, file.path(tmp_dir, "machine_qa.tsv"), sep = "\t", quote = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(tmp_dir, "sessionInfo.txt"))

  report <- c(
    "# Fig3三条局部边full-rank GSEA执行报告 V0.1", "", "日期：2026-08-06", "",
    "状态：`MACHINE_QA_PASS / MANUAL_GATE_PENDING / POST_FREEZE_EXPLORATORY / STAGE3_UNCHANGED`", "",
    "## 执行边界", "",
    "E1和E3新增全排序GSEA；E2复用既有冻结结果并记录MD5。所有单癌BH、双敏感性、gene-set和复现门槛统一；未重跑limma、未使用meta/LOO/跨癌BH。", "",
    "## 汇总", "",
    paste(apply(summary, 1L, function(x) sprintf(
      "%s（%s vs %s）%s × %s：target复现/强支持%s/%s，reference复现/强支持%s/%s。",
      x[["Edge_ID"]], x[["Target_State"]], x[["Reference_State"]],
      x[["Sensitivity"]], x[["Collection"]], x[["N_Recurrent_Target"]],
      x[["N_Strong_Target"]], x[["N_Recurrent_Reference"]],
      x[["N_Strong_Reference"]])), collapse = "\n\n"), "",
    "## 下一步", "",
    "只对NON_DEFINITION_RANK强支持候选进行语义去冗余和definition敏感性审核，冻结10–12行以内的Fig3B候选；本步骤不绘图。"
  )
  writeLines(report, file.path(tmp_dir, "execution_report.md"))

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) fail("无法提交三边GSEA统一汇总结果")
  log_event("DONE", "three-edge GSEA combine complete", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    rows = nrow(gsea), candidates = nrow(candidates),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)), output = final_dir))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)))
  quit(status = 1L)
})
}

# Purpose: Summarise reported RNA-protein concordance and discordance.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/13_run_mrna_rppa_concordance.R
stage_transcriptome_proteome_13_run_mrna_rppa_concordance <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LGG","LIHC","LUSC","OV","PRAD","SARC","SKCM","STAD","UCEC")
if(!dir.exists(paths["mrna_model_dir"]) || file.access(paths["mrna_model_dir"],5L)!=0L) stop("Model directory must be readable/searchable.")
required_files <- c(paths[setdiff(fields,c("mrna_model_dir","output_dir"))],file.path(paths["mrna_model_dir"],paste0(expected_cancers,"_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz")))

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

script_name <- "13_run_mrna_rppa_concordance.R"
start_time <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_012, environment())

collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
  if (length(x)) paste(x, collapse = ";") else ""
}

main <- function() {
  # 独立输出根仅在普通执行创建，源统计与原子提交连续保留。
  mrna_model_dir <- unname(paths["mrna_model_dir"])
  mrna_recurrence_file <- unname(paths["mrna_recurrence"])
  mrna_display_file <- unname(paths["mrna_display"])
  rppa_within_file <- unname(paths["rppa_within"])
  rppa_recurrence_file <- unname(paths["rppa_recurrence"])
  rppa_display_file <- unname(paths["rppa_display"])
  cancer_gate_file <- unname(paths["cancer_gate"])
  dictionary_file <- unname(paths["antibody_dictionary"])
  output_dir <- file.path(out_dir, "01_mrna_rppa_concordance")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(output_dir)) stop("Output subdirectory already exists.")
  tmp_dir <- paste0(output_dir, ".tmp_", Sys.getpid())
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  expected_cancers <- c(
    "BLCA", "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LGG",
    "LIHC", "LUSC", "OV", "PRAD", "SARC", "SKCM", "STAD", "UCEC"
  )
  mrna_model_files <- file.path(
    mrna_model_dir,
    paste0(
      expected_cancers,
      "_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz"
    )
  )
  required <- c(
    mrna_model_files, mrna_recurrence_file, mrna_display_file,
    rppa_within_file, rppa_recurrence_file, rppa_display_file,
    cancer_gate_file, dictionary_file
  )
  if (any(!file.exists(required)) ||
      any(file.info(required)$size <= 0)) {
    stop("存在缺失或空的跨层审计冻结输入")
  }

  log_event("START", "开始mRNA–RPPA C1 vs C4一致性/矛盾审计", list(
    environment = "R4.5.1",
    cancers = length(expected_cancers),
    design = "C1_vs_C4",
    mrna_threshold = "FDR<0.05_AND_ABS_ESTIMATE>1.0",
    rppa_threshold = "FDR<0.05_NO_ABS_EFFECT_CUTOFF",
    cross_cancer_bh = "NOT_PERFORMED",
    meta_loo_i2 = "NOT_PERFORMED",
    output = output_dir
  ))

  cancer_gate <- fread(cancer_gate_file)
  formal_cancers <- cancer_gate[Formal_Cancer == TRUE, Cancer]
  if (!setequal(formal_cancers, expected_cancers) ||
      length(formal_cancers) != 15L) {
    stop("RPPA正式15癌种集合与合同不一致")
  }

  dictionary <- fread(dictionary_file)
  if (nrow(dictionary) != 258L ||
      anyDuplicated(dictionary$standardized_feature)) {
    stop("RPPA冻结字典不是258个唯一feature")
  }
  dictionary[, Gene := toupper(trimws(as.character(gene_symbol)))]
  dictionary[, Base_Direct_Eligible :=
               feature_class == "TOTAL" &
               direct_mrna_pairable == TRUE &
               target_count == 1L &
               !is.na(Gene) & nzchar(Gene)]
  dictionary[
    Base_Direct_Eligible == TRUE,
    N_Direct_Total_Features_Per_Gene := .N,
    by = Gene
  ]
  dictionary[
    is.na(N_Direct_Total_Features_Per_Gene),
    N_Direct_Total_Features_Per_Gene := 0L
  ]
  dictionary[, Cross_Layer_Route := fcase(
    feature_class == "PHOSPHO",
    "EXECUTION_LAYER_CONTEXT_PHOSPHO",
    feature_class == "CLEAVED",
    "EXECUTION_LAYER_CONTEXT_CLEAVED",
    Base_Direct_Eligible == TRUE &
      N_Direct_Total_Features_Per_Gene == 1L,
    "DIRECT_TOTAL_UNIQUE_GENE",
    Base_Direct_Eligible == TRUE &
      N_Direct_Total_Features_Per_Gene > 1L,
    "UNEVALUABLE_DUPLICATE_TOTAL_GENE",
    feature_class == "TOTAL",
    "UNEVALUABLE_TOTAL_MAPPING",
    default = "UNEVALUABLE_OTHER"
  )]

  rppa_recurrence <- fread(rppa_recurrence_file)
  rppa_recurrence[, Gene := toupper(trimws(as.character(Gene_Symbol)))]
  rppa_display <- fread(rppa_display_file)
  rppa_display[, Gene := toupper(trimws(as.character(Gene_Symbol)))]
  rppa_display_small <- rppa_display[, .(
    Feature_ID,
    RPPA_Formal_Display = TRUE,
    RPPA_Display_Order = Display_Order,
    RPPA_TCPA_Pathway_Membership = TCPA_Pathway_Membership,
    RPPA_TCPA_Mapping_Status = TCPA_Mapping_Status
  )]
  feature_routing <- merge(
    dictionary[, .(
      Feature_ID = standardized_feature,
      Feature_Label = raw_feature,
      Feature_Class = feature_class,
      Gene,
      Gene_Symbol = gene_symbol,
      Site = site,
      Target_Count = target_count,
      Direct_mRNA_Pairable = direct_mrna_pairable,
      Dictionary_Mapping_Status = mapping_status,
      Validation_Status = validation_status,
      Dictionary_Cross_Layer_Restriction = cross_layer_restriction,
      N_Direct_Total_Features_Per_Gene,
      Cross_Layer_Route
    )],
    rppa_recurrence[, .(
      Feature_ID, Definition_Overlap,
      RPPA_N_Cancer_Estimable = N_Cancer_Estimable,
      RPPA_N_C1_High_Significant = N_C1_High_Significant,
      RPPA_N_C4_High_Significant = N_C4_High_Significant,
      RPPA_Recurrent_Class = Recurrent_Class,
      RPPA_Formal_Display = Formal_Display
    )],
    by = "Feature_ID", all.x = TRUE
  )
  feature_routing <- merge(
    feature_routing, rppa_display_small,
    by = "Feature_ID", all.x = TRUE,
    suffixes = c("", "_Display_Table")
  )
  feature_routing[
    is.na(RPPA_Formal_Display),
    RPPA_Formal_Display := FALSE
  ]
  feature_routing[
    is.na(RPPA_Formal_Display_Display_Table),
    RPPA_Formal_Display_Display_Table := FALSE
  ]
  setorder(feature_routing, Cross_Layer_Route, Gene, Feature_ID)

  pair_registry <- feature_routing[
    Cross_Layer_Route == "DIRECT_TOTAL_UNIQUE_GENE",
    .(
      Gene, Feature_ID,
      RPPA_Feature_ID = Feature_ID,
      Feature_Label, Feature_Class,
      Gene_Symbol, Definition_Overlap,
      Dictionary_Mapping_Status, Validation_Status,
      RPPA_N_Cancer_Estimable,
      RPPA_N_C1_High_Significant,
      RPPA_N_C4_High_Significant,
      RPPA_Recurrent_Class,
      RPPA_Formal_Display,
      RPPA_Display_Order,
      RPPA_TCPA_Pathway_Membership,
      RPPA_TCPA_Mapping_Status
    )
  ]
  if (anyDuplicated(pair_registry$Gene) ||
      anyDuplicated(pair_registry$Feature_ID)) {
    stop("直接total protein–mRNA配对不是gene与feature双向唯一")
  }

  mrna_keep <- c(
    "Cancer", "Feature_ID", "Feature_Label", "Feature_Class",
    "Definition_Overlap", "Design_ID", "Analysis_Set", "Model_Type",
    "Target", "Reference", "Estimate", "SE", "P_Value", "FDR",
    "Status", "Failure_Reason"
  )
  mrna_rows <- rbindlist(lapply(
    seq_along(mrna_model_files),
    function(i) {
      x <- fread(mrna_model_files[i], select = mrna_keep)
      x <- x[
        Design_ID == "C1_vs_C4" &
          Analysis_Set == "REPRESENTATIVE_MAIN" &
          Model_Type == "UNADJUSTED"
      ]
      x[, Gene := toupper(trimws(as.character(Feature_ID)))]
      x <- x[Gene %chin% pair_registry$Gene]
      if (anyDuplicated(x$Gene)) {
        stop("单癌种mRNA直接配对gene重复: ", expected_cancers[i])
      }
      x
    }
  ), use.names = TRUE, fill = TRUE)

  rppa_within <- fread(rppa_within_file)
  rppa_within <- rppa_within[
    Cancer %chin% expected_cancers &
      Feature_ID %chin% pair_registry$Feature_ID
  ]
  rppa_within <- merge(
    rppa_within,
    pair_registry[, .(Gene, Feature_ID)],
    by = "Feature_ID", all.x = TRUE
  )
  if (anyDuplicated(rppa_within[, .(Cancer, Gene)])) {
    stop("RPPA直接配对Cancer × Gene重复")
  }

  mrna_pair <- mrna_rows[, .(
    Cancer, Gene,
    mRNA_Feature_ID = Feature_ID,
    mRNA_Definition_Overlap = Definition_Overlap,
    mRNA_Estimate = Estimate,
    mRNA_SE = SE,
    mRNA_P_Value = P_Value,
    mRNA_FDR = FDR,
    mRNA_Status = Status,
    mRNA_Failure_Reason = Failure_Reason
  )]
  rppa_pair <- rppa_within[, .(
    Cancer, Gene,
    RPPA_Feature_ID = Feature_ID,
    RPPA_Feature_Label = Feature_Label,
    RPPA_Definition_Overlap = Definition_Overlap,
    RPPA_Estimate = Estimate,
    RPPA_SE = SE,
    RPPA_P_Value = P_Value,
    RPPA_FDR = FDR,
    RPPA_Status = Status,
    RPPA_Failure_Reason = Failure_Reason
  )]
  cancer_gene <- CJ(
    Cancer = expected_cancers,
    Gene = pair_registry$Gene,
    unique = TRUE
  )
  cancer_gene <- merge(
    cancer_gene, mrna_pair,
    by = c("Cancer", "Gene"), all.x = TRUE
  )
  cancer_gene <- merge(
    cancer_gene, rppa_pair,
    by = c("Cancer", "Gene"), all.x = TRUE
  )
  cancer_gene <- merge(
    cancer_gene,
    pair_registry[, .(
      Gene, RPPA_Feature_ID_Registry = Feature_ID,
      RPPA_Gene_Definition_Overlap = Definition_Overlap
    )],
    by = "Gene", all.x = TRUE
  )

  cancer_gene[, `:=`(
    mRNA_Evaluable =
      mRNA_Status == "ESTIMABLE" &
      is.finite(mRNA_Estimate) &
      is.finite(mRNA_FDR),
    RPPA_Evaluable =
      RPPA_Status == "ESTIMABLE" &
      is.finite(RPPA_Estimate) &
      is.finite(RPPA_FDR)
  )]
  cancer_gene[, mRNA_Sig_Direction := fcase(
    mRNA_Evaluable & mRNA_FDR < 0.05 & mRNA_Estimate > 1.0,
    "C1_HIGH",
    mRNA_Evaluable & mRNA_FDR < 0.05 & mRNA_Estimate < -1.0,
    "C4_HIGH",
    default = "NOT_SIGNIFICANT"
  )]
  cancer_gene[, RPPA_Sig_Direction := fcase(
    RPPA_Evaluable & RPPA_FDR < 0.05 & RPPA_Estimate > 0,
    "C1_HIGH",
    RPPA_Evaluable & RPPA_FDR < 0.05 & RPPA_Estimate < 0,
    "C4_HIGH",
    default = "NOT_SIGNIFICANT"
  )]
  cancer_gene[, Cancer_Level_Audit_Class := fcase(
    !mRNA_Evaluable | !RPPA_Evaluable,
    "UNEVALUABLE_COVERAGE",
    mRNA_Sig_Direction != "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction != "NOT_SIGNIFICANT" &
      mRNA_Sig_Direction == RPPA_Sig_Direction,
    "CONCORDANT_SIGNIFICANT",
    mRNA_Sig_Direction != "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction != "NOT_SIGNIFICANT" &
      mRNA_Sig_Direction != RPPA_Sig_Direction,
    "DISCORDANT_SIGNIFICANT",
    mRNA_Sig_Direction != "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction == "NOT_SIGNIFICANT",
    "RNA_ONLY_SIGNIFICANT",
    mRNA_Sig_Direction == "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction != "NOT_SIGNIFICANT",
    "PROTEIN_ONLY_SIGNIFICANT",
    default = "BOTH_NONSIGNIFICANT"
  )]
  cancer_gene[, Effect_Sign_Relation := fcase(
    !mRNA_Evaluable | !RPPA_Evaluable,
    "UNEVALUABLE",
    sign(mRNA_Estimate) == sign(RPPA_Estimate),
    "SAME_DIRECTION",
    default = "OPPOSITE_DIRECTION"
  )]
  cancer_gene[, Cancer_Order := match(Cancer, expected_cancers)]
  setorder(cancer_gene, Gene, Cancer_Order)

  mrna_recurrence <- fread(mrna_recurrence_file)
  mrna_recurrence <- mrna_recurrence[
    Design_ID == "C1_vs_C4"
  ]
  mrna_recurrence[, Gene := toupper(trimws(as.character(Gene)))]
  mrna_recurrence_small <- mrna_recurrence[, .(
    Gene,
    mRNA_Recurrence_Class = Recurrence_Class,
    mRNA_N_Evaluable_Cancers = N_Evaluable_Cancers,
    mRNA_N_C1_High_Significant = N_Sig_Up,
    mRNA_N_C4_High_Significant = N_Sig_Down,
    mRNA_Recurrent_C1_High = Recurrent_Up,
    mRNA_Recurrent_C4_High = Recurrent_Down,
    mRNA_Bidirectional_Recurrent = Bidirectional_Recurrent,
    mRNA_Definition_Overlap_Recurrence = Definition_Overlap
  )]
  rppa_recurrence_small <- rppa_recurrence[, .(
    Gene,
    RPPA_Feature_ID = Feature_ID,
    RPPA_Recurrence_Class = Recurrent_Class,
    RPPA_N_Evaluable_Cancers = N_Cancer_Estimable,
    RPPA_Definition_Overlap_Recurrence = Definition_Overlap
  )]
  gene_summary <- merge(
    pair_registry,
    mrna_recurrence_small,
    by = "Gene", all.x = TRUE
  )
  gene_summary <- merge(
    gene_summary,
    rppa_recurrence_small,
    by = c("Gene", "RPPA_Feature_ID"),
    all.x = TRUE
  )

  per_gene_counts <- cancer_gene[, .(
    N_Formal_Cancers = .N,
    N_Both_Evaluable = sum(mRNA_Evaluable & RPPA_Evaluable),
    N_Concordant_Significant =
      sum(Cancer_Level_Audit_Class == "CONCORDANT_SIGNIFICANT"),
    N_Discordant_Significant =
      sum(Cancer_Level_Audit_Class == "DISCORDANT_SIGNIFICANT"),
    N_RNA_Only_Significant =
      sum(Cancer_Level_Audit_Class == "RNA_ONLY_SIGNIFICANT"),
    N_Protein_Only_Significant =
      sum(Cancer_Level_Audit_Class == "PROTEIN_ONLY_SIGNIFICANT"),
    N_Both_Nonsignificant =
      sum(Cancer_Level_Audit_Class == "BOTH_NONSIGNIFICANT"),
    N_Unevaluable =
      sum(Cancer_Level_Audit_Class == "UNEVALUABLE_COVERAGE"),
    Concordant_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "CONCORDANT_SIGNIFICANT"]
    ),
    Discordant_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "DISCORDANT_SIGNIFICANT"]
    ),
    RNA_Only_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "RNA_ONLY_SIGNIFICANT"]
    ),
    Protein_Only_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "PROTEIN_ONLY_SIGNIFICANT"]
    )
  ), by = Gene]
  gene_summary <- merge(
    gene_summary, per_gene_counts, by = "Gene", all.x = TRUE
  )

  gene_summary[, mRNA_Formal_Direction := fcase(
    mRNA_Recurrence_Class == "RECURRENT_UP", "C1_HIGH",
    mRNA_Recurrence_Class == "RECURRENT_DOWN", "C4_HIGH",
    default = "NO_UNIQUE_FORMAL_DIRECTION"
  )]
  gene_summary[, RPPA_Formal_Direction := fcase(
    RPPA_Recurrence_Class == "C1_HIGH", "C1_HIGH",
    RPPA_Recurrence_Class == "C4_HIGH", "C4_HIGH",
    default = "NO_UNIQUE_FORMAL_DIRECTION"
  )]
  gene_summary[, `:=`(
    mRNA_Formal_Supported =
      mRNA_Formal_Direction %chin% c("C1_HIGH", "C4_HIGH"),
    RPPA_Formal_Supported =
      RPPA_Formal_Direction %chin% c("C1_HIGH", "C4_HIGH")
  )]
  gene_summary[, Cross_Layer_Support_Class := fcase(
    mRNA_Recurrence_Class == "BIDIRECTIONAL_RECURRENT" |
      RPPA_Recurrence_Class == "BIDIRECTIONAL_RECURRENT",
    "BIDIRECTIONAL_OR_UNRESOLVED",
    mRNA_Formal_Supported & RPPA_Formal_Supported &
      mRNA_Formal_Direction == RPPA_Formal_Direction,
    "CONCORDANT_SUPPORTED",
    mRNA_Formal_Supported & RPPA_Formal_Supported &
      mRNA_Formal_Direction != RPPA_Formal_Direction,
    "DISCORDANT_SUPPORTED",
    mRNA_Formal_Supported & !RPPA_Formal_Supported,
    "RNA_ONLY_SUPPORTED",
    !mRNA_Formal_Supported & RPPA_Formal_Supported,
    "PROTEIN_ONLY_SUPPORTED",
    default = "NO_CROSS_CANCER_SUPPORT"
  )]

  mrna_display <- fread(mrna_display_file)
  mrna_display[, Gene := toupper(trimws(as.character(Gene)))]
  mrna_display_small <- mrna_display[, .(
    mRNA_Selected_For_Cytoscape = any(Selected_For_Cytoscape),
    mRNA_Display_State = collapse_values(State),
    mRNA_Display_Direction = collapse_values(Direction),
    mRNA_Primary_Cluster_ID = collapse_values(Primary_Cluster_ID),
    mRNA_Display_Assignment = collapse_values(Display_Assignment)
  ), by = Gene]
  gene_summary <- merge(
    gene_summary, mrna_display_small,
    by = "Gene", all.x = TRUE
  )
  gene_summary[
    is.na(mRNA_Selected_For_Cytoscape),
    mRNA_Selected_For_Cytoscape := FALSE
  ]
  gene_summary[, Definition_Overlap_Any :=
                 Definition_Overlap %in% TRUE |
                 mRNA_Definition_Overlap_Recurrence %in% TRUE |
                 RPPA_Definition_Overlap_Recurrence %in% TRUE]
  gene_summary[, Interpretation_Boundary := fifelse(
    Definition_Overlap_Any,
    "DEFINITION_OVERLAP_AXIS_CONSISTENCY_ONLY",
    "ORTHOGONAL_COMPATIBILITY_NOT_CAUSAL"
  )]
  setorder(
    gene_summary, Cross_Layer_Support_Class,
    -N_Concordant_Significant, -N_Discordant_Significant, Gene
  )

  execution_context <- feature_routing[
    Cross_Layer_Route %chin% c(
      "EXECUTION_LAYER_CONTEXT_PHOSPHO",
      "EXECUTION_LAYER_CONTEXT_CLEAVED"
    )
  ]
  execution_context[, Interpretation_Boundary :=
                      "EXECUTION_LAYER_CONTEXT_NO_REQUIRED_MRNA_DIRECTION"]

  cancer_class_summary <- cancer_gene[, .(
    N_Cancer_Gene_Cells = .N,
    N_Genes = uniqueN(Gene),
    N_Cancers = uniqueN(Cancer)
  ), by = Cancer_Level_Audit_Class]
  support_class_summary <- gene_summary[, .(
    N_Genes = .N,
    N_mRNA_Display_Genes = sum(mRNA_Selected_For_Cytoscape),
    N_RPPA_Display_Features = sum(RPPA_Formal_Display %in% TRUE),
    N_Definition_Overlap = sum(Definition_Overlap_Any)
  ), by = Cross_Layer_Support_Class]
  display_overlap <- gene_summary[
    mRNA_Selected_For_Cytoscape == TRUE |
      RPPA_Formal_Display == TRUE,
    .(
      Gene, RPPA_Feature_ID,
      Cross_Layer_Support_Class,
      mRNA_Formal_Direction, RPPA_Formal_Direction,
      N_Concordant_Significant, N_Discordant_Significant,
      mRNA_Selected_For_Cytoscape,
      mRNA_Display_State, mRNA_Primary_Cluster_ID,
      RPPA_Formal_Display,
      RPPA_TCPA_Pathway_Membership,
      Definition_Overlap_Any,
      Interpretation_Boundary
    )
  ]

  threshold_registry <- data.table(
    Parameter = c(
      "DESIGN_ID", "FORMAL_CANCERS",
      "DIRECT_MAPPING_FEATURE_CLASS",
      "DIRECT_MAPPING_TARGET_COUNT",
      "DIRECT_MAPPING_GENE_MULTIPLICITY",
      "MRNA_WITHIN_CANCER_FDR",
      "MRNA_ABS_EFFECT",
      "MRNA_RECURRENCE",
      "RPPA_WITHIN_CANCER_FDR",
      "RPPA_ABS_EFFECT",
      "RPPA_RECURRENCE",
      "PHOSPHO_CLEAVED_POLICY",
      "CROSS_CANCER_BH", "META_LOO_I2", "PURITY"
    ),
    Value = c(
      "C1_vs_C4", "15",
      "TOTAL_ONLY", "1", "EXACTLY_ONE_TOTAL_FEATURE_PER_GENE",
      "0.05", ">1.0", "AT_LEAST_30_PERCENT_EVALUABLE_CANCERS",
      "0.05", "NONE", "AT_LEAST_3_OF_15_CANCERS",
      "EXECUTION_LAYER_CONTEXT_ONLY",
      "NOT_PERFORMED", "NOT_PERFORMED", "NOT_USED"
    )
  )

  files <- c(
    routing = file.path(
      tmp_dir,
      "stage3_mrna_rppa_feature_routing_v0.1_20260730.tsv"
    ),
    pairs = file.path(
      tmp_dir,
      "stage3_mrna_rppa_unique_total_pair_registry_v0.1_20260730.tsv"
    ),
    cells = file.path(
      tmp_dir,
      "stage3_mrna_rppa_cancer_gene_audit_v0.1_20260730.tsv.gz"
    ),
    genes = file.path(
      tmp_dir,
      "stage3_mrna_rppa_gene_summary_v0.1_20260730.tsv"
    ),
    execution = file.path(
      tmp_dir,
      "stage3_mrna_rppa_execution_layer_context_v0.1_20260730.tsv"
    ),
    cancer_summary = file.path(
      tmp_dir,
      "stage3_mrna_rppa_cancer_class_summary_v0.1_20260730.tsv"
    ),
    support_summary = file.path(
      tmp_dir,
      "stage3_mrna_rppa_support_class_summary_v0.1_20260730.tsv"
    ),
    display = file.path(
      tmp_dir,
      "stage3_mrna_rppa_display_overlap_audit_v0.1_20260730.tsv"
    ),
    thresholds = file.path(
      tmp_dir,
      "stage3_mrna_rppa_threshold_registry_v0.1_20260730.tsv"
    )
  )
  fwrite(feature_routing, files[["routing"]], sep = "\t", na = "NA")
  fwrite(pair_registry, files[["pairs"]], sep = "\t", na = "NA")
  fwrite(cancer_gene, files[["cells"]], sep = "\t", na = "NA")
  fwrite(gene_summary, files[["genes"]], sep = "\t", na = "NA")
  fwrite(execution_context, files[["execution"]], sep = "\t", na = "NA")
  fwrite(
    cancer_class_summary,
    files[["cancer_summary"]], sep = "\t", na = "NA"
  )
  fwrite(
    support_class_summary,
    files[["support_summary"]], sep = "\t", na = "NA"
  )
  fwrite(display_overlap, files[["display"]], sep = "\t", na = "NA")
  fwrite(
    threshold_registry,
    files[["thresholds"]], sep = "\t", na = "NA"
  )

  expected_cell_classes <- c(
    "CONCORDANT_SIGNIFICANT", "DISCORDANT_SIGNIFICANT",
    "RNA_ONLY_SIGNIFICANT", "PROTEIN_ONLY_SIGNIFICANT",
    "BOTH_NONSIGNIFICANT", "UNEVALUABLE_COVERAGE"
  )
  expected_support_classes <- c(
    "CONCORDANT_SUPPORTED", "DISCORDANT_SUPPORTED",
    "RNA_ONLY_SUPPORTED", "PROTEIN_ONLY_SUPPORTED",
    "NO_CROSS_CANCER_SUPPORT", "BIDIRECTIONAL_OR_UNRESOLVED"
  )
  qa <- data.table(
    Check = c(
      "FORMAL_CANCERS_15",
      "RPPA_FEATURES_258_ROUTED",
      "DIRECT_PAIR_GENES_187",
      "DIRECT_PAIR_GENE_FEATURE_UNIQUE",
      "NO_PHOSPHO_CLEAVED_IN_DIRECT_PAIRS",
      "CANCER_GENE_CELLS_15_BY_187",
      "CANCER_LEVEL_CLASSES_EXHAUSTIVE",
      "GENE_SUPPORT_CLASSES_EXHAUSTIVE",
      "ALL_DIRECT_GENES_HAVE_SUMMARY",
      "EXECUTION_LAYER_FEATURES_60",
      "NO_CROSS_CANCER_BH_META_FIELDS",
      "ALL_OUTPUTS_NONEMPTY"
    ),
    Observed = c(
      length(formal_cancers),
      nrow(feature_routing),
      nrow(pair_registry),
      !anyDuplicated(pair_registry$Gene) &&
        !anyDuplicated(pair_registry$Feature_ID),
      all(pair_registry$Feature_Class == "TOTAL"),
      nrow(cancer_gene),
      all(cancer_gene$Cancer_Level_Audit_Class %chin%
            expected_cell_classes),
      all(gene_summary$Cross_Layer_Support_Class %chin%
            expected_support_classes),
      nrow(gene_summary),
      nrow(execution_context),
      !any(grepl(
        "META|LOO|I2|CROSS_CANCER_FDR",
        names(gene_summary), ignore.case = TRUE
      )),
      all(file.exists(files)) &&
        all(file.info(files)$size > 0)
    ),
    Expected = c(
      15L, 258L, 187L, TRUE, TRUE,
      15L * 187L, TRUE, TRUE, 187L, 60L, TRUE, TRUE
    )
  )
  qa[, Pass := as.character(Observed) == as.character(Expected)]
  qa_file <- file.path(
    tmp_dir,
    "stage3_mrna_rppa_machine_qa_v0.1_20260730.tsv"
  )
  fwrite(qa, qa_file, sep = "\t", na = "NA")
  if (!all(qa$Pass)) {
    print(qa)
    stop("mRNA–RPPA一致性/矛盾审计machine QA失败")
  }

  input_manifest <- data.table(
    Role = c(
      rep("MRNA_CANCER_MODEL", length(mrna_model_files)),
      "MRNA_RECURRENCE", "MRNA_DISPLAY",
      "RPPA_WITHIN_CANCER", "RPPA_RECURRENCE", "RPPA_DISPLAY",
      "RPPA_CANCER_GATE", "RPPA_DICTIONARY"
    ),
    Path = required,
    Bytes = file.info(required)$size,
    MD5 = unname(tools::md5sum(required))
  )
  fwrite(
    input_manifest,
    file.path(
      tmp_dir,
      "stage3_mrna_rppa_input_manifest_v0.1_20260730.tsv"
    ),
    sep = "\t", na = "NA"
  )

  dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, output_dir)) {
    stop("无法原子冻结mRNA–RPPA审计目录")
  }
  log_event("DONE", "mRNA–RPPA C1 vs C4一致性/矛盾审计完成", list(
    runtime_sec = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    )), 1),
    direct_pair_genes = nrow(pair_registry),
    cancer_gene_cells = nrow(cancer_gene),
    concordant_supported = gene_summary[
      Cross_Layer_Support_Class == "CONCORDANT_SUPPORTED", .N
    ],
    discordant_supported = gene_summary[
      Cross_Layer_Support_Class == "DISCORDANT_SUPPORTED", .N
    ],
    execution_layer_features = nrow(execution_context),
    output = output_dir
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
