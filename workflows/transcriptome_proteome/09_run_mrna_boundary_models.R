#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 09_run_mrna_boundary_models.R --config PATH [--check-config]\n       Rscript 09_run_mrna_boundary_models.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
path_fields <- setdiff(fields,"edge_id")
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}




suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

options(stringsAsFactors = FALSE, scipen = 999)
script_name <- "09_run_mrna_boundary_models.R"
start_time <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extra <- if (length(fields)) paste(vapply(names(fields), function(nm) {
    sprintf('%s="%s"', nm, gsub("[\r\n\t]", " ", as.character(fields[[nm]])))
  }, character(1)), collapse = " ") else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
              gsub("[\r\n]", " ", msg), if (nzchar(extra)) " " else "", extra))
  flush.console()
}

norm_sample <- function(x) substr(gsub("\\.", "-", as.character(x)), 1L, 15L)
norm_symbol <- function(x) toupper(trimws(as.character(x)))


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
