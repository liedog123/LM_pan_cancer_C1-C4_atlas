#!/usr/bin/env Rscript
# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Task: lncrna_rppa; step: s06
# Inputs (DCF fields): patient_manifest,lncrna_registry,recurrent_candidates,mrna_rdata,lncrna_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lncrna_rppa --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 06_run_lncrna_gene_coexpression.R --config PATH [--check-config]\n       Rscript 06_run_lncrna_gene_coexpression.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_manifest", "lncrna_registry", "recurrent_candidates", "mrna_rdata", "lncrna_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}


suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "06_run_lncrna_gene_coexpression.R"
start_time <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extra <- if (length(fields)) {
    paste(vapply(names(fields), function(nm) {
      sprintf(
        '%s="%s"', nm,
        gsub("[\r\n\t]", " ", as.character(fields[[nm]]))
      )
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
