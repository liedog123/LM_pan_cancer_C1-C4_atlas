#!/usr/bin/env Rscript
# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Task: lncrna_rppa; step: s05
# Inputs (DCF fields): model_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lncrna_rppa --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 05_freeze_lncrna_recurrence.R --config PATH [--check-config]\n       Rscript 05_freeze_lncrna_recurrence.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("model_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}


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
