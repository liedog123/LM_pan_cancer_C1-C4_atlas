#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 07_build_lncrna_edge_recurrence.R --config PATH [--check-config]\n       Rscript 07_build_lncrna_edge_recurrence.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("coexpression_dir", "recurrent_candidates", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}


suppressPackageStartupMessages({
  library(data.table)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "07_build_lncrna_edge_recurrence.R"
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
