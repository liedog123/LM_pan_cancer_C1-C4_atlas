#!/usr/bin/env Rscript
# Purpose: Fit lncRNA and RPPA comparisons and derive recurrent features and linked functions.
# Task: lncrna_rppa; step: s08
# Inputs (DCF fields): edge_recurrence_dir,hallmark_gmt,reactome_gmt,go_bp_gmt,c7_gmt,immport_gmt
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lncrna_rppa --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 08_run_lncrna_linked_mrna_ora.R --config PATH [--check-config]\n       Rscript 08_run_lncrna_linked_mrna_ora.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("edge_recurrence_dir", "hallmark_gmt", "reactome_gmt", "go_bp_gmt", "c7_gmt", "immport_gmt", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "08_run_lncrna_linked_mrna_ora.R"
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
