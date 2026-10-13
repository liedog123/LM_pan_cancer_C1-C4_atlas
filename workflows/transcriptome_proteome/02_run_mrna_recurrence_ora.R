#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 02_run_mrna_recurrence_ora.R --config PATH [--check-config]\n       Rscript 02_run_mrna_recurrence_ora.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("model_dir", "hallmark_gmt", "reactome_gmt", "go_bp_gmt", "c7_gmt", "immport_gmt", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
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
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

suppressPackageStartupMessages(library(data.table))

options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "02_run_mrna_recurrence_ora.R"
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
