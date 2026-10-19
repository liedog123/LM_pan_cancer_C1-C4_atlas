#!/usr/bin/env Rscript
# Purpose: Prepare and test state-adjusted programme-feature associations within cancers.
# Task: programme_feature_associations; step: s05
# Inputs (DCF fields): patient,pilot,mrna,definition,program_genes,program_clusters,rppa,rppa_dict,tcpa,tmb,mutation,cnv,estimate,cibersort,mcpcounter,drug
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task programme_feature_associations --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 05_build_shared_numeric_inputs.R --config PATH [--check-config]\n       Rscript 05_build_shared_numeric_inputs.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { usage(); quit(status = 0L) }
if (!(length(args) %in% c(2L,3L)) || args[1L] != "--config" ||
    !nzchar(args[2L]) || (length(args)==3L && args[3L] != "--check-config")) {
  usage(); stop("Expected --config PATH [--check-config].")
}
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link) && !file.exists(path) && !dir.exists(path)) stop("Dangling symbolic link: ", path)
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
config_lines <- readLines(config_path, warn = FALSE)
config_keys <- trimws(sub(":.*$", "", config_lines[
  grepl("^[^[:space:]#][^:]*:", config_lines)]))
if(anyDuplicated(config_keys)) stop("DCF fields must not repeat.")
cfg <- read.dcf(config_path)
fields <- c("patient", "pilot", "mrna", "definition", "program_genes", "program_clusters", "rppa", "rppa_dict", "tcpa", "tmb", "mutation", "cnv", "estimate", "cibersort", "mcpcounter", "drug", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient", "pilot", "mrna", "definition", "program_genes", "program_clusters", "rppa", "rppa_dict", "tcpa", "tmb", "mutation", "cnv", "estimate", "cibersort", "mcpcounter", "drug")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),config_path,package_root)
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

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

script_name <- "03_build_stage1_shared_numeric_inputs_v0.1_20260810.R"
started <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extras <- if (length(fields)) paste(vapply(names(fields), function(k) {
    sprintf('%s="%s"', k, gsub('"', "'", as.character(fields[[k]]), fixed = TRUE))
  }, character(1)), collapse = " ") else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
              gsub('"', "'", msg, fixed = TRUE),
              if (nzchar(extras)) " " else "", extras))
  flush.console()
}

norm_id <- function(x) substr(toupper(trimws(as.character(x))), 1L, 15L)
norm_gene <- function(x) toupper(trimws(as.character(x)))

z_rows <- function(m) {
  m <- as.matrix(m)
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- rowSds(m, na.rm = TRUE)
  z <- sweep(sweep(m, 1L, mu, "-"), 1L, sdv, "/")
  z[!is.finite(z)] <- NA_real_
  z
}

score_sets <- function(feature_by_sample, sets, min_features = 2L) {
  z <- z_rows(feature_by_sample)
  ans <- lapply(names(sets), function(nm) {
    spec <- sets[[nm]]
    ids <- intersect(names(spec), rownames(z))
    if (length(ids) < min_features) return(rep(NA_real_, ncol(z)))
    signed <- sweep(z[ids, , drop = FALSE], 1L, as.numeric(spec[ids]), "*")
    colMeans(signed, na.rm = TRUE)
  })
  out <- do.call(cbind, ans)
  colnames(out) <- names(sets)
  rownames(out) <- colnames(z)
  out[!is.finite(out)] <- NA_real_
  out
}

left_matrix <- function(ids, dt, id_col, value_cols, prefix = "") {
  idx <- match(ids, dt[[id_col]])
  out <- as.matrix(dt[idx, ..value_cols])
  storage.mode(out) <- "numeric"
  rownames(out) <- ids
  colnames(out) <- paste0(prefix, value_cols)
  out
}

main <- function() {
  # Preserve every shared-input component and source scientific gate.
  out_dir <- unname(paths["output_dir"])
  bundle_dir <- file.path(out_dir, "bundles")
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  p <- as.list(input_paths)
  bad <- names(p)[!file.exists(unlist(p)) | file.info(unlist(p))$size <= 0]
  if (length(bad)) stop("输入缺失或为空: ", paste(bad, collapse = ","))

  log_event("START", "构建Stage 1共用数值输入", list(
    project = "LM_pan_cancer", output = out_dir,
    formal_model = "NONE", rna_topvar = 200L
  ))

  patient <- fread(p$patient)
  pilot <- fread(p$pilot)
  pilot_cancers <- pilot[selected_for_mofa_pilot == TRUE, Cancer]
  if (!setequal(pilot_cancers, c("BRCA", "LGG", "COAD", "PRAD", "OV"))) stop("试点癌种漂移")

  definition <- fread(p$definition)
  definition_genes <- unique(norm_gene(definition$gene_symbol))
  if (length(definition_genes) != 497L) stop("definition gene不是497个唯一symbol")

  pg <- fread(p$program_genes)
  clusters <- fread(p$program_clusters)
  pg[, Gene := norm_gene(Gene)]
  pg <- pg[Definition_Overlap == FALSE & !is.na(Primary_Cluster_ID) & Primary_Cluster_ID != "OTHER"]
  program_ids <- clusters[State %in% c("C1", "C4") & Cluster_ID != "OTHER", Cluster_ID]
  program_sets <- lapply(program_ids, function(id) unique(pg[Primary_Cluster_ID == id, Gene]))
  names(program_sets) <- program_ids
  if (length(program_sets) != 9L || any(lengths(program_sets) < 5L)) stop("9个RNA程序未完整冻结")

  rppa_dt <- fread(p$rppa)
  rppa_id_col <- names(rppa_dt)[1]
  rppa_raw <- as.matrix(rppa_dt[, -1, with = FALSE])
  storage.mode(rppa_raw) <- "numeric"
  rownames(rppa_raw) <- as.character(rppa_dt[[rppa_id_col]])
  colnames(rppa_raw) <- norm_id(colnames(rppa_raw))
  rm(rppa_dt)

  rppa_dict <- fread(p$rppa_dict)
  tcpa <- fread(p$tcpa)
  # TCPA旧表省略磷酸化位点分隔符，按字母数字规范键与冻结字典连接。
  map <- rppa_dict[, .(
    raw_feature,
    Normalized_Token = gsub("[^A-Z0-9]", "", toupper(rppa_key))
  )]
  if (anyDuplicated(map$Normalized_Token)) stop("RPPA规范token不唯一")
  tcpa[, Normalized_Token := gsub("[^A-Z0-9]", "", toupper(Source_Protein_Token))]
  tcpa_map <- merge(tcpa, map, by = "Normalized_Token", all.x = TRUE)
  unmapped_tcpa <- tcpa_map[is.na(raw_feature), .(
    Pathway, Source_Protein_Token, Direction_Sign,
    Reason = "NOT_MEASURED_IN_CURRENT_RPPA_PANEL"
  )]
  tcpa_map <- tcpa_map[!is.na(raw_feature)]
  if (nrow(unmapped_tcpa) != 1L || unmapped_tcpa$Source_Protein_Token != "HER3PY1298") {
    stop("TCPA未映射组件超出冻结的HER3_pY1298平台缺测边界")
  }
  tcpa_sets <- lapply(split(tcpa_map, tcpa_map$Pathway), function(d) {
    x <- d$Direction_Sign
    names(x) <- d$raw_feature
    x[!duplicated(names(x))]
  })
  if (length(tcpa_sets) != 10L || any(lengths(tcpa_sets) < 2L)) {
    stop("TCPA pathway不是10条或存在少于2个可用组件的通路")
  }

  tmb <- fread(p$tmb)
  cnv <- fread(p$cnv)
  mut <- fread(p$mutation, select = c("Patient_ID"))
  mut_count <- mut[, .(N_Nonsilent_Mutated_Genes = .N), by = Patient_ID]
  rm(mut)

  est <- fread(p$estimate)
  cib <- fread(p$cibersort)
  mcp <- fread(p$mcpcounter)
  drug <- fread(p$drug)
  est[, Sample15 := norm_id(Sample)]
  cib[, Sample15 := norm_id(Sample)]
  mcp[, Sample15 := norm_id(Sample)]

  est_cols <- intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"), names(est))
  cib_cols <- setdiff(names(cib), c("Sample", "Sample15", "Cancer", "Subtype", "P-value", "Correlation", "RMSE"))
  mcp_cols <- setdiff(names(mcp), c("Sample", "Sample15", "Cancer", "Subtype"))
  drug_cols <- c("corrected_prediction_row_mean", "corrected_prediction_row_median", "corrected_prediction_PC1")

  env <- new.env(parent = emptyenv())
  load(p$mrna, envir = env)
  if (!exists("pancancer_data_list", envir = env, inherits = FALSE)) stop("mRNA对象缺少pancancer_data_list")
  mrna_list <- env$pancancer_data_list

  manifest_rows <- list()
  program_registry_rows <- list()

  for (cancer in sort(unique(patient$Cancer))) {
    meta <- copy(patient[Cancer == cancer])
    setorder(meta, Patient_ID)
    ids <- meta$Patient_ID
    samples <- norm_id(meta$Sample15)
    if (!cancer %in% names(mrna_list)) stop("mRNA缺少癌种: ", cancer)

    x <- mrna_list[[cancer]]
    expr <- x$expression_data
    rownames(expr) <- norm_gene(rownames(expr))
    colnames(expr) <- norm_id(colnames(expr))
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
    midx <- match(samples, colnames(expr))
    if (anyNA(midx)) stop(cancer, " canonical样本无法完整映射mRNA")
    log_expr <- log2(expr[, midx, drop = FALSE] + 1)
    colnames(log_expr) <- ids

    clin <- as.data.table(x$clinical_data)
    clin[, Sample15 := norm_id(sample)]
    ci <- match(samples, clin$Sample15)
    meta[, Age := suppressWarnings(as.numeric(clin$age_at_initial_pathologic_diagnosis[ci]))]
    meta[, Sex := toupper(as.character(clin$gender[ci]))]
    meta[!Sex %in% c("MALE", "FEMALE"), Sex := NA_character_]

    prog <- matrix(NA_real_, nrow = nrow(meta), ncol = length(program_sets),
                   dimnames = list(ids, names(program_sets)))
    for (program in names(program_sets)) {
      genes <- intersect(program_sets[[program]], rownames(log_expr))
      if (length(genes) < 5L) stop(cancer, " 程序可用基因少于5: ", program)
      zz <- z_rows(log_expr[genes, , drop = FALSE])
      prog[, program] <- colMeans(zz, na.rm = TRUE)
      program_registry_rows[[length(program_registry_rows) + 1L]] <- data.table(
        Cancer = cancer, Program = program, N_Frozen_Nondefinition = length(program_sets[[program]]),
        N_Mapped = length(genes), Genes_Mapped = paste(sort(genes), collapse = ";")
      )
    }

    ridx <- match(samples, colnames(rppa_raw))
    rppa <- t(rppa_raw[, ridx, drop = FALSE])
    rownames(rppa) <- ids
    colnames(rppa) <- rownames(rppa_raw)
    rppa_feature_by_sample <- t(rppa)
    rppa_path <- score_sets(rppa_feature_by_sample, tcpa_sets, min_features = 2L)

    cnv_cols <- c("mean_abs_continuous", "fraction_abs_ge1", "fraction_gain_ge1",
                  "fraction_loss_le_minus1", "fraction_amp_plus2", "fraction_deepdel_minus2")
    cnv[, Patient_ID := as.character(Patient_ID)]
    cnv_mat <- left_matrix(ids, cnv, "Patient_ID", cnv_cols, "CNV_")

    tmb_cols <- c("TMB_Log10", "MANTIS_Score")
    tmb_mat <- left_matrix(ids, tmb, "Patient_ID", tmb_cols, "GENOME_")
    mc_idx <- match(ids, mut_count$Patient_ID)
    mut_log <- log1p(ifelse(is.na(mc_idx), 0, mut_count$N_Nonsilent_Mutated_Genes[mc_idx]))
    genome <- cbind(tmb_mat, GENOME_Log1p_Nonsilent_Mutated_Genes = mut_log)
    rownames(genome) <- ids

    est_mat <- left_matrix(samples, est, "Sample15", est_cols, "ESTIMATE_")
    cib_mat <- left_matrix(samples, cib, "Sample15", cib_cols, "CIBERSORT_")
    mcp_mat <- left_matrix(samples, mcp, "Sample15", mcp_cols, "MCPCOUNTER_")
    rownames(est_mat) <- rownames(cib_mat) <- rownames(mcp_mat) <- ids
    immune <- cbind(est_mat, cib_mat, mcp_mat)

    di <- match(ids, drug$Patient_ID)
    drug_mat <- as.matrix(drug[di, ..drug_cols])
    storage.mode(drug_mat) <- "numeric"
    rownames(drug_mat) <- ids
    colnames(drug_mat) <- paste0("PRED_DRUG_", drug_cols)

    strict_ids <- meta[MOFA_All_Orthogonal == TRUE, Patient_ID]
    topvar <- matrix(numeric(), nrow = length(strict_ids), ncol = 0,
                     dimnames = list(strict_ids, character()))
    rppa_mofa <- matrix(numeric(), nrow = length(strict_ids), ncol = 0,
                        dimnames = list(strict_ids, character()))
    if (cancer %in% pilot_cancers) {
      nondef <- setdiff(rownames(log_expr), definition_genes)
      mm <- log_expr[nondef, match(strict_ids, ids), drop = FALSE]
      vv <- rowVars(mm, na.rm = TRUE)
      valid <- is.finite(vv) & vv > 0 & rowMeans(is.finite(mm)) >= 0.95
      ord <- order(vv[valid], decreasing = TRUE)
      selected_genes <- rownames(mm)[valid][head(ord, 200L)]
      topvar <- t(z_rows(mm[selected_genes, , drop = FALSE]))
      rownames(topvar) <- strict_ids

      rr <- rppa[match(strict_ids, ids), , drop = FALSE]
      obs <- colMeans(is.finite(rr))
      vv_r <- colVars(rr, na.rm = TRUE)
      keep_r <- is.finite(vv_r) & vv_r > 0 & obs >= 0.80
      rppa_mofa <- t(z_rows(t(rr[, keep_r, drop = FALSE])))
      rownames(rppa_mofa) <- strict_ids
    }

    bundle <- list(
      cancer = cancer,
      patient = as.data.frame(meta),
      rna_program = prog,
      rna_topvar_200 = topvar,
      rppa_all = rppa,
      rppa_mofa = rppa_mofa,
      rppa_tcpa_pathway = rppa_path,
      cnv = cnv_mat,
      genomic = genome,
      immune = immune,
      predicted_drug = drug_mat,
      provenance = p
    )
    out_file <- file.path(bundle_dir, paste0(cancer, "_stage1_numeric_bundle_v0.1_20260810.rds"))
    saveRDS(bundle, out_file, compress = "xz")

    manifest_rows[[cancer]] <- data.table(
      Cancer = cancer, N_Canonical = nrow(meta), N_C1 = sum(meta$C_subtype == "C1"),
      N_C4 = sum(meta$C_subtype == "C4"), N_Strict_MOFA = length(strict_ids),
      Is_MOFA_Pilot = cancer %in% pilot_cancers,
      N_RNA_Programs = ncol(prog), N_RNA_TopVar = ncol(topvar),
      N_RPPA_All = ncol(rppa), N_RPPA_MOFA = ncol(rppa_mofa),
      N_TCPA_Pathways = ncol(rppa_path), N_CNV = ncol(cnv_mat),
      N_Genomic = ncol(genome), N_Immune = ncol(immune), N_PredDrug = ncol(drug_mat),
      Bundle = out_file, Bundle_Bytes = file.info(out_file)$size
    )
    log_event("CHECKPOINT", "癌种数值bundle完成", list(
      cancer = cancer, n_samples = nrow(meta), n_strict = length(strict_ids),
      n_programs = ncol(prog), n_rppa = ncol(rppa_mofa)
    ))
  }

  manifest <- rbindlist(manifest_rows)
  program_registry <- rbindlist(program_registry_rows)
  fwrite(manifest, file.path(out_dir, "stage1_numeric_bundle_manifest.tsv"), sep = "\t")
  fwrite(program_registry, file.path(out_dir, "stage1_rna_program_mapping_by_cancer.tsv"), sep = "\t")
  fwrite(unmapped_tcpa, file.path(out_dir, "stage1_tcpa_unmapped_component.tsv"), sep = "\t")

  if (nrow(manifest) != 22L || any(manifest$N_RNA_Programs != 9L) ||
      any(manifest[Is_MOFA_Pilot == TRUE, N_RNA_TopVar] != 200L) ||
      any(manifest[Is_MOFA_Pilot == TRUE, N_RPPA_MOFA] < 100L) ||
      any(manifest$Bundle_Bytes <= 0)) stop("数值bundle QA失败")

  report <- c(
    "# Stage 1共用数值输入构建报告 V0.1", "",
    "日期：2026-08-10  ",
    "状态：`PASS_SHARED_NUMERIC_INPUT_GATE / TRACK_A_READY / MOFA_INPUT_READY_ENGINE_MISSING`", "",
    sprintf("已为22个癌种生成独立bundle；每癌固定9个去定义RNA程序、10个TCPA蛋白通路、6个CNV指标、3个基因组指标和%s个免疫注释。", unique(manifest$N_Immune)),
    sprintf("5个MOFA试点癌种各生成200个无监督高变非定义RNA特征；合格RPPA feature范围为%s–%s。", min(manifest[Is_MOFA_Pilot == TRUE, N_RPPA_MOFA]), max(manifest[Is_MOFA_Pilot == TRUE, N_RPPA_MOFA])),
    "官方TCPA表中HER3_pY1298在当前RPPA平台缺测，已登记并从RTK pathway score中排除；其余通路组件门槛满足。",
    "所有对象按癌种和患者构建，不合并不同癌种；免疫和预测药敏未参与MOFA特征选择。", "",
    "当前R4.5.1缺少MOFA2，因此输入已冻结但未训练MOFA。线路A可直接进入正式关联。"
  )
  writeLines(report, file.path(out_dir, "stage1_shared_numeric_input_report_v0.1.md"), useBytes = TRUE)

  log_event("DONE", "Stage 1共用数值输入完成", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    cancers = nrow(manifest), pilot_cancers = paste(pilot_cancers, collapse = ","),
    output = out_dir
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1L)
})
