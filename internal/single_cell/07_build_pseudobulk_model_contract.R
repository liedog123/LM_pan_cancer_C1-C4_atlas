#!/usr/bin/env Rscript
# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Task: single_cell_pseudobulk; step: s07
# Inputs (DCF fields): sample_sheet,counts_manifest,counts_gate,counts_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 07_build_pseudobulk_model_contract.R --config PATH [--check-config]\n       Rscript 07_build_pseudobulk_model_contract.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("sample_sheet", "counts_manifest", "counts_gate", "counts_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet", "counts_manifest", "counts_gate")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir")])
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),input_dirs,config_path,package_root)
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

# 普通执行仅在原读取位置解析选中产物；元信息检查不枚举输入目录。
relocate_input <- function(recorded, input_dir) {
  if(length(recorded)!=1L || is.na(recorded) || !nzchar(recorded)) stop("Invalid recorded input path.")
  leaf <- basename(recorded)
  if(leaf %in% c(".","..") || !nzchar(leaf)) stop("Invalid input basename.")
  candidate <- resolve_path(file.path(input_dir,leaf), dirname(config_path))
  if(!is_within(candidate,input_dir) || identical(candidate,input_dir)) stop("Input escapes configured directory.")
  if(!file_test("-f",candidate) || file.access(candidate,4L)!=0L || is.na(file.info(candidate)$size) || file.info(candidate)$size<=0) stop("Selected input must be a readable nonempty regular file.")
  candidate
}


Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02b_build_pseudobulk_model_contract.R"

log_event <- function(level, msg, ...) {
  extras <- list(...)
  extra_text <- if (length(extras) == 0L) "" else paste(
    vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)),
    collapse = " "
  )
  cat(sprintf("[%s] %s script=%s msg=%s%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg),
              if (nzchar(extra_text)) " " else "", extra_text))
  flush.console()
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

make_design <- function(dat, include_pre_post = TRUE) {
  dd <- copy(dat)
  dd[, group := factor(analysis_group, levels = c(group2[1], group1[1]))]
  dd[, cohort_block := droplevels(factor(cohort_block))]
  use_pre_post <- include_pre_post && uniqueN(dd$pre_post, na.rm = TRUE) > 1L
  if (use_pre_post) dd[, pre_post := factor(pre_post, levels = c("Pre", "Post"))]
  formula <- if (use_pre_post) ~ 0 + group + cohort_block + pre_post else ~ 0 + group + cohort_block
  mm <- model.matrix(formula, data = dd)
  colnames(mm) <- make.names(colnames(mm))
  list(
    data = dd,
    matrix = mm,
    formula = paste(deparse(formula), collapse = ""),
    rank = qr(mm)$rank,
    n_columns = ncol(mm),
    full_rank = qr(mm)$rank == ncol(mm),
    uses_pre_post = use_pre_post
  )
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  counts_dir <- unname(paths["counts_dir"])
  design_dir <- file.path(outdir, "design_matrices")
  sample_dir <- file.path(outdir, "sample_sets")
  dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  sample_path <- unname(paths["sample_sheet"])
  manifest_path <- unname(paths["counts_manifest"])
  count_gate_path <- unname(paths["counts_gate"])
  if (!all(file.exists(c(sample_path, manifest_path, count_gate_path)))) stop("Pseudobulk count outputs are incomplete")

  log_event("START", "build pseudobulk model and sensitivity contract", input = sample_path, output = outdir)
  fwrite(data.table(
    gate = "G2A_pseudobulk_technical_contract", status = "RUNNING", counts_gate = NA_character_,
    n_primary_core_designs_pass = NA_integer_, n_primary_core_designs_total = NA_integer_,
    smoke_test_status = "RUNNING",
    note = "RUNNING prevents a stale PASS from a prior interrupted execution."
  ), file.path(outdir, "g2a_technical_gate_summary.tsv"), sep = "\t", na = "NA")
  sample_sheet <- fread(sample_path)
  manifest <- fread(manifest_path)
  count_gate <- fread(count_gate_path)
  if (!identical(count_gate$status[[1]], "PASS")) stop("Pseudobulk count gate is not PASS")

  comparisons <- data.table(
    comparison = c("C1_vs_C4", "AG2_IG1_vs_AG2_IG2", "AG2_IG2_vs_AG3_IG2",
                   "C1_vs_C2", "C2_vs_C4", "C3_vs_C2"),
    field = c("sample_C_subtype", "sample_AGIG", "sample_AGIG",
              "sample_C_subtype", "sample_C_subtype", "sample_C_subtype"),
    group1 = c("C1", "AG2_IG1", "AG2_IG2", "C1", "C2", "C3"),
    group2 = c("C4", "AG2_IG2", "AG3_IG2", "C2", "C4", "C2"),
    role = c("primary", "primary", "primary", "secondary", "secondary", "secondary")
  )
  compartments <- manifest$compartment
  standard_threshold <- 20L
  contract_rows <- list()
  smoke_payloads <- list()

  log_event("CHECKPOINT", "build comparison-by-compartment design contracts")
  for (i in seq_len(nrow(comparisons))) {
    con <- comparisons[i]
    for (compartment in compartments) {
      dat <- sample_sheet[project_broad_compartment == compartment & n_cells >= standard_threshold]
      dat[, analysis_group := as.character(get(con$field))]
      dat <- dat[analysis_group %in% c(con$group1, con$group2)]
      dat[, `:=`(group1 = con$group1, group2 = con$group2)]

      n1 <- uniqueN(dat[analysis_group == con$group1, sample_unit_key])
      n2 <- uniqueN(dat[analysis_group == con$group2, sample_unit_key])
      p1 <- uniqueN(dat[analysis_group == con$group1, patient_key])
      p2 <- uniqueN(dat[analysis_group == con$group2, patient_key])
      patient_counts <- dat[, .N, by = patient_key]
      n_repeated_patients <- sum(patient_counts$N > 1L)
      patient_groups <- dat[, .(n_groups = uniqueN(analysis_group)), by = patient_key]
      n_switching_patients <- sum(patient_groups$n_groups > 1L)

      base_design <- adjusted_design <- independent_design <- pre_only_design <- NULL
      status <- "NOT_EVALUABLE"
      reason <- NA_character_
      model_strategy <- "descriptive_only"
      selected_design <- NULL
      selected_dat <- NULL

      if (n1 == 0L || n2 == 0L) {
        reason <- "group_missing"
      } else if (p1 < 2L || p2 < 2L) {
        reason <- "fewer_than_2_patients_in_group"
      } else {
        base_design <- make_design(dat, include_pre_post = FALSE)
        adjusted_design <- make_design(dat, include_pre_post = TRUE)

        selected <- copy(dat)
        selected[, pre_priority := fifelse(pre_post == "Pre", 0L, 1L)]
        setorder(selected, patient_key, pre_priority, -n_cells, -library_size, sample_unit_key)
        selected[, independent_patient_selected := seq_len(.N) == 1L, by = patient_key]
        dat <- merge(dat, selected[, .(unit_id, independent_patient_selected)], by = "unit_id", all.x = TRUE)
        dat[, pre_only_selected := pre_post == "Pre"]
        independent_dat <- dat[independent_patient_selected == TRUE]
        pre_only_dat <- dat[pre_only_selected == TRUE]
        if (uniqueN(independent_dat$analysis_group) == 2L) independent_design <- make_design(independent_dat, include_pre_post = TRUE)
        if (uniqueN(pre_only_dat$analysis_group) == 2L) pre_only_design <- make_design(pre_only_dat, include_pre_post = FALSE)

        if (!base_design$full_rank) {
          reason <- "base_design_rank_deficient"
        } else {
          status <- "PASS_TECHNICAL_DESIGN"
          reason <- "design_identifiable"
          if (n_repeated_patients >= 5L && adjusted_design$full_rank) {
            model_strategy <- "voom_duplicateCorrelation_cohort_prepost_patient_block"
            selected_design <- adjusted_design
            selected_dat <- adjusted_design$data
          } else if (n_repeated_patients == 0L && adjusted_design$full_rank) {
            model_strategy <- "voom_cohort_prepost_no_patient_block"
            selected_design <- adjusted_design
            selected_dat <- adjusted_design$data
          } else if (!is.null(independent_design) && independent_design$full_rank) {
            model_strategy <- "voom_independent_patient_primary_preonly_sensitivity"
            selected_design <- independent_design
            selected_dat <- independent_design$data
          } else {
            model_strategy <- "voom_cohort_base_with_patient_limitation"
            selected_design <- base_design
            selected_dat <- base_design$data
          }
        }

        dat[, comparison := con$comparison]
        dat[, role := con$role]
        dat[, standard_threshold := standard_threshold]
        sample_file <- file.path(sample_dir, paste0(safe_name(con$comparison), "__", safe_name(compartment), ".tsv.gz"))
        fwrite(dat, sample_file, sep = "\t", na = "NA", compress = "gzip")

        if (!is.null(selected_design)) {
          design_dt <- as.data.table(selected_design$matrix)
          design_dt[, unit_id := selected_dat$unit_id]
          setcolorder(design_dt, "unit_id")
          design_file <- file.path(design_dir, paste0(safe_name(con$comparison), "__", safe_name(compartment), "__primary_design.tsv.gz"))
          fwrite(design_dt, design_file, sep = "\t", na = "NA", compress = "gzip")
          if (compartment == "T/NK" && con$comparison %in% c("C1_vs_C4", "AG2_IG1_vs_AG2_IG2")) {
            smoke_payloads[[con$comparison]] <- list(
              comparison = con$comparison, compartment = compartment,
              group1 = con$group1, group2 = con$group2,
              strategy = model_strategy, data = selected_dat,
              design = selected_design$matrix
            )
          }
        }
      }

      contract_rows[[length(contract_rows) + 1L]] <- data.table(
        comparison = con$comparison, role = con$role, field = con$field,
        group1 = con$group1, group2 = con$group2,
        compartment = compartment, threshold = standard_threshold,
        n_samples_group1 = n1, n_samples_group2 = n2,
        n_patients_group1 = p1, n_patients_group2 = p2,
        n_repeated_patients = n_repeated_patients,
        n_patients_switching_groups = n_switching_patients,
        n_pre_samples = sum(dat$pre_post == "Pre", na.rm = TRUE),
        n_post_samples = sum(dat$pre_post == "Post", na.rm = TRUE),
        base_formula = if (is.null(base_design)) NA_character_ else base_design$formula,
        base_rank = if (is.null(base_design)) NA_integer_ else base_design$rank,
        base_columns = if (is.null(base_design)) NA_integer_ else base_design$n_columns,
        base_full_rank = if (is.null(base_design)) FALSE else base_design$full_rank,
        adjusted_formula = if (is.null(adjusted_design)) NA_character_ else adjusted_design$formula,
        adjusted_rank = if (is.null(adjusted_design)) NA_integer_ else adjusted_design$rank,
        adjusted_columns = if (is.null(adjusted_design)) NA_integer_ else adjusted_design$n_columns,
        adjusted_full_rank = if (is.null(adjusted_design)) FALSE else adjusted_design$full_rank,
        independent_patient_full_rank = if (is.null(independent_design)) FALSE else independent_design$full_rank,
        pre_only_full_rank = if (is.null(pre_only_design)) FALSE else pre_only_design$full_rank,
        status = status, reason = reason, recommended_model = model_strategy,
        contrast = paste0("group", make.names(con$group1), " - group", make.names(con$group2))
      )
    }
  }
  contract <- rbindlist(contract_rows, fill = TRUE)
  fwrite(contract, file.path(outdir, "model_design_contract.tsv"), sep = "\t", na = "NA")

  covariates <- data.table(
    variable = c("target_contrast", "cohort_block", "pre_post", "patient_key", "library_size", "n_cells", "treatment_class_broad"),
    primary_role = c("required", "required", "required_when_variable_and_identifiable", "block_when_repeated_support_ge5", "TMM_voom_offset_weight", "QC_not_default_linear_covariate", "excluded_from_primary_model"),
    sensitivity_role = c("fixed", "leave_one_study_out", "pre_only", "one_sample_per_patient", "library_QC", "threshold_20_30_50", "descriptive_confounding_audit"),
    rationale = c(
      "固定亚型或原子组比较，不重新分型",
      "阻断 study 与 cancer type 的组合背景，避免机械加入共线变量",
      "56 对重复样本全部为 Pre/Post，需调整非研究目标的采样时相",
      "细胞不是重复；重复样本相关性由 duplicateCorrelation 或患者级敏感性处理",
      "由 TMM normalization 与 voom precision weight 处理",
      "细胞数用于纳入门槛和敏感性，避免过度校正潜在生物学信号",
      "治疗不是科学问题且类别与队列结构纠缠，不作为默认协变量"
    )
  )
  fwrite(covariates, file.path(outdir, "model_covariate_contract.tsv"), sep = "\t", na = "NA")

  thresholds <- data.table(
    item = c("standard_min_cells_per_unit", "sensitivity_min_cells_per_unit_1", "sensitivity_min_cells_per_unit_2",
             "gene_expression_filter", "minimum_patients_per_group_structural", "formal_evidence_requirement"),
    value = c("20", "30", "50", "edgeR::filterByExpr within each design", "2",
              "overall model plus study-specific direction and leave-one-study-out; technical PASS alone is insufficient")
  )
  fwrite(thresholds, file.path(outdir, "pseudobulk_threshold_contract.tsv"), sep = "\t", na = "NA")

  library_qc <- sample_sheet[, .(
    n_units = .N,
    n_units_ge20 = sum(n_cells >= 20L),
    n_units_ge30 = sum(n_cells >= 30L),
    n_units_ge50 = sum(n_cells >= 50L),
    median_cells = as.numeric(median(n_cells)),
    median_library_size = as.numeric(median(library_size)),
    median_detected_genes = as.numeric(median(detected_genes)),
    spearman_cells_library = as.numeric(cor(n_cells, library_size, method = "spearman")),
    spearman_cells_detected_genes = as.numeric(cor(n_cells, detected_genes, method = "spearman"))
  ), by = project_broad_compartment]
  fwrite(library_qc, file.path(outdir, "pseudobulk_library_qc_summary.tsv"), sep = "\t", na = "NA")

  log_event("CHECKPOINT", "run non-interpretive independent-patient and patient-block smoke tests")
  tnk_file <- manifest[compartment == "T/NK", file][1]
  tnk_se <- readRDS(relocate_input(tnk_file, counts_dir))
  tnk_counts <- assay(tnk_se, "counts")
  smoke <- rbindlist(lapply(smoke_payloads, function(smoke_payload) {
    tryCatch({
      ids <- smoke_payload$data$unit_id
      idx <- match(ids, colnames(tnk_counts))
      if (anyNA(idx)) stop("Smoke-test sample IDs do not match T/NK counts")
      y <- DGEList(counts = tnk_counts[, idx, drop = FALSE])
      design <- smoke_payload$design
      keep <- filterByExpr(y, design = design)
      y <- y[keep, , keep.lib.sizes = FALSE]
      y <- calcNormFactors(y)
      n_genes_after_filter <- nrow(y)
      # 这里只验证模型路径可拟合，不以冒烟测试替代正式的全基因模型。
      # 小型代表集也可避免部分服务器 BLAS 在线程/内存压力下无 R 错误退出。
      smoke_gene_cap <- if (grepl("duplicateCorrelation", smoke_payload$strategy)) 200L else 1000L
      if (is.finite(smoke_gene_cap) && nrow(y) > smoke_gene_cap) {
        mean_log_cpm <- rowMeans(cpm(y, log = TRUE, prior.count = 1))
        keep_smoke <- order(mean_log_cpm, decreasing = TRUE)[seq_len(smoke_gene_cap)]
        y <- y[keep_smoke, , keep.lib.sizes = FALSE]
      }
      v <- voom(y, design, plot = FALSE)
      block_used <- grepl("duplicateCorrelation", smoke_payload$strategy)
      consensus <- NA_real_
      if (block_used) {
        log_event("CHECKPOINT", "patient-block smoke: estimate within-patient correlation",
                  comparison = smoke_payload$comparison, genes = nrow(y), samples = ncol(y))
        block <- factor(smoke_payload$data$patient_key)
        corfit <- duplicateCorrelation(v, design, block = block)
        consensus <- corfit$consensus.correlation
        log_event("CHECKPOINT", "patient-block smoke: fit blocked model",
                  comparison = smoke_payload$comparison, consensus_correlation = signif(consensus, 4))
        fit <- lmFit(v, design, block = block, correlation = consensus)
      } else {
        fit <- lmFit(v, design)
      }
      contrast_text <- paste0("group", make.names(smoke_payload$group1), "-group", make.names(smoke_payload$group2))
      contrast_matrix <- makeContrasts(contrasts = contrast_text, levels = design)
      fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
      if (nrow(fit2$coefficients) != nrow(y)) stop("Smoke-test fitted gene count mismatch")
      data.table(
        comparison = smoke_payload$comparison,
        compartment = smoke_payload$compartment,
        status = "PASS",
        strategy = smoke_payload$strategy,
        n_samples = ncol(y), n_patients = uniqueN(smoke_payload$data$patient_key),
        n_genes_input = length(keep), n_genes_after_filter = n_genes_after_filter,
        n_genes_in_smoke_fit = nrow(y),
        design_rows = nrow(design), design_columns = ncol(design), design_rank = qr(design)$rank,
        patient_block_used = block_used, consensus_correlation = consensus,
        note = "Fit completed only to validate the pipeline; no DEG statistics were exported or interpreted."
      )
    }, error = function(e) {
      data.table(
        comparison = smoke_payload$comparison, compartment = smoke_payload$compartment, status = "FAIL",
        strategy = smoke_payload$strategy, n_samples = NA_integer_, n_patients = NA_integer_,
        n_genes_input = NA_integer_, n_genes_after_filter = NA_integer_,
        n_genes_in_smoke_fit = NA_integer_,
        design_rows = NA_integer_, design_columns = NA_integer_, design_rank = NA_integer_,
        patient_block_used = NA, consensus_correlation = NA_real_, note = conditionMessage(e)
      )
    })
  }), fill = TRUE)
  fwrite(smoke, file.path(outdir, "technical_smoke_test_summary.tsv"), sep = "\t", na = "NA")

  core <- c("T/NK", "Myeloid", "Malignant", "Stromal/CAF")
  primary_core <- contract[role == "primary" & compartment %in% core]
  technical_pass <- all(primary_core$status == "PASS_TECHNICAL_DESIGN") && nrow(smoke) == 2L && all(smoke$status == "PASS")
  gate <- data.table(
    gate = "G2A_pseudobulk_technical_contract",
    status = if (technical_pass) "PASS" else "FAIL",
    counts_gate = count_gate$status[[1]],
    n_primary_core_designs_pass = sum(primary_core$status == "PASS_TECHNICAL_DESIGN"),
    n_primary_core_designs_total = nrow(primary_core),
    smoke_test_status = if (all(smoke$status == "PASS")) "PASS" else "FAIL",
    note = "This authorizes formal pseudobulk execution planning only; it is not G2 biological evidence."
  )
  fwrite(gate, file.path(outdir, "g2a_technical_gate_summary.tsv"), sep = "\t", na = "NA")

  run_manifest <- data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version",
             "edgeR_version", "limma_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string,
              as.character(packageVersion("edgeR")), as.character(packageVersion("limma")), outdir)
  )
  fwrite(run_manifest, file.path(outdir, "02b_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pseudobulk model contract completed",
            runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
            status = gate$status[[1]], output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
