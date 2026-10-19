#!/usr/bin/env Rscript
# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Task: single_cell_pseudobulk; step: s08
# Inputs (DCF fields): technical_gate,model_contract,counts_manifest,exclusion_genes,counts_dir,sample_sets_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 08_run_primary_pseudobulk_models.R --config PATH [--check-config]\n       Rscript 08_run_primary_pseudobulk_models.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("technical_gate", "model_contract", "counts_manifest", "exclusion_genes", "counts_dir", "sample_sets_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("technical_gate", "model_contract", "counts_manifest", "exclusion_genes")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir", "sample_sets_dir")])
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

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02c_run_primary_pseudobulk_models.R"

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

select_one_per_patient <- function(dat) {
  dd <- copy(dat)
  dd[, pre_priority := fifelse(pre_post == "Pre", 0L, 1L)]
  setorder(dd, patient_key, pre_priority, -n_cells, -library_size, sample_unit_key)
  dd[, .SD[1L], by = patient_key]
}

make_design <- function(dat, group1, group2, include_pre_post = TRUE) {
  dd <- copy(dat)
  # sample set 中也有 group1/group2 列，使用局部别名避免 data.table 同名解析。
  g1 <- as.character(group1)[1]
  g2 <- as.character(group2)[1]
  dd[, group := factor(analysis_group, levels = c(g2, g1))]
  dd[, cohort_block := droplevels(factor(cohort_block))]
  use_cohort <- uniqueN(dd$cohort_block) > 1L
  use_pre_post <- include_pre_post && uniqueN(dd$pre_post, na.rm = TRUE) > 1L
  if (use_pre_post) dd[, pre_post := factor(pre_post, levels = c("Pre", "Post"))]

  rhs <- c("0 + group", if (use_cohort) "cohort_block", if (use_pre_post) "pre_post")
  formula <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  mm <- model.matrix(formula, data = dd)
  colnames(mm) <- make.names(colnames(mm))

  # pre/post 在稀疏子集中可能与 cohort 完全共线；只降级这一非研究目标协变量。
  if (qr(mm)$rank < ncol(mm) && use_pre_post) {
    use_pre_post <- FALSE
    rhs <- c("0 + group", if (use_cohort) "cohort_block")
    formula <- as.formula(paste("~", paste(rhs, collapse = " + ")))
    mm <- model.matrix(formula, data = dd)
    colnames(mm) <- make.names(colnames(mm))
  }
  if (qr(mm)$rank < ncol(mm)) stop("Design matrix remains rank deficient")
  list(data = dd, design = mm, formula = paste(deparse(formula), collapse = ""),
       uses_cohort = use_cohort, uses_pre_post = use_pre_post)
}

fit_voom_model <- function(counts, dat, group1, group2, block_patient = FALSE,
                           consensus_gene_cap = 200L) {
  g1 <- as.character(group1)[1]
  g2 <- as.character(group2)[1]
  idx <- match(dat$unit_id, colnames(counts))
  if (anyNA(idx)) stop("Selected sample IDs do not match pseudobulk counts")
  payload <- make_design(dat, g1, g2, include_pre_post = TRUE)
  dat <- payload$data
  design <- payload$design
  y <- DGEList(counts = counts[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = design)
  if (sum(keep) < 100L) stop("Fewer than 100 genes passed filterByExpr")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)

  consensus <- NA_real_
  correlation_genes <- 0L
  if (block_patient) {
    block <- factor(dat$patient_key)
    mean_log_cpm <- rowMeans(cpm(y, log = TRUE, prior.count = 1))
    correlation_genes <- min(consensus_gene_cap, nrow(y))
    correlation_idx <- order(mean_log_cpm, decreasing = TRUE)[seq_len(correlation_genes)]
    # consensus correlation 是技术参数；使用高表达代表集估计并在全基因模型中固定。
    corfit <- duplicateCorrelation(v[correlation_idx, ], design, block = block)
    consensus <- corfit$consensus.correlation
    fit <- lmFit(v, design, block = block, correlation = consensus)
  } else {
    fit <- lmFit(v, design)
  }

  contrast_text <- paste0("group", make.names(g1), "-group", make.names(g2))
  contrast_matrix <- makeContrasts(contrasts = contrast_text, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
  tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"), keep.rownames = "gene")
  setnames(tt, c("P.Value", "adj.P.Val"), c("p_value", "fdr"))
  tt[, direction := fifelse(logFC > 0, paste0(g1, "_higher"),
                            fifelse(logFC < 0, paste0(g2, "_higher"), "zero"))]
  list(
    result = tt,
    summary = data.table(
      n_samples = nrow(dat), n_patients = uniqueN(dat$patient_key),
      n_group1_samples = sum(dat$analysis_group == g1),
      n_group2_samples = sum(dat$analysis_group == g2),
      n_group1_patients = uniqueN(dat[analysis_group == g1, patient_key]),
      n_group2_patients = uniqueN(dat[analysis_group == g2, patient_key]),
      n_genes_tested = nrow(tt), design_columns = ncol(design), design_rank = qr(design)$rank,
      formula = payload$formula, uses_cohort = payload$uses_cohort,
      uses_pre_post = payload$uses_pre_post, patient_block_used = block_patient,
      consensus_correlation = consensus, correlation_genes = correlation_genes,
      n_fdr_005 = sum(tt$fdr < 0.05, na.rm = TRUE),
      n_fdr_010 = sum(tt$fdr < 0.10, na.rm = TRUE),
      median_abs_logFC = median(abs(tt$logFC), na.rm = TRUE)
    )
  )
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  counts_dir <- unname(paths["counts_dir"])
  sample_sets_dir <- unname(paths["sample_sets_dir"])
  result_dir <- file.path(outdir, "gene_results")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  gate_path <- unname(paths["technical_gate"])
  contract_path <- unname(paths["model_contract"])
  manifest_path <- unname(paths["counts_manifest"])
  exclusion_path <- unname(paths["exclusion_genes"])
  required <- c(gate_path, contract_path, manifest_path, exclusion_path)
  if (!all(file.exists(required))) stop("Required G2A inputs are incomplete")
  gate <- fread(gate_path)
  if (gate$status[[1]] != "PASS") stop("G2A technical gate is not PASS")

  fwrite(data.table(
    gate = "G2B_primary_overall_models", status = "RUNNING",
    n_main_pass = NA_integer_, n_main_total = 12L,
    n_sensitivity_pass = NA_integer_, n_sensitivity_total = 48L,
    note = "RUNNING prevents stale biological-model status after interruption."
  ), file.path(outdir, "g2b_overall_model_gate.tsv"), sep = "\t", na = "NA")

  contract <- fread(contract_path)
  manifest <- fread(manifest_path)
  exclusion <- fread(exclusion_path)[exclude_from_independent_discovery == TRUE, unique(gene)]
  core <- c("Malignant", "T/NK", "Myeloid", "Stromal/CAF")
  targets <- contract[role == "primary" & compartment %in% core & status == "PASS_TECHNICAL_DESIGN"]
  if (nrow(targets) != 12L) stop("Expected 12 primary core models")

  main_summaries <- list()
  sensitivity_summaries <- list()
  run_status <- list()

  log_event("START", "run formal primary pseudobulk overall models", output = outdir, n_targets = nrow(targets))
  for (i in seq_len(nrow(targets))) {
    target <- targets[i]
    prefix <- paste0(safe_name(target$comparison), "__", safe_name(target$compartment))
    log_event("CHECKPOINT", "fit target", comparison = target$comparison,
              compartment = target$compartment, strategy = target$recommended_model)
    tryCatch({
      count_file <- manifest[compartment == target$compartment, file][1]
      se <- readRDS(relocate_input(count_file, counts_dir))
      counts <- assay(se, "counts")
      sample_file <- relocate_input(paste0(prefix, ".tsv.gz"), sample_sets_dir)
      dat_all <- fread(sample_file)
      dat_all <- dat_all[n_cells >= target$threshold]

      if (grepl("duplicateCorrelation", target$recommended_model)) {
        main_dat <- copy(dat_all)
        main_block <- TRUE
        main_label <- "all_samples_patient_block"
      } else if (grepl("no_patient_block", target$recommended_model)) {
        main_dat <- copy(dat_all)
        main_block <- FALSE
        main_label <- "all_samples_no_patient_block"
      } else {
        main_dat <- select_one_per_patient(dat_all)
        main_block <- FALSE
        main_label <- "independent_patient"
      }

      main_fit <- fit_voom_model(counts, main_dat, target$group1, target$group2,
                                 block_patient = main_block, consensus_gene_cap = 200L)
      main_res <- main_fit$result
      main_res[, `:=`(
        comparison = target$comparison, compartment = target$compartment,
        analysis = "main", model_label = main_label,
        definition_gene = gene %in% exclusion
      )]
      setcolorder(main_res, c("comparison", "compartment", "analysis", "model_label", "gene",
                              "definition_gene", "logFC", "AveExpr", "t", "p_value", "fdr", "B", "direction"))
      fwrite(main_res, file.path(result_dir, paste0(prefix, "__main.tsv.gz")),
             sep = "\t", na = "NA", compress = "gzip")
      main_summaries[[length(main_summaries) + 1L]] <- cbind(
        data.table(comparison = target$comparison, compartment = target$compartment,
                   model_label = main_label, recommended_model = target$recommended_model,
                   status = "PASS"),
        main_fit$summary
      )

      sensitivity_defs <- list(
        independent_patient = list(threshold = 20L, pre_only = FALSE),
        pre_only = list(threshold = 20L, pre_only = TRUE),
        independent_patient_ge30 = list(threshold = 30L, pre_only = FALSE),
        independent_patient_ge50 = list(threshold = 50L, pre_only = FALSE)
      )
      for (sname in names(sensitivity_defs)) {
        spec <- sensitivity_defs[[sname]]
        sdat <- dat_all[n_cells >= spec$threshold]
        if (spec$pre_only) sdat <- sdat[pre_post == "Pre"]
        sdat <- select_one_per_patient(sdat)
        if (uniqueN(sdat$analysis_group) < 2L ||
            min(sdat[, uniqueN(patient_key), by = analysis_group]$V1) < 2L) {
          sensitivity_summaries[[length(sensitivity_summaries) + 1L]] <- data.table(
            comparison = target$comparison, compartment = target$compartment,
            analysis = sname, status = "NOT_EVALUABLE", reason = "fewer_than_2_patients_in_group"
          )
          next
        }
        sfit <- fit_voom_model(counts, sdat, target$group1, target$group2, block_patient = FALSE)
        sres <- sfit$result
        sres[, `:=`(
          comparison = target$comparison, compartment = target$compartment,
          analysis = sname, model_label = "independent_patient",
          definition_gene = gene %in% exclusion
        )]
        setcolorder(sres, c("comparison", "compartment", "analysis", "model_label", "gene",
                            "definition_gene", "logFC", "AveExpr", "t", "p_value", "fdr", "B", "direction"))
        fwrite(sres, file.path(result_dir, paste0(prefix, "__", sname, ".tsv.gz")),
               sep = "\t", na = "NA", compress = "gzip")
        sensitivity_summaries[[length(sensitivity_summaries) + 1L]] <- cbind(
          data.table(comparison = target$comparison, compartment = target$compartment,
                     analysis = sname, status = "PASS", reason = NA_character_),
          sfit$summary
        )
      }
      run_status[[length(run_status) + 1L]] <- data.table(
        comparison = target$comparison, compartment = target$compartment, status = "PASS", error = NA_character_)
      rm(se, counts, main_fit, main_res)
      gc(verbose = FALSE)
    }, error = function(e) {
      log_event("ERROR", "target failed", comparison = target$comparison,
                compartment = target$compartment, error = conditionMessage(e))
      run_status[[length(run_status) + 1L]] <<- data.table(
        comparison = target$comparison, compartment = target$compartment,
        status = "FAIL", error = conditionMessage(e))
    })
  }

  main_summary <- rbindlist(main_summaries, fill = TRUE)
  sensitivity_summary <- rbindlist(sensitivity_summaries, fill = TRUE)
  status_table <- rbindlist(run_status, fill = TRUE)
  fwrite(main_summary, file.path(outdir, "main_model_summary.tsv"), sep = "\t", na = "NA")
  fwrite(sensitivity_summary, file.path(outdir, "sensitivity_model_summary.tsv"), sep = "\t", na = "NA")
  fwrite(status_table, file.path(outdir, "target_run_status.tsv"), sep = "\t", na = "NA")

  n_main_pass <- sum(status_table$status == "PASS")
  n_sens_pass <- sum(sensitivity_summary$status == "PASS")
  n_sens_total <- nrow(sensitivity_summary)
  overall_pass <- n_main_pass == 12L && n_sens_total == 48L &&
    all(sensitivity_summary$status %in% c("PASS", "NOT_EVALUABLE"))
  out_gate <- data.table(
    gate = "G2B_primary_overall_models",
    status = if (overall_pass) "PASS" else "FAIL",
    n_main_pass = n_main_pass, n_main_total = 12L,
    n_sensitivity_pass = n_sens_pass, n_sensitivity_total = n_sens_total,
    note = "Overall and sample-selection models only; cohort direction, LOSO and legacy survival audit remain required for G2."
  )
  fwrite(out_gate, file.path(outdir, "g2b_overall_model_gate.tsv"), sep = "\t", na = "NA")

  manifest <- data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version",
             "edgeR_version", "limma_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string,
              as.character(packageVersion("edgeR")), as.character(packageVersion("limma")), outdir)
  )
  fwrite(manifest, file.path(outdir, "02c_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "formal primary pseudobulk overall models completed",
            runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
            status = out_gate$status[[1]], output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
