#!/usr/bin/env Rscript
# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Task: single_cell_pseudobulk; step: s09
# Inputs (DCF fields): primary_gate,model_contract,counts_manifest,counts_dir,sample_sets_dir,primary_gene_results_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 09_run_pseudobulk_stability.R --config PATH [--check-config]\n       Rscript 09_run_pseudobulk_stability.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("primary_gate", "model_contract", "counts_manifest", "counts_dir", "sample_sets_dir", "primary_gene_results_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("primary_gate", "model_contract", "counts_manifest")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir", "sample_sets_dir", "primary_gene_results_dir")])
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
script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02d_run_pseudobulk_stability.R"

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

make_design <- function(dat, group1, group2, include_cohort = TRUE, include_pre_post = TRUE) {
  dd <- copy(dat)
  g1 <- as.character(group1)[1]
  g2 <- as.character(group2)[1]
  dd[, group := factor(analysis_group, levels = c(g2, g1))]
  dd[, cohort_block := droplevels(factor(cohort_block))]
  use_cohort <- include_cohort && uniqueN(dd$cohort_block) > 1L
  use_pre_post <- include_pre_post && uniqueN(dd$pre_post, na.rm = TRUE) > 1L
  if (use_pre_post) dd[, pre_post := factor(pre_post, levels = c("Pre", "Post"))]
  rhs <- c("0 + group", if (use_cohort) "cohort_block", if (use_pre_post) "pre_post")
  formula <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  mm <- model.matrix(formula, data = dd)
  colnames(mm) <- make.names(colnames(mm))
  if (qr(mm)$rank < ncol(mm) && use_pre_post) {
    use_pre_post <- FALSE
    rhs <- c("0 + group", if (use_cohort) "cohort_block")
    formula <- as.formula(paste("~", paste(rhs, collapse = " + ")))
    mm <- model.matrix(formula, data = dd)
    colnames(mm) <- make.names(colnames(mm))
  }
  if (qr(mm)$rank < ncol(mm)) stop("rank_deficient_design")
  list(data = dd, design = mm, formula = paste(deparse(formula), collapse = ""),
       uses_cohort = use_cohort, uses_pre_post = use_pre_post)
}

fit_independent <- function(counts, dat, group1, group2, include_cohort = TRUE) {
  g1 <- as.character(group1)[1]
  g2 <- as.character(group2)[1]
  idx <- match(dat$unit_id, colnames(counts))
  if (anyNA(idx)) stop("sample_id_mismatch")
  payload <- make_design(dat, g1, g2, include_cohort = include_cohort, include_pre_post = TRUE)
  y <- DGEList(counts = counts[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = payload$design)
  if (sum(keep) < 100L) stop("too_few_expressed_genes")
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  v <- voom(y, payload$design, plot = FALSE)
  fit <- lmFit(v, payload$design)
  contrast <- makeContrasts(contrasts = paste0("group", make.names(g1), "-group", make.names(g2)),
                            levels = payload$design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"), keep.rownames = "gene")
  setnames(tt, c("P.Value", "adj.P.Val"), c("p_value", "fdr"))
  tt[, standard_error := fifelse(is.finite(t) & abs(t) > 1e-12, abs(logFC / t), NA_real_)]
  list(result = tt, formula = payload$formula, n_samples = nrow(dat),
       n_patients = uniqueN(dat$patient_key))
}

meta_gene <- function(beta, se) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  beta <- beta[ok]
  se <- se[ok]
  k <- length(beta)
  if (k < 2L) return(list(k = k, beta_fixed = NA_real_, p_fixed = NA_real_,
                         beta_random = NA_real_, p_random = NA_real_, i2 = NA_real_))
  w <- 1 / se^2
  b_fixed <- sum(w * beta) / sum(w)
  se_fixed <- sqrt(1 / sum(w))
  q <- sum(w * (beta - b_fixed)^2)
  df <- k - 1
  cval <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - df) / cval)
  wr <- 1 / (se^2 + tau2)
  b_random <- sum(wr * beta) / sum(wr)
  se_random <- sqrt(1 / sum(wr))
  list(k = k, beta_fixed = b_fixed, p_fixed = 2 * pnorm(-abs(b_fixed / se_fixed)),
       beta_random = b_random, p_random = 2 * pnorm(-abs(b_random / se_random)),
       i2 = if (q > 0) max(0, (q - df) / q) * 100 else 0)
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  primary_gene_results_dir <- unname(paths["primary_gene_results_dir"])
  counts_dir <- unname(paths["counts_dir"])
  sample_sets_dir <- unname(paths["sample_sets_dir"])
  detail_dir <- file.path(outdir, "detail")
  dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)

  gate <- fread(unname(paths["primary_gate"]))
  if (gate$status[[1]] != "PASS") stop("G2B overall-model gate is not PASS")
  fwrite(data.table(
    gate = "G2C_stability_execution", status = "RUNNING", n_targets_pass = NA_integer_, n_targets_total = 12L,
    note = "RUNNING prevents stale stability status after interruption."
  ), file.path(outdir, "g2c_stability_gate.tsv"), sep = "\t", na = "NA")

  contract <- fread(unname(paths["model_contract"]))
  manifest <- fread(unname(paths["counts_manifest"]))
  core <- c("Malignant", "T/NK", "Myeloid", "Stromal/CAF")
  targets <- contract[role == "primary" & compartment %in% core & status == "PASS_TECHNICAL_DESIGN"]
  if (nrow(targets) != 12L) stop("Expected 12 primary core targets")

  candidate_rows <- list()
  sensitivity_rows <- list()
  cohort_rows <- list()
  cohort_run_rows <- list()
  loso_rows <- list()
  loso_run_rows <- list()
  target_rows <- list()

  log_event("START", "run pseudobulk sensitivity, cohort-direction and LOSO stability", n_targets = 12L)
  for (i in seq_len(nrow(targets))) {
    target <- targets[i]
    prefix <- paste0(safe_name(target$comparison), "__", safe_name(target$compartment))
    log_event("CHECKPOINT", "stability target", comparison = target$comparison, compartment = target$compartment)
    tryCatch({
      main_file <- relocate_input(paste0(prefix, "__main.tsv.gz"), primary_gene_results_dir)
      main_res <- fread(main_file)
      eligible <- main_res[definition_gene == FALSE][order(p_value)]
      sig <- eligible[fdr < 0.10]
      selected <- unique(c(head(sig$gene, 1000L), head(eligible$gene, 200L)))
      cand <- main_res[gene %in% selected, .(comparison, compartment, gene, definition_gene,
                                             main_logFC = logFC, main_p_value = p_value, main_fdr = fdr)]
      cand[, selection_basis := fifelse(main_fdr < 0.10, "main_fdr_lt_0.10", "top200_fallback")]
      candidate_rows[[length(candidate_rows) + 1L]] <- cand

      sensitivity_names <- c("independent_patient", "pre_only",
                             "independent_patient_ge30", "independent_patient_ge50")
      sens_long <- list()
      for (sname in sensitivity_names) {
        sfile <- relocate_input(paste0(prefix, "__", sname, ".tsv.gz"), primary_gene_results_dir)
        sr <- fread(sfile)[gene %in% selected, .(gene, sensitivity_logFC = logFC,
                                                 sensitivity_p_value = p_value, sensitivity_fdr = fdr)]
        sr[, analysis := sname]
        sr <- merge(cand[, .(gene, main_logFC, main_fdr)], sr, by = "gene", all.x = TRUE)
        sr[, `:=`(
          comparison = target$comparison, compartment = target$compartment,
          direction_concordant = sign(main_logFC) == sign(sensitivity_logFC),
          retained_fdr_010 = main_fdr < 0.10 & sensitivity_fdr < 0.10
        )]
        sens_long[[length(sens_long) + 1L]] <- sr
      }
      sens_dt <- rbindlist(sens_long, fill = TRUE)
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- sens_dt

      count_file <- manifest[compartment == target$compartment, file][1]
      counts <- assay(readRDS(relocate_input(count_file, counts_dir)), "counts")
      dat_all <- fread(relocate_input(paste0(prefix, ".tsv.gz"), sample_sets_dir))[n_cells >= 20L]
      dat_independent <- select_one_per_patient(dat_all)

      # cohort_block = study × cancer；只在每组至少两个独立患者时估计方向。
      cohort_effects <- list()
      for (cohort_name in unique(dat_independent$cohort_block)) {
        cdat <- dat_independent[cohort_block == cohort_name]
        group_n <- cdat[, .(n_patients = uniqueN(patient_key)), by = analysis_group]
        evaluable <- nrow(group_n) == 2L && min(group_n$n_patients) >= 2L
        if (!evaluable) {
          cohort_run_rows[[length(cohort_run_rows) + 1L]] <- data.table(
            comparison = target$comparison, compartment = target$compartment,
            cohort_block = cohort_name, status = "NOT_EVALUABLE", reason = "fewer_than_2_patients_in_group"
          )
          next
        }
        cfit <- fit_independent(counts, cdat, target$group1, target$group2, include_cohort = FALSE)
        cr <- cfit$result[gene %in% selected, .(gene, logFC, standard_error, p_value, fdr)]
        cr[, `:=`(comparison = target$comparison, compartment = target$compartment,
                  cohort_block = cohort_name, n_samples = cfit$n_samples, n_patients = cfit$n_patients)]
        cohort_effects[[length(cohort_effects) + 1L]] <- cr
        cohort_run_rows[[length(cohort_run_rows) + 1L]] <- data.table(
          comparison = target$comparison, compartment = target$compartment,
          cohort_block = cohort_name, status = "PASS", reason = NA_character_,
          n_samples = cfit$n_samples, n_patients = cfit$n_patients
        )
      }
      cohort_dt <- rbindlist(cohort_effects, fill = TRUE)
      if (nrow(cohort_dt) > 0L) {
        fwrite(cohort_dt, file.path(detail_dir, paste0(prefix, "__cohort_effects.tsv.gz")),
               sep = "\t", na = "NA", compress = "gzip")
        meta_dt <- cohort_dt[, {
          meta <- meta_gene(logFC, standard_error)
          list(n_cohorts = meta$k, cohort_sign_fraction = mean(sign(logFC) == sign(cand[gene == .BY$gene, main_logFC][1])),
               beta_fixed = meta$beta_fixed, p_fixed = meta$p_fixed,
               beta_random = meta$beta_random, p_random = meta$p_random, i2 = meta$i2)
        }, by = gene]
        meta_dt[, `:=`(comparison = target$comparison, compartment = target$compartment)]
        cohort_rows[[length(cohort_rows) + 1L]] <- meta_dt
      }

      # LOSO 始终采用独立患者设计，避免把近似相关参数变化混入 study 稳健性。
      loso_effects <- list()
      for (study_name in unique(dat_all$study)) {
        ldat <- select_one_per_patient(dat_all[study != study_name])
        group_n <- ldat[, .(n_patients = uniqueN(patient_key)), by = analysis_group]
        if (nrow(group_n) != 2L || min(group_n$n_patients) < 2L) {
          loso_run_rows[[length(loso_run_rows) + 1L]] <- data.table(
            comparison = target$comparison, compartment = target$compartment,
            omitted_study = study_name, status = "NOT_EVALUABLE", reason = "fewer_than_2_patients_in_group"
          )
          next
        }
        lfit <- tryCatch(fit_independent(counts, ldat, target$group1, target$group2, include_cohort = TRUE),
                         error = function(e) e)
        if (inherits(lfit, "error")) {
          loso_run_rows[[length(loso_run_rows) + 1L]] <- data.table(
            comparison = target$comparison, compartment = target$compartment,
            omitted_study = study_name, status = "NOT_EVALUABLE", reason = conditionMessage(lfit)
          )
          next
        }
        lr <- lfit$result[gene %in% selected, .(gene, logFC, p_value, fdr)]
        lr[, `:=`(comparison = target$comparison, compartment = target$compartment,
                  omitted_study = study_name, n_samples = lfit$n_samples, n_patients = lfit$n_patients)]
        loso_effects[[length(loso_effects) + 1L]] <- lr
        loso_run_rows[[length(loso_run_rows) + 1L]] <- data.table(
          comparison = target$comparison, compartment = target$compartment,
          omitted_study = study_name, status = "PASS", reason = NA_character_,
          n_samples = lfit$n_samples, n_patients = lfit$n_patients
        )
      }
      loso_dt <- rbindlist(loso_effects, fill = TRUE)
      if (nrow(loso_dt) > 0L) {
        fwrite(loso_dt, file.path(detail_dir, paste0(prefix, "__loso_effects.tsv.gz")),
               sep = "\t", na = "NA", compress = "gzip")
        loso_sum <- merge(loso_dt, cand[, .(gene, main_logFC)], by = "gene", all.x = TRUE)[,
          .(n_loso = .N, loso_sign_fraction = mean(sign(logFC) == sign(main_logFC)),
            loso_median_logFC = median(logFC), loso_min_abs_logFC = min(abs(logFC))), by = gene]
        loso_sum[, `:=`(comparison = target$comparison, compartment = target$compartment)]
        loso_rows[[length(loso_rows) + 1L]] <- loso_sum
      }

      target_rows[[length(target_rows) + 1L]] <- data.table(
        comparison = target$comparison, compartment = target$compartment, status = "PASS",
        n_main_fdr_010_nondefinition = nrow(sig), n_candidates = length(selected),
        n_evaluable_cohorts = uniqueN(cohort_dt$cohort_block), n_loso_runs = uniqueN(loso_dt$omitted_study)
      )
      rm(counts, main_res, sens_dt, cohort_dt, loso_dt)
      gc(verbose = FALSE)
    }, error = function(e) {
      log_event("ERROR", "stability target failed", comparison = target$comparison,
                compartment = target$compartment, error = conditionMessage(e))
      target_rows[[length(target_rows) + 1L]] <<- data.table(
        comparison = target$comparison, compartment = target$compartment,
        status = "FAIL", error = conditionMessage(e))
    })
  }

  candidates <- rbindlist(candidate_rows, fill = TRUE)
  sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)
  cohort_meta <- rbindlist(cohort_rows, fill = TRUE)
  loso_summary <- rbindlist(loso_rows, fill = TRUE)
  target_summary <- rbindlist(target_rows, fill = TRUE)
  cohort_runs <- rbindlist(cohort_run_rows, fill = TRUE)
  loso_runs <- rbindlist(loso_run_rows, fill = TRUE)

  fwrite(candidates, file.path(outdir, "stability_candidate_genes.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
  fwrite(sensitivity, file.path(outdir, "sensitivity_gene_concordance.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
  fwrite(cohort_meta, file.path(outdir, "cohort_meta_summary.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
  fwrite(loso_summary, file.path(outdir, "loso_gene_summary.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
  fwrite(cohort_runs, file.path(outdir, "cohort_run_status.tsv"), sep = "\t", na = "NA")
  fwrite(loso_runs, file.path(outdir, "loso_run_status.tsv"), sep = "\t", na = "NA")
  fwrite(target_summary, file.path(outdir, "target_stability_summary.tsv"), sep = "\t", na = "NA")

  n_pass <- sum(target_summary$status == "PASS")
  out_gate <- data.table(
    gate = "G2C_stability_execution", status = if (n_pass == 12L) "PASS" else "FAIL",
    n_targets_pass = n_pass, n_targets_total = 12L,
    note = "Execution gate only; biological stability is gene- and target-specific and legacy survival audit remains required."
  )
  fwrite(out_gate, file.path(outdir, "g2c_stability_gate.tsv"), sep = "\t", na = "NA")
  fwrite(data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version",
             "edgeR_version", "limma_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string,
              as.character(packageVersion("edgeR")), as.character(packageVersion("limma")), outdir)
  ), file.path(outdir, "02d_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pseudobulk stability execution completed",
            runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
            status = out_gate$status[[1]], output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
