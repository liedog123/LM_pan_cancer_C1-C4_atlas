# Purpose: signature acid analysis stages.
# Configuration: config/states.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Select the nine-gene lactate-associated signature using the adopted multi-cancer filters.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/signature_selection/01_select_signature_all_sample.R
stage_signature_selection_01_select_signature_all_sample <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Select the nine-gene lactate-associated signature using the adopted multi-cancer filters.
# Task: acid_signature_selection; step: s01
# Inputs (DCF fields): candidate_genes,deg_results,tumor_scores,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_signature_selection --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
candidate_file <- unname(paths["candidate_genes"])
deg_file <- unname(paths["deg_results"])
score_file <- unname(paths["tumor_scores"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
# Inspect path metadata only; actual project keys are validated during an ordinary run.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
tpm_metadata_files <- list.files(tpm_root, pattern = "_mrna_expr_tpm\\.rdata$",
                                 recursive = TRUE, full.names = TRUE)
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L),
                         winslash = "/", mustWork = TRUE)
required_files <- c(candidate_file, deg_file, score_file, tpm_real_paths)

# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

script_name <- "01_select_signature_all_sample.R"
start_time <- Sys.time()

log_event <- function(level, msg, ...) {
  extra <- list(...)
  fields <- if (length(extra)) {
    paste(vapply(names(extra), function(nm) {
      val <- gsub('"', "'", as.character(extra[[nm]]), fixed = TRUE)
      paste0(nm, '=\"', val, '\"')
    }, character(1)), collapse = " ")
  } else ""
  cat(sprintf(
    "[%s] %s script=%s msg=\"%s\"%s%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
    gsub('"', "'", msg, fixed = TRUE), if (nzchar(fields)) " " else "", fields
  ))
  flush.console()
}

bind_rows_base <- function(x) {
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}
`write_tsv` <- bind_stage_helper(shared_write_tsv_052, environment())

load_object <- function(path, expected = NULL) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (!is.null(expected)) {
    if (!expected %in% loaded) stop("Missing object ", expected, " in ", path)
    return(env[[expected]])
  }
  if (length(loaded) != 1L) stop("Expected one object in ", path, "; found: ", paste(loaded, collapse = ","))
  env[[loaded]]
}

safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3L) return(list(n = n, r = NA_real_, p = NA_real_, reason = "n_complete_lt_3"))
  if (length(unique(x[ok])) < 2L) return(list(n = n, r = NA_real_, p = NA_real_, reason = "zero_variance_gene"))
  if (length(unique(y[ok])) < 2L) return(list(n = n, r = NA_real_, p = NA_real_, reason = "zero_variance_score"))
  fit <- tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = "pearson")), error = identity)
  if (inherits(fit, "error")) return(list(n = n, r = NA_real_, p = NA_real_, reason = paste0("cor_test_error:", conditionMessage(fit))))
  list(n = n, r = unname(fit$estimate), p = fit$p.value, reason = "evaluable")
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3L) return(list(n = n, r = NA_real_, reason = "n_complete_lt_3"))
  if (length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(list(n = n, r = NA_real_, reason = "zero_variance_pair"))
  }
  r <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  if (!is.finite(r)) return(list(n = n, r = NA_real_, reason = "nonfinite_correlation"))
  list(n = n, r = unname(r), reason = "evaluable")
}

# Format-only bridge: reuse historical table ordering and Count >= 6 selection.
build_final_selected_genes <- function(stage3, projects) {
  gene_postive_list <- setNames(lapply(projects, function(project) {
    stage3$gene[stage3$project == project & stage3$mode == "all_sample" & stage3$pass]
  }), projects)
all_genes <- unlist(gene_postive_list, use.names = FALSE)
gene_counts <- as.data.frame(table(all_genes))
colnames(gene_counts) <- c("Gene", "Count")

final_selected_genes <- subset(gene_counts, Count >= 6)
  final_selected_genes
}

main <- function() {
  result_dir <- out_dir
  qa_dir <- file.path(out_dir, "qa")
  for (path in c(result_dir, qa_dir)) {
    if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(path)) {
      stop("Cannot create output directory.", call. = FALSE)
    }
  }
  log_event("START", "all-sample selection begin", input = candidate_file, output = result_dir)

  candidate <- read.delim(candidate_file, check.names = FALSE)
  stopifnot(nrow(candidate) == 142L, length(unique(candidate$gene)) == 142L, !anyNA(candidate$gene))
  genes <- candidate$gene
  degs <- load_object(deg_file, "DEGs_list")
  scores <- load_object(score_file, "tumor_gsva_list")

  projects <- sort(intersect(names(degs), names(scores)))
  stopifnot(length(projects) == 21L, setequal(names(degs), names(scores)))

  input_inventory <- list()
  stage0 <- list()
  stage1 <- list()
  stage2 <- list()
  stage3 <- list()
  pairs <- list()

  static_inputs <- c(candidate_file, deg_file, score_file)
  for (f in static_inputs) {
    input_inventory[[length(input_inventory) + 1L]] <- data.frame(
      input_type = "static", project = NA_character_, path = f,
      size_bytes = file.info(f)$size, md5 = unname(tools::md5sum(f))
    )
  }

  for (project in projects) {
    tpm_file <- file.path(tpm_root, project, paste0(project, "_mrna_expr_tpm.rdata"))
    if (!file.exists(tpm_file)) stop("Missing TPM: ", tpm_file)
    tpm <- load_object(tpm_file, "mrna_expr_tpm")
    tpm <- as.matrix(tpm)
    storage.mode(tpm) <- "double"
    if (is.null(rownames(tpm)) || is.null(colnames(tpm))) stop("TPM dimnames missing: ", project)
    if (anyDuplicated(rownames(tpm))) stop("Duplicated TPM gene names: ", project)
    if (anyDuplicated(colnames(tpm))) stop("Duplicated TPM sample IDs: ", project)

    score_mat <- as.matrix(scores[[project]])
    if (nrow(score_mat) != 1L || is.null(colnames(score_mat))) stop("Invalid score matrix: ", project)
    score_vec <- as.numeric(score_mat[1, ])
    names(score_vec) <- colnames(score_mat)
    if (anyDuplicated(names(score_vec))) stop("Duplicated score sample IDs: ", project)

    all_samples <- intersect(colnames(tpm), names(score_vec))
    sample_code <- suppressWarnings(as.integer(substr(all_samples, 14, 15)))
    if (anyNA(sample_code)) stop("Unparseable TCGA sample code: ", project)
    tumor_samples <- all_samples[sample_code < 10L]
    normal_samples <- all_samples[sample_code >= 10L]
    if (!length(all_samples) || !length(tumor_samples)) stop("Empty sample scope: ", project)

    quant <- quantile(as.numeric(tpm), probs = c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
    stage0[[length(stage0) + 1L]] <- data.frame(
      project = project, n_tpm_samples = ncol(tpm), n_score_samples = length(score_vec),
      n_all_intersection = length(all_samples), n_tumor_intersection = length(tumor_samples),
      n_normal_intersection = length(normal_samples), n_tpm_genes = nrow(tpm),
      tpm_min = quant[1], tpm_q25 = quant[2], tpm_median = quant[3],
      tpm_q75 = quant[4], tpm_q99 = quant[5], tpm_max = quant[6]
    )
    input_inventory[[length(input_inventory) + 1L]] <- data.frame(
      input_type = "tpm", project = project, path = tpm_file,
      size_bytes = file.info(tpm_file)$size, md5 = unname(tools::md5sum(tpm_file))
    )

    deg <- as.data.frame(degs[[project]])
    required_deg <- c("row", "baseMean", "log2FoldChange", "padj")
    if (!all(required_deg %in% names(deg))) stop("DEG columns missing: ", project)
    if (anyDuplicated(as.character(deg$row))) stop("Duplicated DEG genes: ", project)
    idx <- match(genes, as.character(deg$row))
    s1 <- data.frame(
      project = project, gene = genes, in_deg = !is.na(idx),
      baseMean = deg$baseMean[idx], log2FoldChange = deg$log2FoldChange[idx], padj = deg$padj[idx]
    )
    s1$evaluable <- s1$in_deg & is.finite(s1$baseMean) & is.finite(s1$log2FoldChange) & is.finite(s1$padj)
    s1$pass <- s1$evaluable & s1$baseMean > 100 & s1$log2FoldChange >= 0.5 & s1$padj < 0.05
    s1$reason <- ifelse(!s1$in_deg, "gene_absent_from_deg",
      ifelse(!s1$evaluable, "nonfinite_deg_stat",
        ifelse(s1$pass, "pass", "threshold_not_met")))
    stage1[[length(stage1) + 1L]] <- s1

    for (mode in c("all_sample")) {
      use_samples <- if (mode == "tumor_only") tumor_samples else all_samples
      s1_genes <- s1$gene[s1$pass]
      cor_rows <- vector("list", length(s1_genes))
      for (i in seq_along(s1_genes)) {
        gene <- s1_genes[i]
        if (!gene %in% rownames(tpm)) {
          cor_rows[[i]] <- data.frame(project = project, mode = mode, gene = gene,
            n_complete = 0L, correlation = NA_real_, pvalue = NA_real_, reason = "gene_absent_from_tpm")
        } else {
          ans <- safe_pearson(as.numeric(tpm[gene, use_samples]), score_vec[use_samples])
          cor_rows[[i]] <- data.frame(project = project, mode = mode, gene = gene,
            n_complete = ans$n, correlation = ans$r, pvalue = ans$p, reason = ans$reason)
        }
      }
      s2 <- bind_rows_base(cor_rows)
      if (nrow(s2)) {
        s2$fdr <- NA_real_
        eval_idx <- which(s2$reason == "evaluable" & is.finite(s2$pvalue))
        s2$fdr[eval_idx] <- p.adjust(s2$pvalue[eval_idx], method = "BH")
        s2$pass <- s2$reason == "evaluable" & s2$correlation > 0.25 & s2$fdr < 0.05
      } else {
        s2 <- data.frame(project = character(), mode = character(), gene = character(),
          n_complete = integer(), correlation = numeric(), pvalue = numeric(), reason = character(),
          fdr = numeric(), pass = logical())
      }
      stage2[[length(stage2) + 1L]] <- s2

      s2_genes <- s2$gene[s2$pass]
      pair_rows <- list()
      if (length(s2_genes) >= 2L) {
        cmb <- combn(s2_genes, 2, simplify = FALSE)
        pair_rows <- lapply(cmb, function(gp) {
          ans <- safe_spearman(as.numeric(tpm[gp[1], use_samples]), as.numeric(tpm[gp[2], use_samples]))
          data.frame(project = project, mode = mode, gene1 = gp[1], gene2 = gp[2],
            n_complete = ans$n, spearman_r = ans$r, reason = ans$reason)
        })
      }
      pair_df <- bind_rows_base(pair_rows)
      if (!nrow(pair_df)) {
        pair_df <- data.frame(project = character(), mode = character(), gene1 = character(),
          gene2 = character(), n_complete = integer(), spearman_r = numeric(), reason = character())
      }
      pairs[[length(pairs) + 1L]] <- pair_df

      synergy_rows <- lapply(s2_genes, function(gene) {
        gp <- pair_df[pair_df$gene1 == gene | pair_df$gene2 == gene, , drop = FALSE]
        valid <- gp$spearman_r[gp$reason == "evaluable" & is.finite(gp$spearman_r)]
        positive <- valid[valid > 0]
        positive_ratio <- if (length(valid)) sum(valid > 0) / length(valid) else NA_real_
        geo <- if (length(positive)) exp(mean(log(positive))) else NA_real_
        arithmetic <- if (length(valid)) mean(valid) else NA_real_
        evaluable <- length(valid) == (length(s2_genes) - 1L) && length(valid) > 0L
        pass <- evaluable && positive_ratio >= 1 && is.finite(geo) && geo > 0.25
        data.frame(project = project, mode = mode, gene = gene,
          n_other_genes = max(length(s2_genes) - 1L, 0L), n_evaluable_pairs = length(valid),
          positive_ratio = positive_ratio, geometric_mean_r = geo,
          arithmetic_mean_r = arithmetic, evaluable = evaluable, pass = pass,
          reason = if (pass) "pass" else if (!evaluable) "incomplete_or_no_evaluable_pairs" else "threshold_not_met")
      })
      s3 <- bind_rows_base(synergy_rows)
      if (!nrow(s3)) {
        s3 <- data.frame(project = character(), mode = character(), gene = character(),
          n_other_genes = integer(), n_evaluable_pairs = integer(), positive_ratio = numeric(),
          geometric_mean_r = numeric(), arithmetic_mean_r = numeric(), evaluable = logical(),
          pass = logical(), reason = character())
      }
      stage3[[length(stage3) + 1L]] <- s3
    }

    log_event("CHECKPOINT", "project reconstructed", project = project,
      n_all = length(all_samples), n_tumor = length(tumor_samples), n_normal = length(normal_samples),
      stage1_pass = sum(s1$pass))
    rm(tpm)
    gc(verbose = FALSE)
  }

  input_inventory <- bind_rows_base(input_inventory)
  stage0 <- bind_rows_base(stage0)
  stage1 <- bind_rows_base(stage1)
  stage2 <- bind_rows_base(stage2)
  stage3 <- bind_rows_base(stage3)
  pairwise <- bind_rows_base(pairs)

  stage1_summary <- do.call(rbind, lapply(projects, function(p) {
    x <- stage1[stage1$project == p, ]
    data.frame(project = p, n_candidates = nrow(x), n_evaluable = sum(x$evaluable), n_pass = sum(x$pass))
  }))
  stage2_summary <- do.call(rbind, lapply(c("all_sample"), function(m) {
    do.call(rbind, lapply(projects, function(p) {
      x <- stage2[stage2$project == p & stage2$mode == m, ]
      data.frame(project = p, mode = m, n_tested = nrow(x), n_evaluable = sum(x$reason == "evaluable"), n_pass = sum(x$pass))
    }))
  }))
  stage3_summary <- do.call(rbind, lapply(c("all_sample"), function(m) {
    do.call(rbind, lapply(projects, function(p) {
      x <- stage3[stage3$project == p & stage3$mode == m, ]
      data.frame(project = p, mode = m, n_tested = nrow(x), n_evaluable = sum(x$evaluable), n_pass = sum(x$pass))
    }))
  }))

  stage4 <- do.call(rbind, lapply(c("all_sample"), function(m) {
    counts <- table(factor(stage3$gene[stage3$mode == m & stage3$pass], levels = genes))
    data.frame(mode = m, gene = genes, count = as.integer(counts), pass = as.integer(counts) >= 6L)
  }))

  retention <- expand.grid(gene = genes, project = projects, mode = c("all_sample"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key3 <- paste(stage3$gene, stage3$project, stage3$mode, sep = "||")
  retention$stage3_pass <- stage3$pass[match(paste(retention$gene, retention$project, retention$mode, sep = "||"), key3)]
  retention$stage3_pass[is.na(retention$stage3_pass)] <- FALSE

  s1_unique <- length(unique(stage1$gene[stage1$pass]))
  flow <- do.call(rbind, lapply(c("all_sample"), function(m) {
    final <- stage4$gene[stage4$mode == m & stage4$pass]
    data.frame(mode = m,
      stage = c("candidate_universe", "tumor_up_any_cancer", "tumor_program_correlated_any_cancer",
        "within_cancer_synergy_any_cancer", "cross_cancer_count_ge_6"),
      value = c(142L, s1_unique, length(unique(stage2$gene[stage2$mode == m & stage2$pass])),
        length(unique(stage3$gene[stage3$mode == m & stage3$pass])), length(final)))
  }))

  output_map <- c(
    input_inventory = file.path(result_dir, "stage0_input_inventory.tsv"),
    stage0 = file.path(result_dir, "stage0_sample_and_scale_inventory.tsv"),
    stage1 = file.path(result_dir, "stage1_tumor_up_status_long.tsv"),
    stage1_summary = file.path(result_dir, "stage1_cancer_summary.tsv"),
    stage2 = file.path(result_dir, "stage2_tumor_program_correlation_long.tsv"),
    stage2_summary = file.path(result_dir, "stage2_cancer_summary.tsv"),
    stage3 = file.path(result_dir, "stage3_within_cancer_synergy_long.tsv"),
    pairwise = file.path(result_dir, "stage3_pairwise_correlations_long.tsv"),
    stage3_summary = file.path(result_dir, "stage3_cancer_summary.tsv"),
    stage4 = file.path(result_dir, "stage4_cross_cancer_counts.tsv"),
    retention = file.path(result_dir, "gene_by_cancer_retention_long.tsv"),
    flow = file.path(result_dir, "reconstructed_flow_counts.tsv")
  )
  objects <- list(input_inventory, stage0, stage1, stage1_summary, stage2, stage2_summary,
    stage3, pairwise, stage3_summary, stage4, retention, flow)
  for (i in seq_along(output_map)) write_tsv(objects[[i]], output_map[[i]])

  final_all <- stage4$gene[stage4$mode == "all_sample" & stage4$pass]
  # Rederive counts from long-format retention without sparse-array dimension changes.
  rederived_all <- vapply(genes, function(g) {
    sum(retention$stage3_pass[retention$mode == "all_sample" & retention$gene == g])
  }, integer(1))
  qa <- data.frame(
    check = c("candidate_142_unique", "project_sets_exact_21", "stage1_complete_grid",
      "all_scope_contains_normal_each_cancer",
      "sample_intersections_nonempty", "stage4_counts_rederive_all",
      "outputs_nonempty"),
    pass = c(
      nrow(candidate) == 142L && length(unique(genes)) == 142L,
      length(projects) == 21L && setequal(names(degs), names(scores)),
      nrow(stage1) == 142L * 21L,
      all(stage0$n_normal_intersection > 0L),
      all(stage0$n_all_intersection > 0L & stage0$n_tumor_intersection > 0L),
      all(stage4$count[stage4$mode == "all_sample"] == rederived_all),
      all(file.exists(output_map) & file.info(output_map)$size > 0)
    ),
    detail = c("142 unique genes", "21 matched cancer objects", paste0("rows=", nrow(stage1)),
      paste0("normal range=", min(stage0$n_normal_intersection), "-", max(stage0$n_normal_intersection)),
      "all intersections >0",
      "counts equal retention matrix", paste0(length(output_map), " outputs"))
  )
  write_tsv(qa, file.path(qa_dir, "reconstruction_machine_checks.tsv"))
  if (!all(qa$pass)) stop("One or more reconstruction QA checks failed")

  final_selected_genes <- build_final_selected_genes(stage3, projects)
  selected_stage4 <- stage4[stage4$mode == "all_sample" & stage4$pass, , drop = FALSE]
  stopifnot(setequal(as.character(final_selected_genes$Gene), selected_stage4$gene),
            identical(as.integer(final_selected_genes$Count),
                      as.integer(selected_stage4$count[match(as.character(final_selected_genes$Gene),
                                                           selected_stage4$gene)])))
  save(final_selected_genes, file = file.path(result_dir, "final_selected_genes.rdata"))
  write_tsv(final_selected_genes, file.path(result_dir, "final_selected_genes.tsv"))
  output_map <- c(output_map,
    final_selected_genes_rdata = file.path(result_dir, "final_selected_genes.rdata"),
    final_selected_genes_tsv = file.path(result_dir, "final_selected_genes.tsv"))
  manifest <- data.frame(
    file = unname(output_map), size_bytes = file.info(output_map)$size,
    md5 = unname(tools::md5sum(output_map)), stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(result_dir, "reconstruction_manifest.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(qa_dir, "reconstruction_sessionInfo.txt"))

  report <- c(
    "# All-sample signature selection",
    paste0("Final selected genes: ", length(final_all)),
    paste(final_all, collapse = ", "),
    "The selected count is calculated, not constrained to nine.",
    "This continuation does not regenerate DESeq2 results or tumour-program scores.",
    "Independent output copies do not replace frozen signatures or Figure source data."
  )
  writeLines(report, file.path(qa_dir, "reconstruction_report.md"))
  log_event("DONE", "all-sample selection complete",
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    all_final_n = length(final_all), output = result_dir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1L)
})
}

# Purpose: Calculate the reported tumour-programme and within-signature gene correlations.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/signature_selection/02_correlate_final_nine_tumor_program.R
stage_signature_selection_02_correlate_final_nine_tumor_program <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Calculate the reported tumour-programme and within-signature gene correlations.
# Task: signature_correlations; step: s02
# Inputs (DCF fields): final_genes,tumor_scores,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task signature_correlations --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
final_file <- unname(paths["final_genes"])
score_file <- unname(paths["tumor_scores"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
# Inspect path metadata only; actual project keys are validated during an ordinary run.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
tpm_metadata_files <- list.files(tpm_root, pattern = "_mrna_expr_tpm\\.rdata$",
                                 recursive = TRUE, full.names = TRUE)
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L),
                         winslash = "/", mustWork = TRUE)
required_files <- c(final_file, score_file, tpm_real_paths)

# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# Preserve the source raw TPM and complete TPM/score sample intersection.
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
`write_tsv` <- bind_stage_helper(shared_write_tsv_048, environment())
`load_object` <- bind_stage_helper(shared_load_object_024, environment())
safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(c(n = sum(ok), r = NA, p = NA))
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  c(n = sum(ok), r = unname(z$estimate), p = z$p.value)
}
formal <- read.delim(final_file, check.names = FALSE)
gene_col <- if ("Gene" %in% names(formal)) "Gene" else names(formal)[1]
genes <- unique(as.character(formal[[gene_col]]))
stopifnot(length(genes) == 9L)

scores <- load_object(score_file, "tumor_gsva_list")
projects <- sort(names(scores))
stopifnot(length(projects) == 21L)

# Validate identifiers without changing the source sample definition.
if (anyNA(genes) || any(!nzchar(genes))) stop("Final gene names must be nonmissing and nonempty")
if (anyNA(projects) || any(!nzchar(projects)) || anyDuplicated(names(scores))) stop("Invalid project names")
cor_rows <- list()
for (project in projects) {
  tpm <- load_object(file.path(tpm_root, project, paste0(project, "_mrna_expr_tpm.rdata")), "mrna_expr_tpm")
  tpm <- as.matrix(tpm)
  storage.mode(tpm) <- "double"
  score_mat <- as.matrix(scores[[project]])
  score_vec <- as.numeric(score_mat[1, ])
  names(score_vec) <- colnames(score_mat)
  samples <- intersect(colnames(tpm), names(score_vec))
  if (is.null(colnames(tpm)) || is.null(rownames(tpm)) || anyDuplicated(colnames(tpm)) ||
      anyDuplicated(rownames(tpm)) || is.null(colnames(score_mat)) || anyDuplicated(colnames(score_mat)) ||
      nrow(score_mat) < 1L) stop("Input matrices require unique sample/gene identifiers in ", project)
  if (!all(genes %in% rownames(tpm))) stop("Final gene absent in ", project)

  project_cor <- do.call(rbind, lapply(genes, function(gene) {
    z <- safe_pearson(as.numeric(tpm[gene, samples]), score_vec[samples])
    data.frame(project = project, gene = gene, n_complete = as.integer(z["n"]),
               pearson_r = z["r"], pvalue = z["p"])
  }))
  project_cor$fdr <- p.adjust(project_cor$pvalue, method = "BH")
  project_cor$significance <- ifelse(project_cor$fdr < 0.001, "***",
    ifelse(project_cor$fdr < 0.01, "**", ifelse(project_cor$fdr < 0.05, "*", "")))
  cor_rows[[project]] <- project_cor

}
cor_long <- do.call(rbind, cor_rows)
cor_long$gene_order <- match(cor_long$gene, genes)
cor_long$cancer_order <- match(cor_long$project, projects)
cor_long <- cor_long[order(cor_long$gene_order, cor_long$cancer_order), ]

write_tsv(cor_long, file.path(out_dir, "panel_a_final_nine_complete_pearson.tsv"))
checks <- data.frame(check = c("genes_9", "cancers_21", "complete_a_189", "no_missing_a"),
  pass = c(length(genes) == 9L, length(projects) == 21L, nrow(cor_long) == 189L, !anyNA(cor_long$pearson_r)))
write_tsv(checks, file.path(out_dir, "machine_checks.tsv"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
if (!all(checks$pass)) stop("Final-nine correlation QA failed")
message("Final-nine correlation table exported.")
}

# Purpose: Calculate the reported tumour-programme and within-signature gene correlations.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/signature_selection/03_correlate_final_nine_gene_pairs.R
stage_signature_selection_03_correlate_final_nine_gene_pairs <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Calculate the reported tumour-programme and within-signature gene correlations.
# Task: signature_correlations; step: s03
# Inputs (DCF fields): final_genes,tumor_scores,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task signature_correlations --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
final_file <- unname(paths["final_genes"])
score_file <- unname(paths["tumor_scores"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
# Inspect path metadata only; actual project keys are validated during an ordinary run.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
tpm_metadata_files <- list.files(tpm_root, pattern = "_mrna_expr_tpm\\.rdata$",
                                 recursive = TRUE, full.names = TRUE)
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L),
                         winslash = "/", mustWork = TRUE)
required_files <- c(final_file, score_file, tpm_real_paths)

# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# Preserve the source raw TPM and complete TPM/score sample intersection.
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
`write_tsv` <- bind_stage_helper(shared_write_tsv_048, environment())
`load_object` <- bind_stage_helper(shared_load_object_024, environment())
safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(c(n = n, rho = NA, p = NA))
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(n = n, rho = unname(z$estimate), p = z$p.value)
}

formal <- read.delim(final_file, check.names = FALSE)
gene_col <- if ("Gene" %in% names(formal)) "Gene" else names(formal)[1]
genes <- unique(as.character(formal[[gene_col]]))
stopifnot(length(genes) == 9L)

scores <- load_object(score_file, "tumor_gsva_list")
projects <- sort(names(scores))
stopifnot(length(projects) == 21L)

# Validate identifiers without changing the source sample definition.
if (anyNA(genes) || any(!nzchar(genes))) stop("Final gene names must be nonmissing and nonempty")
if (anyNA(projects) || any(!nzchar(projects)) || anyDuplicated(names(scores))) stop("Invalid project names")
pairs <- combn(genes, 2, simplify = FALSE)
pair_levels <- vapply(pairs, function(gp) paste(gp, collapse = " / "), character(1))

rows <- list()
for (project in projects) {
  tpm <- load_object(file.path(tpm_root, project, paste0(project, "_mrna_expr_tpm.rdata")), "mrna_expr_tpm")
  tpm <- as.matrix(tpm)
  storage.mode(tpm) <- "double"
  score_mat <- as.matrix(scores[[project]])
  samples <- intersect(colnames(tpm), colnames(score_mat))
  if (is.null(colnames(tpm)) || is.null(rownames(tpm)) || anyDuplicated(colnames(tpm)) ||
      anyDuplicated(rownames(tpm)) || is.null(colnames(score_mat)) || anyDuplicated(colnames(score_mat)))
    stop("Input matrices require unique sample/gene identifiers in ", project)
  if (!all(genes %in% rownames(tpm))) stop("Final gene absent in ", project)
  d <- do.call(rbind, lapply(pairs, function(gp) {
    z <- safe_spearman(as.numeric(tpm[gp[1], samples]), as.numeric(tpm[gp[2], samples]))
    data.frame(project = project, gene1 = gp[1], gene2 = gp[2], gene_pair = paste(gp, collapse = " / "),
               n_complete = as.integer(z["n"]), spearman_rho = z["rho"], pvalue = z["p"])
  }))
  d$fdr <- p.adjust(d$pvalue, method = "BH")
  d$significance <- ifelse(d$fdr < 0.001, "***", ifelse(d$fdr < 0.01, "**", ifelse(d$fdr < 0.05, "*", "")))
  rows[[project]] <- d
}
pair_long <- do.call(rbind, rows)
pair_long$pair_order <- match(pair_long$gene_pair, pair_levels)
pair_long$cancer_order <- match(pair_long$project, projects)
pair_long <- pair_long[order(pair_long$pair_order, pair_long$cancer_order), ]
write_tsv(pair_long, file.path(out_dir, "panel_b_final_nine_pairwise_spearman_complete.tsv"))

checks <- data.frame(check = c("genes_9", "gene_pairs_36", "cancers_21", "complete_pair_grid_756", "no_missing_rho"),
  pass = c(length(genes) == 9L, length(pair_levels) == 36L, length(projects) == 21L,
           nrow(pair_long) == 756L, !anyNA(pair_long$spearman_rho)))
write_tsv(checks, file.path(out_dir, "machine_checks.tsv"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
if (!all(checks$pass)) stop("Final-nine correlation QA failed")
message("Final-nine correlation table exported.")
}

# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/acid_ag_preprocessing/01_score_frozen_acid_signature.R
stage_acid_ag_preprocessing_01_score_frozen_acid_signature <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Task: acid_ag_groups; step: s01
# Inputs (DCF fields): final_genes,clinical_table,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_ag_groups --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
final_genes_path <- unname(paths["final_genes"])
clinical_table_path <- unname(paths["clinical_table"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- c(final_genes_path, clinical_table_path)
# Inspect TPM path metadata only; no serialized objects are loaded.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
expected_projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP",  "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)
tpm_metadata_files <- file.path(tpm_root, expected_projects, paste0(expected_projects, "_mrna_expr_tpm.rdata"))
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM input files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
required_files <- c(required_files, tpm_real_paths)


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_027, environment())
suppressPackageStartupMessages({
  library(GSVA)
  library(data.table)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
final_selected_genes <- load_named_object(final_genes_path, "final_selected_genes")

acid_gene_set <- list( "acid_genes" = as.matrix(final_selected_genes$Gene))
#获取临床数据
all_clinical_data <- fread(clinical_table_path, sep = "\t", header = TRUE)
all_clinical_data <- all_clinical_data[!duplicated(all_clinical_data$sample), ]

#选择癌症类型
projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP",  "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)
`get_tcga_data_tpm` <- bind_stage_helper(shared_get_tcga_data_tpm_029, environment())

#整理临床数据和表达数据
common_sample_select <- function(all_clinical_data){
  tpm_matrix <- get_tcga_data_tpm(project)
  colnames(tpm_matrix) <- gsub("\\.", "-", substr(colnames(tpm_matrix), 1, 15))

  tumor_samples <- which(as.numeric(substr(colnames(tpm_matrix), 14, 15)) < 10)
  tpm_matrix <- tpm_matrix[, tumor_samples]

  common_samples <- intersect(all_clinical_data$sample, colnames(tpm_matrix))
  clinical_data <- all_clinical_data[all_clinical_data$sample %in% common_samples, ]

  return(clinical_data)
}

clinical_data_list <- list()
for (project in projects) {
  clinical_data <- common_sample_select(all_clinical_data)
  clinical_data_list[[project]] <- clinical_data
}
#保留含有表达矩阵的临床样本
save(clinical_data_list, file = file.path(out_dir, "clinical_data_list.rdata"))

#遍历癌症类型计算酸性得分
acid_scores_list <- list()
for (project in projects) {
  # 获取数据
  tpm_matrix <- get_tcga_data_tpm(project)
  colnames(tpm_matrix) <- gsub("\\.", "-", substr(colnames(tpm_matrix), 1, 15))
  all_clinical_data <- clinical_data_list[[project]]

  # 筛选样本
  tumor_samples <- which(as.numeric(substr(colnames(tpm_matrix), 14, 15)) < 10)
  tpm_matrix <- tpm_matrix[, tumor_samples]
  common_samples <- intersect(all_clinical_data$sample, colnames(tpm_matrix))
  tpm_matrix <- tpm_matrix[, common_samples]
  tpm_matrix <- log2(tpm_matrix+1)

  # ssgsea分析：酸性基因集得分
  ssgsea_param <- ssgseaParam(exprData = as.matrix(tpm_matrix),
                              geneSets = acid_gene_set,
                              normalize = TRUE)
  gsva_results <- gsva(ssgsea_param, verbose = FALSE)
  acid_scores_list[[project]] <- gsva_results
}

save(acid_scores_list, file = file.path(out_dir, "per_acid_scores_list.rdata"))
}

# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/acid_ag_preprocessing/02_assign_acid_groups.R
stage_acid_ag_preprocessing_02_assign_acid_groups <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Task: acid_ag_groups; step: s02
# Inputs (DCF fields): acid_scores,final_genes,clinical_data,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_ag_groups --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
acid_scores_path <- unname(paths["acid_scores"])
final_genes_path <- unname(paths["final_genes"])
clinical_data_path <- unname(paths["clinical_data"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- c(acid_scores_path, final_genes_path, clinical_data_path)
# Inspect TPM path metadata only; no serialized objects are loaded.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
expected_projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)
tpm_metadata_files <- file.path(tpm_root, expected_projects, paste0(expected_projects, "_mrna_expr_tpm.rdata"))
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM input files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
required_files <- c(required_files, tpm_real_paths)


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_027, environment())
suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
acid_scores_list <- load_named_object(acid_scores_path, "acid_scores_list")
final_selected_genes <- load_named_object(final_genes_path, "final_selected_genes")
clinical_data_list <- load_named_object(clinical_data_path, "clinical_data_list")
output_root <- out_dir

#选择癌症类型
projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)
gene_list <- as.factor(final_selected_genes$Gene)
`get_tcga_data_tpm` <- bind_stage_helper(shared_get_tcga_data_tpm_029, environment())

# 聚类热图函数
cluster_heatmap <- function(project, gene_list, clinical_data, acid_scores) {
  # 获取数据
  tpm_matrix <- get_tcga_data_tpm(project)
  colnames(tpm_matrix) <- gsub("\\.", "-", substr(colnames(tpm_matrix), 1, 15))

  # 匹配临床数据和表达数据的样本（只保留两者都有数据的样本）
  common_samples <- intersect(clinical_data$sample, colnames(tpm_matrix))
  clinical_data <- clinical_data[clinical_data$sample %in% common_samples, ]
  tpm_matrix <- tpm_matrix[, common_samples]

  # 提取酸性相关的基因并标准化
  tpm_matrix <- log2(tpm_matrix + 1)
  tpm_matrix_selected <- tpm_matrix[rownames(tpm_matrix) %in% gene_list, ]
  acid_expr <- t(scale(t(tpm_matrix_selected)))

  # 对样本进行层次聚类
  distance_matrix <- dist(t(acid_expr), method = "euclidean")
  hclust_result <- hclust(distance_matrix, method = "ward.D")
  clusters <- cutree(hclust_result, k = 3)
  table(clusters)

  # 将 acid_scores 转化为一维命名向量
  acid_scores <- as.numeric(acid_scores[1, ])
  names(acid_scores) <- colnames(acid_scores)

  # 根据酸性评分均值对每个聚类打标签，高中低三类
  grouped_scores <- split(acid_scores, clusters)
  group_means <- sapply(grouped_scores, mean)
  ranked_groups <- rank(-group_means)  # 按均值降序排列
  group_labels <- c("High", "Medium", "Low")[ranked_groups]

  group_mapping <- setNames(group_labels, names(group_means))
  new_clusters <- group_mapping[as.character(clusters)]

  # 将聚类标签添加到临床数据
  clinical_data$Group <- new_clusters[match(clinical_data$sample, colnames(tpm_matrix_selected))]

  # 创建列注释（聚类信息和酸性评分）
  cluster_colors <- list(Cluster = c("Low" = "blue", "Medium" = "gray", "High" = "red"))
  acid_score_colors <- colorRamp2(c(min(acid_scores), mean(acid_scores), max(acid_scores)),
                                  c("blue", "white", "red"))

  cluster_annotation <- HeatmapAnnotation(
    Cluster = factor(new_clusters, levels = c("Low", "Medium", "High")),
    Acid_Score = anno_simple(acid_scores, col = acid_score_colors),
    col = cluster_colors,
    annotation_name_side = "left"
  )

  # 确保保存热图的文件夹存在
  output_dir <- file.path(output_root, "Heatmap")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 保存聚类结果的热图
  png(file = file.path(output_dir, paste0(project, "_clustered_heatmap.png")),
      width = 1000, height = 800, res = 150)
  ht_cluster <- Heatmap(
    acid_expr,
    name = "Z-score Expression",
    row_order = rownames(acid_expr),
    top_annotation = cluster_annotation,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = factor(new_clusters, levels = c("Low", "Medium", "High")),
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    column_title = project
  )
  draw(ht_cluster)
  dev.off()

  # 返回临床数据
  return(clinical_data)
}

# 遍历所有项目并绘制热图，将返回的临床数据整合到列表中
clinical_data_cluster_list <- list()
for (project in projects) {
  clinical_data <- clinical_data_list[[project]]
  acid_scores <- acid_scores_list[[project]]

  clinical_data <- cluster_heatmap(project, gene_list, clinical_data, acid_scores)
  clinical_data_cluster_list[[project]] <- clinical_data
}

#保存
save(clinical_data_cluster_list, file = file.path(out_dir, "clinical_data_cluster_list.rdata"))
}

# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/acid_ag_preprocessing/03_export_processed_input.R
stage_acid_ag_preprocessing_03_export_processed_input <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Score the fixed signature and construct ordered acidification groups and downstream inputs.
# Task: acid_ag_groups; step: s03
# Inputs (DCF fields): acid_scores,clinical_groups,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_ag_groups --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
acid_scores_path <- unname(paths["acid_scores"])
clinical_groups_path <- unname(paths["clinical_groups"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- c(acid_scores_path, clinical_groups_path)
# Inspect TPM path metadata only; no serialized objects are loaded.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
tpm_metadata_files <- list.files(tpm_root, pattern = "_mrna_expr_tpm\\.rdata$", recursive = TRUE, full.names = TRUE)
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM input files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
required_files <- c(required_files, tpm_real_paths)


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_027, environment())
# Adapt only data.table clinical items for the source's data.frame column selection.
as_merge_clinical_tables <- function(clinical_data_cluster_list) {
  lapply(clinical_data_cluster_list, function(x) {
    if (inherits(x, "data.table")) as.data.frame(x) else x
  })
}
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
acid_scores_list <- load_named_object(acid_scores_path, "acid_scores_list")
clinical_data_cluster_list <- load_named_object(clinical_groups_path, "clinical_data_cluster_list")
clinical_data_cluster_list <- as_merge_clinical_tables(clinical_data_cluster_list)

# 获取并预处理 TCGA 表达数据的函数
get_tcga_data_tpm <- function(project) {
  project_dir <- file.path(tpm_root, project)
  if (!dir.exists(project_dir)) {
    stop("项目文件夹不存在：", project_dir)
  }
  tpm_file <- file.path(project_dir, paste0(project, "_mrna_expr_tpm.rdata"))
  if (!file.exists(tpm_file)) {
    stop("TPM 表达数据文件不存在：", tpm_file)
  }
  mrna_expr_tpm <- load_named_object(tpm_file, "mrna_expr_tpm")
  # 标准化样本ID
  colnames(mrna_expr_tpm) <- gsub("\\.", "-", substr(colnames(mrna_expr_tpm), 1, 15))
  # 筛选肿瘤样本：TCGA 样本ID的第14-15位 < 10 代表肿瘤样本
  tumor_samples_indices <- which(as.numeric(substr(colnames(mrna_expr_tpm), 14, 15)) < 10)
  tpm_tumor <- mrna_expr_tpm[, tumor_samples_indices, drop = FALSE]
  return(tpm_tumor)
}

# 筛选样本数 > 30 的癌症类型
eligible_projects <- c()
for (project in names(clinical_data_cluster_list)) {
  clinical_data <- clinical_data_cluster_list[[project]]
  if (!"Group" %in% colnames(clinical_data)) {
    warning(paste("项目", project, "的临床数据中缺少 'Group' 列，已跳过。"))
    next
  }
  group_counts <- table(clinical_data$Group)
  if (all(c("High", "Low") %in% names(group_counts)) &&
      group_counts["High"] > 30 && group_counts["Low"] > 30) {
    eligible_projects <- c(eligible_projects, project)
  }
}

# 需要保留的临床数据列
columns_to_keep <- c(
  "sample", "_PATIENT", "cancer type abbreviation",
  "age_at_initial_pathologic_diagnosis", "gender", "race",
  "ajcc_pathologic_tumor_stage", "clinical_stage",
  "histological_grade", "OS", "OS.time", "DSS", "DSS.time",
  "DFI", "DFI.time", "PFI", "PFI.time", "Group"
)

# 整合数据
pancancer_data_list <- list()
for (project in eligible_projects) {
  message(paste("正在处理项目:", project))

  # 获取当前项目的所有原始数据
  current_clinical <- clinical_data_cluster_list[[project]]
  current_scores <- acid_scores_list[[project]]
  current_expression <- get_tcga_data_tpm(project)

  colnames(current_scores) <- gsub("\\.", "-", substr(colnames(current_scores), 1, 15))

  # 寻找三大数据集共有的样本
  common_samples <- intersect(current_clinical$sample, colnames(current_expression))
  common_samples <- intersect(common_samples, colnames(current_scores))
  common_samples <- sort(common_samples)

  if (length(common_samples) == 0) {
    warning(paste("项目", project, "没有找到任何共同样本，已跳过。"))
    next
  }
  message(paste("...找到", length(common_samples), "个共同样本。"))

  # 根据共同样本过滤并对齐所有数据
  # 过滤临床数据
  clinical_final <- current_clinical[current_clinical$sample %in% common_samples, ]
  clinical_final <- clinical_final[match(common_samples, clinical_final$sample), ]

  # 只保留指定的列
  existing_columns <- intersect(columns_to_keep, colnames(clinical_final))
  clinical_final <- clinical_final[, existing_columns, drop = FALSE]

  # 过滤表达矩阵和得分矩阵
  expression_final <- current_expression[, common_samples, drop = FALSE]
  scores_final <- current_scores[, common_samples, drop = FALSE]

  # 存入最终列表
  project_name_cleaned <- sub("TCGA-", "", project)

  pancancer_data_list[[project_name_cleaned]] <- list(
    clinical_data = clinical_final,
    expression_data = as.matrix(expression_final),
    acid_scores = scores_final
  )
}

save(pancancer_data_list, file = file.path(out_dir, "pancancer_full_data_list.rdata"))
}

# Purpose: Compare Acid_Score with the fixed pH-related transcript panel.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/acid_ag_preprocessing/04_correlate_formal_ph_panel.R
stage_acid_ag_preprocessing_04_correlate_formal_ph_panel <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Compare Acid_Score with the fixed pH-related transcript panel.
# Task: acid_ph_correlations; step: s04
# Inputs (DCF fields): acid_scores,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_ph_correlations --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
acid_scores_path <- unname(paths["acid_scores"])
tpm_root <- unname(paths["tpm_root"])
out_dir <- unname(paths["output_dir"])
required_files <- acid_scores_path
# Inspect TPM path metadata only; no serialized objects are loaded.
if (!dir.exists(tpm_root) || file.access(tpm_root, mode = 5L) != 0L) {
  stop("tpm_root must be a readable searchable directory.", call. = FALSE)
}
expected_projects <- c(
  "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
  "TCGA-COAD",  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP",  "TCGA-LGG",
  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM"
)
tpm_metadata_files <- file.path(tpm_root, expected_projects, paste0(expected_projects, "_mrna_expr_tpm.rdata"))
if (length(tpm_metadata_files) == 0L || any(!file_test("-f", tpm_metadata_files)) ||
    any(file.access(tpm_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(tpm_metadata_files)$size) | file.info(tpm_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular TPM input files.", call. = FALSE)
}
tpm_real_paths <- vapply(tpm_metadata_files, normalizePath, character(1L), winslash = "/", mustWork = TRUE)
required_files <- c(required_files, tpm_real_paths)


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(tpm_root, dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_027, environment())
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
acid_scores_list <- load_named_object(acid_scores_path, "acid_scores_list")
projects <- expected_projects

# 正式成员由作者显式固定；源文件随后覆盖为7成员的语句不纳入公开副本。
phGene <- c("PDK1", "ACSS2", "DHCR7", "ELOVL5", "FDFT1", "HMGCR",
            "HMGCS1", "IDI1", "LDLR", "MSMO1", "MTHFR", "SQLE")
`get_tcga_data_tpm` <- bind_stage_helper(shared_get_tcga_data_tpm_029, environment())

#计算酸性得分和ph基因的相关性
# 初始化一个空矩阵来存储所有癌症类型的相关性结果
correlation_matrix <- matrix(NA, nrow = length(phGene), ncol = length(projects))
rownames(correlation_matrix) <- phGene
colnames(correlation_matrix) <- projects

# 遍历每个癌症类型
for (project in projects) {
  # 提取数据
  expr <- get_tcga_data_tpm(project)
  colnames(expr) <- gsub("\\.", "-", substr(colnames(expr), 1, 15))
  acid_scores <- acid_scores_list[[project]]

  # 确保表达矩阵和酸性评分有共同的样本
  common_samples <- intersect(colnames(expr), colnames(acid_scores))
  expr <- expr[, common_samples, drop = FALSE]
  acid_scores <- acid_scores[, common_samples]

  # 提取 pH 相关基因的表达数据
  ph_matrix <- expr[rownames(expr) %in% phGene, , drop = FALSE]

  # 遍历每个 pH 基因，计算与酸性评分的相关性
  for (i in 1:nrow(ph_matrix)) {
    ph_gene_expression <- as.numeric(ph_matrix[i, ])
    cor_test <- cor.test(acid_scores, ph_gene_expression, method = "spearman")  # 使用 Spearman 相关性
    correlation_matrix[rownames(ph_matrix)[i], project] <- cor_test$estimate    # 存储相关性系数
  }
}
colnames(correlation_matrix) <- sub("^TCGA-", "", colnames(correlation_matrix))


# 仅导出rho矩阵；不保存P值、不新增校正或绘图。
save(correlation_matrix, file = file.path(out_dir, "ph_acid_correlation_matrix.rdata"))
write.table(correlation_matrix, file = file.path(out_dir, "ph_acid_correlation_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
}

# Purpose: Summarise Acid_Score across cancers and tumour-normal comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/acid_ag_preprocessing/05_score_integrated_acid_matrix.R
stage_acid_ag_preprocessing_05_score_integrated_acid_matrix <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Summarise Acid_Score across cancers and tumour-normal comparisons.
# Task: acid_distributions; step: s05
# Inputs (DCF fields): clinical_groups,final_genes,integrated_expression
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_distributions --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
required_files <- unname(paths[setdiff(names(paths), "output_dir")])
out_dir <- unname(paths["output_dir"])


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_027, environment())
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
suppressPackageStartupMessages({library(data.table); library(dplyr); library(GSVA)})
clinical_data_cluster_list <- load_named_object(paths["clinical_groups"], "clinical_data_cluster_list")
final_selected_genes <- load_named_object(paths["final_genes"], "final_selected_genes")
expr_file <- unname(paths["integrated_expression"])
output_rdata <- file.path(out_dir, "pancancer_integrated_ssgsea_acid_scores.rdata")
output_csv <- file.path(out_dir, "pancancer_integrated_ssgsea_acid_scores_by_sample.csv")

all_clinical_data <- bind_rows(clinical_data_cluster_list) %>%
  select(sample, Group, `cancer type abbreviation`) %>%
  rename(Cancer = `cancer type abbreviation`)

message("Loading integrated pan-cancer matrix...")
expr <- fread(expr_file, sep = "\t", header = TRUE, data.table = FALSE)
if (sum(duplicated(expr$sample)) > 0) {
  expr <- expr[!duplicated(expr$sample), ]
}
rownames(expr) <- expr$sample
expr <- expr[, -1, drop = FALSE]

common_samples <- intersect(all_clinical_data$sample, colnames(expr))
all_clinical_data <- all_clinical_data %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, .keep_all = TRUE)

expr_matrix <- as.matrix(expr[, all_clinical_data$sample, drop = FALSE])

gene_list <- as.character(unique(final_selected_genes$Gene))
alias_map <- c(PKM = "PKM2")
mapped_gene_list <- vapply(
  gene_list,
  function(gene) {
    if (gene %in% rownames(expr_matrix)) {
      gene
    } else if (!is.na(alias_map[gene]) && alias_map[gene] %in% rownames(expr_matrix)) {
      unname(alias_map[gene])
    } else {
      NA_character_
    }
  },
  character(1)
)

missing_genes <- gene_list[is.na(mapped_gene_list)]
if (length(missing_genes) > 0) {
  stop("Integrated matrix is missing genes: ", paste(missing_genes, collapse = ", "))
}

if (any(mapped_gene_list != gene_list)) {
  message(
    "Gene alias mapping applied: ",
    paste(gene_list[mapped_gene_list != gene_list], mapped_gene_list[mapped_gene_list != gene_list], sep = " -> ", collapse = "; ")
  )
}

message("Running ssGSEA on integrated matrix...")
ssgsea_param <- ssgseaParam(
  exprData = expr_matrix,
  geneSets = list(Acid_Score = unname(mapped_gene_list)),
  normalize = TRUE
)
ssgsea_scores <- gsva(ssgsea_param, verbose = FALSE)

sample_scores_df <- data.frame(
  sample = colnames(ssgsea_scores),
  Acid_Score = as.numeric(ssgsea_scores[1, ]),
  stringsAsFactors = FALSE
) %>%
  left_join(all_clinical_data, by = "sample")

stopifnot(!any(is.na(sample_scores_df$Acid_Score)))
stopifnot(!any(is.na(sample_scores_df$Cancer)))
stopifnot(!any(is.na(sample_scores_df$Group)))

pancancer_integrated_ssgsea_acid_scores <- list(
  score_matrix = ssgsea_scores,
  sample_scores = sample_scores_df,
  metadata = list(
    method = "ssGSEA",
    normalize = TRUE,
    input_matrix = expr_file,
    input_scale = "EB++AdjustPANCAN integrated matrix used directly",
    gene_count = length(gene_list),
    gene_mapping = stats::setNames(as.list(mapped_gene_list), gene_list),
    sample_count = nrow(sample_scores_df),
    cancer_count = dplyr::n_distinct(sample_scores_df$Cancer),
    mode = "full"
  )
)

save(pancancer_integrated_ssgsea_acid_scores, file = output_rdata)
write.csv(sample_scores_df, output_csv, row.names = FALSE)

summary_df <- sample_scores_df %>%
  group_by(Cancer) %>%
  summarise(
    N = n(),
    Median_Acid = median(Acid_Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Median_Acid))

message("Done.")
message(
  "Range: ",
  sprintf("%.4f", min(sample_scores_df$Acid_Score, na.rm = TRUE)),
  " to ",
  sprintf("%.4f", max(sample_scores_df$Acid_Score, na.rm = TRUE))
)
print(summary_df)

# 沿用panel-source的癌种及High/Low中位摘要；不添加纳入门槛。
top_order_df <- sample_scores_df %>%
  group_by(Cancer) %>%
  summarise(
    Tumor_Median = median(Acid_Score, na.rm = TRUE),
    High_Median = median(Acid_Score[Group == "High"], na.rm = TRUE),
    Low_Median = median(Acid_Score[Group == "Low"], na.rm = TRUE),
    Tumor_N = n(),
    High_N = sum(Group == "High"),
    Low_N = sum(Group == "Low"),
    .groups = "drop"
  ) %>%
  arrange(desc(Tumor_Median)) %>%
  mutate(Order = row_number())

write.csv(summary_df, file.path(out_dir, "pancancer_integrated_ssgsea_acid_summary.csv"), row.names = FALSE)
write.csv(top_order_df, file.path(out_dir, "Fig1C_top_cancer_group_medians.csv"), row.names = FALSE)
}

# Purpose: Summarise Acid_Score across cancers and tumour-normal comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/acid_ag_preprocessing/06_compare_acid_tumor_normal.R
stage_acid_ag_preprocessing_06_compare_acid_tumor_normal <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Summarise Acid_Score across cancers and tumour-normal comparisons.
# Task: acid_distributions; step: s06
# Inputs (DCF fields): validation_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_distributions --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
required_files <- unname(paths[setdiff(names(paths), "output_dir")])
out_dir <- unname(paths["output_dir"])


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_027, environment())
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
full_validation_data <- load_named_object(paths["validation_data"], "full_validation_data")
# 明确截取Fig. 1C Acid_Score；旧五模块及Stage/Grade不属于此入口。
score_cols <- "Acid_Score"

tvn_stats_list <- list()

cancer_list <- unique(full_validation_data$Cancer_Type)

for (cancer in cancer_list) {
  sub_df <- full_validation_data %>% filter(Cancer_Type == cancer)

  # 检查是否同时有 T 和 N
  if (length(unique(sub_df$Sample_Type)) < 2) next

  # 对每个模块 (包括 Acid Score) 计算差异
  for (mod in score_cols) {
    # Wilcoxon Test
    res <- wilcox.test(sub_df[[mod]] ~ sub_df$Sample_Type) # Normal vs Tumor (default factor order)

    # 计算均值差 (Tumor - Normal)
    mean_t <- mean(sub_df[[mod]][sub_df$Sample_Type == "Tumor"], na.rm = TRUE)
    mean_n <- mean(sub_df[[mod]][sub_df$Sample_Type == "Normal"], na.rm = TRUE)
    diff_val <- mean_t - mean_n # 正值表示 Tumor 高

    # 存入列表
    tvn_stats_list[[paste(cancer, mod, sep = "_")]] <- data.frame(
      Cancer = cancer,
      Module = mod,
      Diff_Mean = diff_val,
      P_Value = res$p.value
    )
  }
}
# 合并列表
tvn_stats <- bind_rows(tvn_stats_list)

# 标记显著性星号
tvn_stats$Significance <- case_when(
  tvn_stats$P_Value < 0.001 ~ "***",
  tvn_stats$P_Value < 0.01 ~ "**",
  tvn_stats$P_Value < 0.05 ~ "*",
  TRUE ~ ""
)

write.csv(tvn_stats, file.path(out_dir, "Stats_Tumor_vs_Normal.csv"), row.names = FALSE)

tvn_acid <- tvn_stats %>%
  filter(Module == "Acid_Score") %>%
  select(Cancer, Diff_Mean, P_Value, Significance)

# 中位数按源panel数据摘要另表导出；Diff_Mean仍是均值差。
bottom_summary <- full_validation_data %>%
  select(Cancer = Cancer_Type, Sample_Type, Acid_Score) %>%
  group_by(Cancer, Sample_Type) %>%
  summarise(Median_Acid_Score = median(Acid_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample_Type, values_from = Median_Acid_Score) %>%
  left_join(tvn_acid, by = "Cancer")
write.csv(bottom_summary, file.path(out_dir, "Fig1C_bottom_tumor_normal_medians.csv"), row.names = FALSE)
}
