#!/usr/bin/env Rscript
# Purpose: Select the nine-gene lactate-associated signature using the adopted multi-cancer filters.
# Task: acid_signature_selection; step: s01
# Inputs (DCF fields): candidate_genes,deg_results,tumor_scores,tpm_root
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task acid_signature_selection --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 01_select_signature_all_sample.R --config PATH [--check-config]\n",
      "       Rscript 01_select_signature_all_sample.R --help\n\n",
      "DCF fields: candidate_genes, deg_results, tumor_scores, tpm_root, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation recalculates all-sample selection stages and exports a compatible signature object.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (!(length(args) %in% c(2L, 3L)) || args[1L] != "--config" || !nzchar(args[2L]) ||
    (length(args) == 3L && args[3L] != "--check-config")) {
  usage()
  stop("Expected --config PATH, optionally followed by --check-config.", call. = FALSE)
}
flags <- if (length(args) == 3L) args[3L] else character()

# Resolve DCF-relative paths through existing ancestors and symlinks.
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
config_path <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path, mode = 4L) != 0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size <= 0) {
  stop("Configuration must be a readable nonempty regular file.", call. = FALSE)
}
cfg <- read.dcf(config_path)
required_fields <- c("candidate_genes", "deg_results", "tumor_scores", "tpm_root", "output_dir")
if (nrow(cfg) != 1L || !setequal(colnames(cfg), required_fields) ||
    ncol(cfg) != length(required_fields)) {
  stop("Configuration must contain one DCF record with exactly: ",
       paste(required_fields, collapse = ", "), call. = FALSE)
}
cfg <- setNames(as.character(cfg[1L, required_fields]), required_fields)
if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME", cfg, fixed = TRUE))) {
  stop("Replace all example placeholders with explicit paths.", call. = FALSE)
}
paths <- vapply(cfg, resolve_path, character(1L), base_dir = dirname(config_path))
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
if ("--check-config" %in% flags) {
  cat("Configuration paths passed metadata checks. No inputs loaded; no analysis or outputs created.\n")
  quit(status = 0L)
}

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

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

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
