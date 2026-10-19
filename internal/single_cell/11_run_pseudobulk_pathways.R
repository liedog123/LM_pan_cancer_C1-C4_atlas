#!/usr/bin/env Rscript
# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Task: single_cell_pseudobulk; step: s11
# Inputs (DCF fields): evidence_gate,integrated_gene_evidence,target_gene_evidence_summary,hallmark_gmt,reactome_gmt,gobp_gmt,c7_gmt,legacy_c7,primary_gene_results_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task single_cell_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 11_run_pseudobulk_pathways.R --config PATH [--check-config]\n       Rscript 11_run_pseudobulk_pathways.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("evidence_gate", "integrated_gene_evidence", "target_gene_evidence_summary", "hallmark_gmt", "reactome_gmt", "gobp_gmt", "c7_gmt", "legacy_c7", "primary_gene_results_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("evidence_gate", "integrated_gene_evidence", "target_gene_evidence_summary", "hallmark_gmt", "reactome_gmt", "gobp_gmt", "c7_gmt", "legacy_c7")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("primary_gene_results_dir")])
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
  library(fgsea)
})

options(stringsAsFactors = FALSE)
set.seed(20260716)
script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02f_run_pseudobulk_pathways_and_finalize_g2.R"

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

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  names(sets) <- vapply(fields, `[`, character(1), 1L)
  sets[lengths(sets) > 0L]
}

run_fgsea <- function(stats, pathways) {
  stats <- stats[is.finite(stats)]
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  res <- fgseaMultilevel(pathways = pathways, stats = stats, minSize = 10L, maxSize = 500L,
                         eps = 0, nproc = 4L)
  res <- as.data.table(res)
  if (nrow(res) == 0L) return(res)
  res[, leadingEdge_genes := vapply(leadingEdge, paste, collapse = ",", character(1))]
  res[, leadingEdge := NULL]
  res
}

run_ora <- function(selected, universe, pathways) {
  selected <- intersect(unique(selected), universe)
  if (length(selected) < 3L) return(data.table())
  rows <- lapply(names(pathways), function(term) {
    members <- intersect(pathways[[term]], universe)
    m <- length(members)
    if (m < 10L || m > 500L) return(NULL)
    overlap <- intersect(selected, members)
    k <- length(overlap)
    if (k == 0L) return(NULL)
    data.table(pathway = term, overlap = k, pathway_size = m, selected_size = length(selected),
               universe_size = length(universe), p_value = phyper(k - 1L, m, length(universe) - m,
                                                                  length(selected), lower.tail = FALSE),
               overlap_genes = paste(overlap, collapse = ","))
  })
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) > 0L) out[, fdr := p.adjust(p_value, method = "BH")]
  out
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  primary_gene_results_dir <- unname(paths["primary_gene_results_dir"])
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (fread(unname(paths["evidence_gate"]))$status[[1]] != "PASS") stop("G2D gate is not PASS")
  fwrite(data.table(gate = "G2E_pathway_and_final_decision", status = "RUNNING",
                    note = "RUNNING prevents stale final G2 status after interruption."),
         file.path(outdir, "g2_final_gate_summary.tsv"), sep = "\t", na = "NA")

  collection_files <- c(
    Hallmark = unname(paths["hallmark_gmt"]),
    Reactome = unname(paths["reactome_gmt"]),
    GO_BP = unname(paths["gobp_gmt"]),
    C7 = unname(paths["c7_gmt"])
  )
  if (!all(file.exists(collection_files))) stop("One or more GMT resources are missing")
  collections <- lapply(collection_files, read_gmt)
  evidence <- fread(unname(paths["integrated_gene_evidence"]))
  target_gene_summary <- fread(unname(paths["target_gene_evidence_summary"]))

  targets <- unique(evidence[, .(comparison, compartment)])
  gsea_path <- file.path(outdir, "gsea_all_results.tsv.gz")
  manifest_path <- file.path(outdir, "gsea_manifest.tsv")
  ora_path <- file.path(outdir, "tier_ab_ora_results.tsv.gz")
  reuse_completed <- all(file.exists(c(gsea_path, manifest_path, ora_path)))
  if (reuse_completed) {
    cached_manifest <- fread(manifest_path)
    reuse_completed <- nrow(cached_manifest) == nrow(targets) * length(collections) * 2L &&
      all(cached_manifest$status == "PASS")
  }

  if (reuse_completed) {
    log_event("CHECKPOINT", "reuse completed pathway results", n_manifest_rows = nrow(cached_manifest))
    gsea <- fread(gsea_path)
    manifest <- cached_manifest
    ora <- fread(ora_path)
  } else {
    gsea_rows <- list()
    gsea_manifest <- list()
    ora_rows <- list()
    log_event("START", "run full-rank and definition-excluded pathway analysis", n_targets = nrow(targets))

    for (i in seq_len(nrow(targets))) {
      target <- targets[i]
      prefix <- paste0(safe_name(target$comparison), "__", safe_name(target$compartment))
      main_file <- relocate_input(paste0(prefix, "__main.tsv.gz"), primary_gene_results_dir)
      main_res <- fread(main_file)
      tier <- evidence[comparison == target$comparison & compartment == target$compartment]
      tier_ab <- tier[evidence_tier %in% c("A_cross_cohort_stable", "B_model_stable_context_dependent")]

      for (collection_name in names(collections)) {
        pathways <- collections[[collection_name]]
        for (rank_scope in c("all_tested", "nondefinition")) {
          rank_df <- if (rank_scope == "all_tested") main_res else main_res[definition_gene == FALSE]
          rank_df <- rank_df[is.finite(t)][order(-abs(t), gene)][!duplicated(gene)]
          stats <- rank_df$t
          names(stats) <- rank_df$gene
          log_event("CHECKPOINT", "fgsea", comparison = target$comparison, compartment = target$compartment,
                    collection = collection_name, rank_scope = rank_scope, n_genes = length(stats))
          fg <- run_fgsea(stats, pathways)
          if (nrow(fg) > 0L) {
            fg[, `:=`(comparison = target$comparison, compartment = target$compartment,
                      collection = collection_name, rank_scope = rank_scope,
                      enriched_side = fifelse(NES > 0, sub("_vs_.*$", "", target$comparison),
                                               sub("^.*_vs_", "", target$comparison)),
                      n_rank_genes = length(stats))]
            gsea_rows[[length(gsea_rows) + 1L]] <- fg
          }
          gsea_manifest[[length(gsea_manifest) + 1L]] <- data.table(
            comparison = target$comparison, compartment = target$compartment,
            collection = collection_name, rank_scope = rank_scope,
            n_rank_genes = length(stats), n_pathways_tested = nrow(fg),
            n_significant_fdr005 = sum(fg$padj < 0.05, na.rm = TRUE), status = "PASS"
          )
        }

        universe <- main_res[definition_gene == FALSE, gene]
        for (side in c("group1_higher", "group2_higher")) {
          selected <- if (side == "group1_higher") tier_ab[main_logFC > 0, gene] else tier_ab[main_logFC < 0, gene]
          ora <- run_ora(selected, universe, pathways)
          if (nrow(ora) > 0L) {
            ora[, `:=`(comparison = target$comparison, compartment = target$compartment,
                       collection = collection_name, side = side)]
            ora_rows[[length(ora_rows) + 1L]] <- ora
          }
        }
      }
    }

    gsea <- rbindlist(gsea_rows, fill = TRUE)
    manifest <- rbindlist(gsea_manifest, fill = TRUE)
    ora <- rbindlist(ora_rows, fill = TRUE)
    setcolorder(gsea, c("comparison", "compartment", "collection", "rank_scope", "enriched_side",
                        "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "n_rank_genes",
                        "leadingEdge_genes"))
    fwrite(gsea, gsea_path, sep = "\t", na = "NA", compress = "gzip")
    fwrite(gsea[padj < 0.05], file.path(outdir, "gsea_significant_results.tsv.gz"),
           sep = "\t", na = "NA", compress = "gzip")
    fwrite(manifest, manifest_path, sep = "\t", na = "NA")
    fwrite(ora, ora_path, sep = "\t", na = "NA", compress = "gzip")
  }

  all_rank <- gsea[rank_scope == "all_tested", .(
    comparison, compartment, collection, pathway, all_NES = NES, all_padj = padj,
    all_leadingEdge = leadingEdge_genes
  )]
  nondef_rank <- gsea[rank_scope == "nondefinition", .(
    comparison, compartment, collection, pathway, nondef_NES = NES, nondef_padj = padj,
    nondef_leadingEdge = leadingEdge_genes
  )]
  pathway_evidence <- merge(all_rank, nondef_rank,
                            by = c("comparison", "compartment", "collection", "pathway"), all = TRUE)
  pathway_evidence[, direction_concordant := sign(all_NES) == sign(nondef_NES)]
  pathway_evidence[, definition_robust := is.finite(nondef_padj) & nondef_padj < 0.05 &
                     is.finite(all_padj) & all_padj < 0.05 & direction_concordant]
  # 仅对双排序均显著且同向的通路计算 leading-edge 重叠；按 target 缓存 Tier A/B 基因。
  pathway_evidence[, leading_tier_ab_n := 0L]
  candidate_idx <- which(pathway_evidence$definition_robust)
  if (length(candidate_idx) > 0L) {
    stable_map <- split(
      evidence[evidence_tier %in% c("A_cross_cohort_stable", "B_model_stable_context_dependent"), gene],
      paste(evidence[evidence_tier %in% c("A_cross_cohort_stable", "B_model_stable_context_dependent"), comparison],
            evidence[evidence_tier %in% c("A_cross_cohort_stable", "B_model_stable_context_dependent"), compartment],
            sep = "||")
    )
    leading_genes <- strsplit(
      fifelse(is.na(pathway_evidence$nondef_leadingEdge[candidate_idx]), "",
              pathway_evidence$nondef_leadingEdge[candidate_idx]),
      ",", fixed = TRUE
    )
    target_keys <- paste(pathway_evidence$comparison[candidate_idx],
                         pathway_evidence$compartment[candidate_idx], sep = "||")
    overlap_n <- mapply(function(gs, key) {
      stable_genes <- stable_map[[key]]
      if (is.null(stable_genes)) stable_genes <- character()
      length(intersect(gs, stable_genes))
    }, leading_genes, target_keys)
    set(pathway_evidence, candidate_idx, "leading_tier_ab_n", as.integer(overlap_n))
  }
  pathway_evidence[, stable_pathway := definition_robust & leading_tier_ab_n >= 3L]
  pathway_evidence[, evidence_class := fcase(
    stable_pathway, "stable_independent",
    is.finite(nondef_padj) & nondef_padj < 0.05 & direction_concordant, "independent_rank_only",
    is.finite(all_padj) & all_padj < 0.05 & (!is.finite(nondef_padj) | nondef_padj >= 0.05), "definition_sensitive",
    is.finite(all_padj) & all_padj < 0.05 & !direction_concordant, "direction_changed",
    default = "not_significant"
  )]
  fwrite(pathway_evidence, file.path(outdir, "integrated_pathway_evidence.tsv.gz"),
         sep = "\t", na = "NA", compress = "gzip")

  curated <- c("Hallmark", "Reactome", "GO_BP")
  pathway_summary <- pathway_evidence[, .(
    n_stable_pathways = sum(stable_pathway, na.rm = TRUE),
    n_stable_curated = sum(stable_pathway & collection %in% curated, na.rm = TRUE),
    n_stable_c7 = sum(stable_pathway & collection == "C7", na.rm = TRUE),
    n_curated_collections_supported = uniqueN(collection[stable_pathway & collection %in% curated]),
    n_definition_sensitive = sum(evidence_class == "definition_sensitive", na.rm = TRUE)
  ), by = .(comparison, compartment)]

  final_decision <- merge(target_gene_summary, pathway_summary,
                          by = c("comparison", "compartment"), all.x = TRUE)
  for (nm in c("n_stable_pathways", "n_stable_curated", "n_stable_c7",
               "n_curated_collections_supported", "n_definition_sensitive")) {
    set(final_decision, which(is.na(final_decision[[nm]])), nm, 0L)
  }
  final_decision[, final_g2_status := fcase(
    target_evidence_status == "CROSS_COHORT_SUPPORTED" & n_stable_curated >= 3L &
      n_curated_collections_supported >= 1L, "SUPPORTED",
    target_evidence_status == "CROSS_COHORT_SUPPORTED", "CONTEXT_DEPENDENT_PATHWAY_LIMITED",
    target_evidence_status == "NO_MAIN_FDR_SIGNAL" & n_stable_curated >= 3L,
      "CONTEXT_DEPENDENT_DISTRIBUTED_PATHWAY_ONLY",
    default = "NO_STABLE_SIGNAL"
  )]
  final_decision[, downstream_role := fcase(
    final_g2_status == "SUPPORTED", "eligible_for_program_association",
    grepl("^CONTEXT_DEPENDENT", final_g2_status), "sensitivity_context_only",
    default = "negative_evidence_only"
  )]
  fwrite(final_decision, file.path(outdir, "final_g2_target_decision.tsv"), sep = "\t", na = "NA")
  fwrite(pathway_summary, file.path(outdir, "target_pathway_summary.tsv"), sep = "\t", na = "NA")

  # 旧 C1/C4 pooled-cell C7 仅与新 C1/C4 直接映射。
  legacy_file <- unname(paths["legacy_c7"])
  legacy <- fread(legacy_file)[comparison == "C1_vs_C4" &
                                project_broad_compartment %in% c("Malignant", "T/NK", "Myeloid", "Stromal/CAF") &
                                padj < 0.05, .(comparison, compartment = project_broad_compartment,
                                               pathway, legacy_NES = NES, legacy_padj = padj)]
  new_c7 <- nondef_rank[comparison == "C1_vs_C4" & collection == "C7",
                         .(comparison, compartment, pathway, new_NES = nondef_NES, new_padj = nondef_padj)]
  legacy_c7 <- merge(legacy, new_c7, by = c("comparison", "compartment", "pathway"), all.x = TRUE)
  legacy_c7[, direction_concordant := sign(legacy_NES) == sign(new_NES)]
  legacy_c7[, survival_class := fcase(
    is.na(new_NES), "lost_not_tested",
    direction_concordant & new_padj < 0.05, "retained",
    direction_concordant, "attenuated",
    !direction_concordant & new_padj < 0.05, "reversed",
    default = "lost"
  )]
  fwrite(legacy_c7, file.path(outdir, "legacy_c1c4_c7_survival.tsv.gz"),
         sep = "\t", na = "NA", compress = "gzip")
  legacy_c7_summary <- legacy_c7[, .(
    n_legacy_significant = .N, n_retained = sum(survival_class == "retained"),
    n_attenuated = sum(survival_class == "attenuated"), n_reversed = sum(survival_class == "reversed"),
    n_lost = sum(survival_class %in% c("lost", "lost_not_tested")),
    retained_fraction = mean(survival_class == "retained"),
    direction_concordant_fraction = mean(direction_concordant, na.rm = TRUE)
  ), by = .(comparison, compartment)]
  fwrite(legacy_c7_summary, file.path(outdir, "legacy_c1c4_c7_survival_summary.tsv"), sep = "\t", na = "NA")

  n_assessed <- nrow(final_decision)
  n_supported <- sum(final_decision$final_g2_status == "SUPPORTED")
  n_context <- sum(grepl("^CONTEXT_DEPENDENT", final_decision$final_g2_status))
  n_negative <- sum(final_decision$final_g2_status == "NO_STABLE_SIGNAL")
  c1c4 <- final_decision[comparison == "C1_vs_C4"]
  gate_status <- if (n_assessed == 12L && all(c1c4$final_g2_status %in% c("SUPPORTED", "CONTEXT_DEPENDENT_PATHWAY_LIMITED"))) {
    "PASS_TARGET_SPECIFIC"
  } else {
    "FAIL"
  }
  final_gate <- data.table(
    gate = "G2_final_pseudobulk", status = gate_status,
    n_targets_assessed = n_assessed, n_supported = n_supported,
    n_context_dependent = n_context, n_no_stable_signal = n_negative,
    note = "G2 is target-specific sample/patient-level association evidence; it does not establish mechanism."
  )
  fwrite(final_gate, file.path(outdir, "g2_final_gate_summary.tsv"), sep = "\t", na = "NA")
  fwrite(data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version",
             "fgsea_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string,
              as.character(packageVersion("fgsea")), outdir)
  ), file.path(outdir, "02f_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pathway analysis and final G2 decision completed",
            runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
            status = gate_status, output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
