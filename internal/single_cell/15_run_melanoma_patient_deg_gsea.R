#!/usr/bin/env Rscript
# Purpose: Fit melanoma patient-level pseudobulk comparisons and programme enrichment.
# Task: melanoma_pseudobulk; step: s15
# Inputs (DCF fields): sample_sheet,counts_manifest,program_catalog,program_membership,fig6_pathways,hallmark_gmt,reactome_gmt,gobp_gmt,counts_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task melanoma_pseudobulk --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 15_run_melanoma_patient_deg_gsea.R --config PATH [--check-config]\n       Rscript 15_run_melanoma_patient_deg_gsea.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("sample_sheet", "counts_manifest", "program_catalog", "program_membership", "fig6_pathways", "hallmark_gmt", "reactome_gmt", "gobp_gmt", "counts_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet", "counts_manifest", "program_catalog", "program_membership", "fig6_pathways", "hallmark_gmt", "reactome_gmt", "gobp_gmt")]
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

# 普通执行定位所选counts；NA和普通缺失仍由源missing-count stop处理。
relocate_input <- function(recorded, input_dir, allow_missing = FALSE) {
  if(length(recorded)!=1L) stop("Invalid recorded input path.")
  if(is.na(recorded)) {
    if(isTRUE(allow_missing)) return(NA_character_)
    stop("Invalid recorded input path.")
  }
  if(!nzchar(recorded)) stop("Invalid recorded input path.")
  leaf <- basename(recorded)
  if(leaf %in% c(".","..") || !nzchar(leaf)) stop("Invalid input basename.")
  candidate <- resolve_path(file.path(input_dir,leaf), dirname(config_path))
  if(!is_within(candidate,input_dir) || identical(candidate,input_dir)) stop("Input escapes configured directory.")
  if(isTRUE(allow_missing) && !file.exists(candidate) && !dir.exists(candidate)) return(candidate)
  if(!file_test("-f",candidate) || file.access(candidate,4L)!=0L || is.na(file.info(candidate)$size) || file.info(candidate)$size<=0) stop("Selected input must be a readable nonempty regular file.")
  candidate
}

cfg_paths <- paths

# ---- Retained source body; only documented path/read adapters ----
cfg_paths <- paths


Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
options(stringsAsFactors = FALSE, width = 220)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
})

SCRIPT_NAME <- "7.3scRNA_c1c4_state_closure/scripts/08_run_fig7_melanoma_patient_level_deg_gsea.R"
CONTRACT_VERSION <- "FIG7_MELANOMA_PATIENT_LEVEL_DEG_GSEA_V0.1_20260812"
CELL_TYPES <- c("CD8 T", "CD4 T", "Monocyte/Macrophage", "cDC", "Fibroblast", "Tumor")

log_event <- function(level, msg, fields = list()) {
  clean <- function(x) gsub('"', "'", as.character(x), fixed = TRUE)
  extra <- if (length(fields) == 0L) "" else paste0(
    " ", paste(sprintf('%s="%s"', names(fields), vapply(fields, clean, character(1L))), collapse = " ")
  )
  cat(sprintf("[%s] %s script=\"%s\" msg=\"%s\"%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, SCRIPT_NAME, clean(msg), extra))
  flush.console()
}

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = TRUE, na = "NA")

read_gmt <- function(path, collection) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  out <- rbindlist(lapply(fields, function(x) {
    if (length(x) < 3L) return(NULL)
    data.table(pathway = x[[1L]], description = x[[2L]], gene = unique(x[-c(1L, 2L)]))
  }), fill = TRUE)
  out[, collection := collection]
  unique(out[nzchar(pathway) & nzchar(gene), .(collection, pathway, description, gene)])
}

select_one_per_patient <- function(dat) {
  dd <- copy(dat)
  dd[, pre_priority := fifelse(pre_post == "Pre", 0L, 1L)]
  setorder(dd, patient, pre_priority, -n_cells, -library_size, sample_ID_long)
  dd <- dd[, .SD[1L], by = patient]
  dd[, pre_priority := NULL]
  dd
}

make_design <- function(dat, ref_state, allow_pre_post = TRUE) {
  dd <- copy(dat)
  dd[, group := factor(sample_C_subtype, levels = c(ref_state, "C1"))]
  use_pre_post <- allow_pre_post && uniqueN(dd$pre_post) > 1L
  if (use_pre_post) dd[, pre_post := factor(pre_post, levels = c("Pre", "Post"))]
  form <- if (use_pre_post) ~ 0 + group + pre_post else ~ 0 + group
  design <- model.matrix(form, dd)
  colnames(design) <- make.names(colnames(design))
  if (qr(design)$rank < ncol(design) && use_pre_post) {
    use_pre_post <- FALSE
    form <- ~ 0 + group
    design <- model.matrix(form, dd)
    colnames(design) <- make.names(colnames(design))
  }
  if (qr(design)$rank < ncol(design)) stop("rank_deficient_group_design")
  list(data = dd, design = design, formula = paste(deparse(form), collapse = ""), uses_pre_post = use_pre_post)
}

prepare_dge <- function(counts, dat, design) {
  idx <- match(dat$unit_id, colnames(counts))
  if (anyNA(idx)) stop("count/sample mismatch")
  y <- DGEList(counts = counts[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = design, min.count = 5)
  if (sum(keep) < 1000L) stop("fewer_than_1000_expressed_genes")
  y <- y[keep, , keep.lib.sizes = FALSE]
  calcNormFactors(y)
}

fit_one <- function(counts, dat, ref_state) {
  payload <- make_design(dat, ref_state, TRUE)
  y <- prepare_dge(counts, payload$data, payload$design)
  block <- factor(payload$data$patient)
  has_repeats <- any(table(block) > 1L)
  use_block <- FALSE
  correlation <- NA_real_
  block_reason <- "no_repeated_patient"

  v0 <- voom(y, payload$design, plot = FALSE)
  if (has_repeats) {
    dc <- tryCatch(duplicateCorrelation(v0, payload$design, block = block), error = function(e) e)
    if (!inherits(dc, "error") && is.finite(dc$consensus.correlation)) {
      use_block <- TRUE
      correlation <- dc$consensus.correlation
      block_reason <- "patient_block_estimated"
      v <- voom(y, payload$design, block = block, correlation = correlation, plot = FALSE)
      fit <- lmFit(v, payload$design, block = block, correlation = correlation)
    } else {
      # 相关性无法估计时按合同退回每患者一个样本，并重新完成过滤和voom。
      independent <- select_one_per_patient(payload$data)
      payload <- make_design(independent, ref_state, TRUE)
      y <- prepare_dge(counts, payload$data, payload$design)
      v <- voom(y, payload$design, plot = FALSE)
      fit <- lmFit(v, payload$design)
      block_reason <- "independent_patient_fallback"
    }
  } else {
    v <- v0
    fit <- lmFit(v, payload$design)
  }

  contrast <- makeContrasts(
    contrasts = paste0("groupC1-group", make.names(ref_state)),
    levels = payload$design
  )
  fit2 <- eBayes(contrasts.fit(fit, contrast), robust = TRUE)
  tt <- as.data.table(topTable(fit2, number = Inf, sort.by = "none"), keep.rownames = "gene")
  setnames(tt, c("P.Value", "adj.P.Val", "logFC"), c("p_value", "fdr", "log2fc_c1_minus_ref"))
  tt[, moderated_t := t]

  summary <- data.table(
    n_samples = nrow(payload$data), n_patients = uniqueN(payload$data$patient),
    n_c1_samples = sum(payload$data$sample_C_subtype == "C1"),
    n_ref_samples = sum(payload$data$sample_C_subtype == ref_state),
    n_c1_patients = uniqueN(payload$data[sample_C_subtype == "C1", patient]),
    n_ref_patients = uniqueN(payload$data[sample_C_subtype == ref_state, patient]),
    n_repeated_patients = payload$data[, .N, by = patient][N > 1L, .N],
    n_patients_spanning_groups = payload$data[, .(n_groups = uniqueN(sample_C_subtype)), by = patient][n_groups > 1L, .N],
    n_genes_tested = nrow(tt), formula = payload$formula,
    pre_post_adjusted = payload$uses_pre_post, patient_block_used = use_block,
    consensus_correlation = correlation, block_reason = block_reason,
    design_columns = ncol(payload$design), design_rank = qr(payload$design)$rank
  )
  list(de = tt, summary = summary)
}

run_fgsea <- function(de, pathway_membership, cell_compartment) {
  pathway_membership <- pathway_membership[
    is.na(pathway_compartment) | pathway_compartment == cell_compartment
  ]
  ranks <- de$moderated_t
  names(ranks) <- de$gene
  ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  sets <- split(pathway_membership$gene, pathway_membership$pathway_key)
  sets <- lapply(sets, unique)
  # 固定单线程与随机种子，保证multilevel P值在正式重跑中可复现。
  set.seed(20260812L)
  fg <- fgseaMultilevel(
    pathways = sets, stats = ranks, minSize = 3L, maxSize = 500L,
    eps = 0, nproc = 1L
  )
  as.data.table(fg)[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1L))]
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  deg_dir <- file.path(outdir, "gene_level_deg")
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- as.list(cfg_paths[c("sample_sheet", "counts_manifest", "program_catalog", "program_membership", "fig6_pathways", "hallmark_gmt", "reactome_gmt", "gobp_gmt")])

  missing <- names(paths)[!file.exists(unlist(paths))]
  if (length(missing)) stop("missing inputs: ", paste(missing, collapse = ", "))

  log_event("START", "Fig7 melanoma patient-level pseudobulk DEG/GSEA", list(contract = CONTRACT_VERSION))
  sample_sheet <- fread(paths$sample_sheet)
  counts_manifest <- fread(paths$counts_manifest)
  catalog <- fread(paths$program_catalog)
  project_membership <- fread(paths$program_membership)
  fig6_selection <- fread(paths$fig6_pathways)

  hallmark <- read_gmt(paths$hallmark_gmt, "MSigDB_Hallmark_2025.1")
  reactome <- read_gmt(paths$reactome_gmt, "MSigDB_Reactome_2025.1")
  gobp <- read_gmt(paths$gobp_gmt, "MSigDB_GOBP_2024.1")
  canonical <- rbindlist(list(hallmark, reactome, gobp), use.names = TRUE)

  # 固定库：全部Hallmark + Fig6冻结通路 + 7.3 canonical program通路。
  fig6_ids <- unique(fig6_selection$pathway)
  stage1_ids <- unique(catalog[source_type == "GMT", source_id])
  canonical_keep <- canonical[
    collection == "MSigDB_Hallmark_2025.1" | pathway %in% union(fig6_ids, stage1_ids)
  ]
  canonical_keep[, `:=`(
    pathway_key = paste0("CANONICAL::", pathway),
    pathway_source = fifelse(pathway %in% fig6_ids, "FIG6_FROZEN_PATHWAY",
                      fifelse(collection == "MSigDB_Hallmark_2025.1", "HALLMARK_ALL", "STAGE1_CANONICAL_PROGRAM")),
    program_id = NA_character_,
    pathway_compartment = NA_character_
  )]

  project_sets <- merge(
    project_membership[, .(program_id, gene)],
    catalog[, .(program_id, compartment, biological_family, source_type, source_id)],
    by = "program_id", all.x = TRUE
  )
  project_sets[, `:=`(
    collection = "LM_pan_cancer_nondefinition_program",
    pathway = program_id,
    description = biological_family,
    pathway_key = paste0("PROJECT::", program_id),
    pathway_source = "STAGE1_NONDEFINITION_PROGRAM",
    pathway_compartment = compartment
  )]
  membership <- unique(rbindlist(list(
    canonical_keep[, .(pathway_key, pathway, collection, description, gene, pathway_source, program_id, pathway_compartment)],
    project_sets[, .(pathway_key, pathway, collection, description, gene, pathway_source, program_id, pathway_compartment)]
  ), use.names = TRUE, fill = TRUE))
  membership[, n_genes_source := uniqueN(gene), by = pathway_key]
  membership <- membership[n_genes_source >= 3L]
  write_tsv(membership, file.path(outdir, "fixed_gsea_pathway_membership.tsv.gz"))

  coverage <- data.table(pathway = fig6_ids)[, present_in_fixed_library := pathway %in% membership$pathway]
  coverage <- merge(coverage, unique(fig6_selection[, .(pathway, pathway_display)]), by = "pathway", all.x = TRUE)
  if (!all(coverage$present_in_fixed_library)) stop("Fig6 pathway coverage gate failed")
  write_tsv(coverage, file.path(outdir, "fig6_pathway_coverage.tsv"))

  comparisons <- data.table(
    comparison = c("C1_vs_C4", "C1_vs_C3", "C1_vs_C2"),
    ref_state = c("C4", "C3", "C2"),
    comparison_role = c("PRIMARY_ENDPOINT", "SECONDARY_CONTEXT", "DESCRIPTIVE_LOW_SUPPORT")
  )
  deg_all <- list(); fg_all <- list(); design_all <- list()

  for (cell_type_value in CELL_TYPES) {
    count_path <- counts_manifest[cell_type == cell_type_value, file][1L]
    if (!is.na(count_path)) count_path <- relocate_input(count_path, unname(cfg_paths["counts_dir"]), allow_missing = TRUE)
    if (is.na(count_path) || !file.exists(count_path)) stop("missing count bundle: ", cell_type_value)
    bundle <- readRDS(count_path)
    counts <- bundle$counts
    for (i in seq_len(nrow(comparisons))) {
      cmp <- comparisons[i]
      dat <- sample_sheet[cell_type == cell_type_value & eligible_ge20 == TRUE & sample_C_subtype %in% c("C1", cmp$ref_state)]
      support <- dat[, .(n_patients = uniqueN(patient)), by = sample_C_subtype]
      evaluable <- nrow(support) == 2L && min(support$n_patients) >= 3L
      log_event("CHECKPOINT", "cell type comparison", list(cell_type = cell_type_value, comparison = cmp$comparison, evaluable = evaluable))
      if (!evaluable) {
        design_all[[length(design_all) + 1L]] <- data.table(
          cell_type = cell_type_value, comparison = cmp$comparison,
          comparison_role = cmp$comparison_role, status = "NOT_EVALUABLE",
          reason = "fewer_than_3_patients_in_group"
        )
        next
      }
      fit <- fit_one(counts, dat, cmp$ref_state)
      de <- fit$de
      de[, `:=`(cell_type = cell_type_value, comparison = cmp$comparison,
                comparison_role = cmp$comparison_role, reference_state = cmp$ref_state)]
      out_stub <- gsub("[^A-Za-z0-9]+", "_", tolower(paste(cell_type_value, cmp$comparison, sep = "__")))
      write_tsv(de, file.path(deg_dir, paste0(out_stub, "__gene_level_deg.tsv.gz")))
      deg_all[[length(deg_all) + 1L]] <- de

      cell_compartment <- unique(dat$compartment)
      if (length(cell_compartment) != 1L) stop("cell type has non-unique compartment")
      fg <- run_fgsea(de, membership, cell_compartment)
      setnames(fg, "pathway", "pathway_key")
      meta <- unique(membership[, .(pathway_key, pathway, collection, description, pathway_source, program_id, pathway_compartment, n_genes_source)])
      fg <- merge(fg, meta, by = "pathway_key", all.x = TRUE)
      fg[, `:=`(cell_type = cell_type_value, comparison = cmp$comparison,
                comparison_role = cmp$comparison_role, reference_state = cmp$ref_state,
                direction = fifelse(NES > 0, "C1_higher", fifelse(NES < 0, paste0(cmp$ref_state, "_higher"), "zero")))]
      fg_all[[length(fg_all) + 1L]] <- fg

      ds <- fit$summary
      ds[, `:=`(cell_type = cell_type_value, comparison = cmp$comparison,
                comparison_role = cmp$comparison_role, reference_state = cmp$ref_state,
                status = "PASS", reason = NA_character_)]
      design_all[[length(design_all) + 1L]] <- ds
    }
  }

  deg_combined <- rbindlist(deg_all, use.names = TRUE, fill = TRUE)
  fg_combined <- rbindlist(fg_all, use.names = TRUE, fill = TRUE)
  design_summary <- rbindlist(design_all, use.names = TRUE, fill = TRUE)
  write_tsv(deg_combined, file.path(outdir, "melanoma_patient_level_gene_deg_all.tsv.gz"))
  write_tsv(fg_combined, file.path(outdir, "melanoma_patient_level_gsea_all.tsv.gz"))
  write_tsv(design_summary, file.path(outdir, "melanoma_patient_level_design_support.tsv"))

  # 绘图handoff完整保留固定库，不按当前结果筛通路；Fig7 preview脚本只可按后续合同确定性取子集。
  handoff <- fg_combined[, .(
    cell_type, comparison, comparison_role, reference_state,
    pathway_key, pathway, description, collection, pathway_source, program_id,
    size, NES, pval, padj, direction, leadingEdge
  )]
  write_tsv(handoff, file.path(outdir, "fig7d_patient_level_gsea_plot_handoff.tsv.gz"))

  checks <- data.table(
    check_id = c("fig6_pathways_complete", "six_cell_types", "three_comparisons", "all_18_targets_accounted",
                 "only_predeclared_low_support_target_not_evaluable",
                 "all_designs_full_rank", "all_primary_prepost_adjusted", "gsea_nonempty",
                 "project_programs_match_cell_compartment", "direction_is_c1_minus_ref"),
    expected = c("TRUE", "6", "3", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"),
    observed = vapply(list(
      all(coverage$present_in_fixed_library), uniqueN(design_summary$cell_type), uniqueN(design_summary$comparison),
      nrow(design_summary) == length(CELL_TYPES) * nrow(comparisons),
      identical(
        design_summary[status == "NOT_EVALUABLE", paste(cell_type, comparison, reason, sep = "||")],
        "Tumor||C1_vs_C2||fewer_than_3_patients_in_group"
      ),
      all(design_summary[status == "PASS", design_rank == design_columns]),
      all(design_summary[comparison == "C1_vs_C4", pre_post_adjusted] %in% TRUE),
      nrow(fg_combined) > 0L,
      nrow(fg_combined[
        pathway_source == "STAGE1_NONDEFINITION_PROGRAM" &
          !((cell_type %in% c("CD8 T", "CD4 T") & pathway_compartment == "T/NK") |
            (cell_type %in% c("Monocyte/Macrophage", "cDC") & pathway_compartment == "Myeloid") |
            (cell_type == "Fibroblast" & pathway_compartment == "Stromal/CAF") |
            (cell_type == "Tumor" & pathway_compartment == "Malignant"))
      ]) == 0L,
      all(deg_combined$comparison %in% comparisons$comparison)
    ), as.character, character(1L))
  )
  checks[, pass := expected == observed]
  write_tsv(checks, file.path(outdir, "fig7_patient_level_deg_gsea_qa.tsv"))
  if (!all(checks$pass)) stop("Fig7 DEG/GSEA QA gate failed")

  manifest <- data.table(
    field = c("contract_version", "analysis_unit", "comparisons", "cell_types", "model",
              "ranking_metric", "gsea_library", "fig6_pathway_coverage", "started", "finished"),
    value = c(CONTRACT_VERSION, "sample x native cell type pseudobulk; patient biological unit",
              paste(comparisons$comparison, collapse = ";"), paste(CELL_TYPES, collapse = ";"),
              "edgeR TMM + limma-voom; pre/post covariate when full rank; duplicateCorrelation patient block",
              "moderated t statistic", "all Hallmark + Fig6 frozen GO/Reactome + Stage1 canonical + nondefinition programs",
              paste0(sum(coverage$present_in_fixed_library), "/", nrow(coverage)),
              format(started, "%Y-%m-%d %H:%M:%S %Z"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  write_tsv(manifest, file.path(outdir, "run_manifest.tsv"))
  capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
  log_event("DONE", "Fig7 melanoma patient-level pseudobulk DEG/GSEA passed", list(outdir = outdir))
}

main()
