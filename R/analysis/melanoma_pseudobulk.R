# Purpose: melanoma pseudobulk analysis stages.
# Configuration: config/single_cell.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Fit melanoma patient-level pseudobulk comparisons and programme enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/14_run_melanoma_pseudobulk_programs.R
stage_single_cell_14_run_melanoma_pseudobulk_programs <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("contract", "stage3_gate", "source_object", "frozen_cell_metadata", "sample_labels", "readiness", "catalog", "membership_full", "membership_nondefinition")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# ---- Retained source body; only documented path/read adapters ----
cfg_paths <- paths


Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(edgeR)
  library(GSVA)
})

options(stringsAsFactors = FALSE)

SCRIPT_NAME <- "7.3scRNA_c1c4_state_closure/scripts/04a_build_melanoma_pseudobulk_and_program_scores.R"
CONTRACT_VERSION <- "STAGE4_MELANOMA_C3_VALIDATION_FROZEN_V0.1_20260722"
CELLTYPE_TO_COMPARTMENT <- c(
  "CD8 T" = "T/NK", "CD4 T" = "T/NK",
  "Monocyte/Macrophage" = "Myeloid", "cDC" = "Myeloid",
  "Fibroblast" = "Stromal/CAF", "Tumor" = "Malignant"
)

log_event <- function(level, msg, fields = list()) {
  clean <- function(x) gsub('"', "'", as.character(x), fixed = TRUE)
  extras <- if (length(fields) == 0L) "" else paste0(
    " ", paste(sprintf('%s="%s"', names(fields), vapply(fields, clean, character(1L))), collapse = " ")
  )
  cat(sprintf('[%s] %s script="%s" msg="%s"%s\n',
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, SCRIPT_NAME, clean(msg), extras))
  flush.console()
}
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())
`write_tsv` <- bind_stage_helper(shared_write_tsv_047, environment())
`safe_md5` <- bind_stage_helper(shared_md5_or_na_020, environment())

file_record <- function(role, path, critical = TRUE) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  data.table(
    role = role, path = path, critical = critical, exists = exists,
    bytes = if (exists) as.numeric(info$size) else NA_real_,
    nonempty = exists && !dir.exists(path) && isTRUE(info$size > 0),
    modified_time = if (exists) format(info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
    md5 = safe_md5(path)
  )
}

get_counts <- function(object) {
  mat <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", layer = "counts"),
    error = function(e) tryCatch(
      suppressWarnings(Seurat::GetAssayData(object, assay = "RNA", slot = "counts")),
      error = function(e2) NULL
    )
  )
  if (is.null(mat) || nrow(mat) == 0L || ncol(mat) == 0L) stop("RNA counts matrix is unavailable")
  if (!inherits(mat, "sparseMatrix")) mat <- as(mat, "dgCMatrix")
  mat
}

zscore_rows <- function(mat) {
  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, sd, na.rm = TRUE)
  out <- sweep(mat, 1L, means, "-")
  out <- sweep(out, 1L, sds, "/")
  out[!is.finite(out)] <- NA_real_
  out
}

aggregate_celltype_counts <- function(counts, meta, cell_type) {
  cell_type_value <- as.character(cell_type)[1L]
  md <- meta[get("cell_type") == cell_type_value]
  if (nrow(md) == 0L) stop("No cells for cell type: ", cell_type)
  idx <- match(md$cell_id, colnames(counts))
  if (anyNA(idx)) stop("Frozen cell metadata does not match counts for: ", cell_type)
  sample_levels <- sort(unique(md$sample_ID_long))
  sample_factor <- factor(md$sample_ID_long, levels = sample_levels)
  design <- sparse.model.matrix(~ 0 + sample_factor)
  colnames(design) <- sample_levels
  pb <- counts[, idx, drop = FALSE] %*% design
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sample_levels
  pb
}

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  count_dir <- file.path(outdir, "counts")
  dir.create(count_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- as.list(cfg_paths[c("contract", "stage3_gate", "source_object", "frozen_cell_metadata", "sample_labels", "readiness", "catalog", "membership_full", "membership_nondefinition")])

  log_event("START", "build melanoma sample-by-native-cell-type pseudobulk and frozen program scores", list(output = outdir))
  input_manifest <- rbindlist(Map(function(role, path) file_record(role, path, TRUE), names(paths), unname(unlist(paths))))
  write_tsv(input_manifest, file.path(outdir, "input_file_manifest.tsv"))
  if (nrow(input_manifest[!exists | !nonempty]) > 0L) stop("Stage 4 score inputs are incomplete")

  contract <- fread(paths$contract)
  if (!any(contract$contract_item == "contract_version" & contract$frozen_value == CONTRACT_VERSION)) {
    stop("Stage 4 contract version mismatch")
  }
  stage3 <- fread(paths$stage3_gate)
  if (stage3$status[[1L]] != "PASS" || !isTRUE(stage3$stage4_authorized[[1L]])) {
    stop("Stage 3 did not authorize Stage 4")
  }

  cell_meta <- fread(paths$frozen_cell_metadata)
  cell_meta <- cell_meta[flag_keep_for_scoring %in% TRUE & cell_type %in% names(CELLTYPE_TO_COMPARTMENT)]
  labels <- fread(paths$sample_labels)
  label_cols <- c("sample_ID_long", "patient", "sample_C_subtype", "state", "response_short",
                  "treatment_group", "subtype", "tissue", "day")
  if (!all(label_cols %in% names(labels))) stop("Melanoma label columns are incomplete")
  labels <- unique(labels[, ..label_cols], by = "sample_ID_long")

  log_event("CHECKPOINT", "load frozen melanoma Seurat object", list(path = paths$source_object))
  object <- readRDS(paths$source_object)
  counts <- get_counts(object)
  rm(object)
  gc(verbose = FALSE)
  if (!all(cell_meta$cell_id %in% colnames(counts))) stop("Frozen metadata cells are missing from source object")

  catalog <- fread(paths$catalog)
  membership_full <- fread(paths$membership_full)
  membership_nondefinition <- fread(paths$membership_nondefinition)
  score_rows <- list()
  coverage_rows <- list()
  sample_rows <- list()
  count_manifest_rows <- list()

  for (cell_type in names(CELLTYPE_TO_COMPARTMENT)) {
    cell_type_value <- as.character(cell_type)[1L]
    compartment_value <- unname(CELLTYPE_TO_COMPARTMENT[[cell_type_value]])
    log_event("CHECKPOINT", "aggregate and score native cell type", list(cell_type = cell_type_value, compartment = compartment_value))
    pb <- aggregate_celltype_counts(counts, cell_meta, cell_type_value)
    cell_counts <- cell_meta[get("cell_type") == cell_type_value, .(n_cells = .N), by = sample_ID_long]
    sample_sheet <- merge(data.table(sample_ID_long = colnames(pb)), cell_counts,
                          by = "sample_ID_long", all.x = TRUE)
    sample_sheet <- merge(sample_sheet, labels, by = "sample_ID_long", all.x = TRUE)
    sample_sheet[, `:=`(
      cell_type = cell_type_value,
      compartment = compartment_value,
      unit_id = paste(safe_name(cell_type_value), sample_ID_long, sep = "||"),
      library_size = as.numeric(colSums(pb)),
      pre_post = fifelse(state == "PRE", "Pre", "Post"),
      eligible_ge20 = n_cells >= 20L,
      eligible_ge30 = n_cells >= 30L,
      eligible_ge50 = n_cells >= 50L
    )]
    if (anyNA(sample_sheet$patient) || anyNA(sample_sheet$sample_C_subtype)) {
      stop("Sample labels failed to join for: ", cell_type_value)
    }
    idx_all <- match(sample_sheet$sample_ID_long, colnames(pb))
    pb <- pb[, idx_all, drop = FALSE]
    colnames(pb) <- sample_sheet$unit_id

    count_path <- file.path(count_dir, paste0(safe_name(cell_type_value), "_pseudobulk_counts.rds"))
    saveRDS(list(counts = pb, sample_sheet = sample_sheet), count_path, compress = TRUE)
    count_manifest_rows[[length(count_manifest_rows) + 1L]] <- data.table(
      cell_type = cell_type_value, compartment = compartment_value, file = count_path,
      n_genes = nrow(pb), n_units = ncol(pb), n_units_ge20 = sum(sample_sheet$eligible_ge20)
    )
    sample_rows[[length(sample_rows) + 1L]] <- sample_sheet

    eligible <- sample_sheet[eligible_ge20 == TRUE]
    idx <- match(eligible$unit_id, colnames(pb))
    y <- DGEList(counts = pb[, idx, drop = FALSE])
    y <- calcNormFactors(y)
    logcpm <- cpm(y, log = TRUE, prior.count = 2)

    program_ids <- catalog[get("compartment") == compartment_value, program_id]
    full_sets <- lapply(program_ids, function(id) intersect(membership_full[program_id == id, unique(gene)], rownames(logcpm)))
    nondef_sets <- lapply(program_ids, function(id) intersect(membership_nondefinition[program_id == id, unique(gene)], rownames(logcpm)))
    names(full_sets) <- paste(program_ids, "FULL", sep = "||")
    names(nondef_sets) <- paste(program_ids, "NONDEFINITION", sep = "||")
    gene_sets <- c(full_sets, nondef_sets)
    if (any(lengths(gene_sets) < 3L)) stop("Fewer than 3 scored genes for a frozen program in: ", cell_type_value)

    score_mat <- gsva(
      ssgseaParam(exprData = logcpm, geneSets = gene_sets, minSize = 1L,
                  alpha = 0.25, normalize = TRUE),
      verbose = FALSE
    )
    score_z <- zscore_rows(score_mat)
    raw_long <- melt(as.data.table(t(score_mat), keep.rownames = "unit_id"), id.vars = "unit_id",
                     variable.name = "score_id", value.name = "score_raw")
    z_long <- melt(as.data.table(t(score_z), keep.rownames = "unit_id"), id.vars = "unit_id",
                   variable.name = "score_id", value.name = "score_z_global")
    long <- merge(raw_long, z_long, by = c("unit_id", "score_id"), all = TRUE)
    long[, c("program_id", "score_layer") := tstrsplit(score_id, "||", fixed = TRUE)]
    long[, score_id := NULL]
    long <- merge(long, eligible, by = "unit_id", all.x = TRUE, sort = FALSE)
    long <- merge(long, catalog[, .(program_id, program_order, biological_family, evidence_use,
                                    higher_score_interpretation, allowed_claim)],
                  by = "program_id", all.x = TRUE, sort = FALSE)
    long[, contract_version := CONTRACT_VERSION]
    setorder(long, program_order, score_layer, sample_ID_long)
    score_rows[[length(score_rows) + 1L]] <- long

    coverage_rows[[length(coverage_rows) + 1L]] <- rbindlist(lapply(names(gene_sets), function(id) {
      parts <- strsplit(id, "||", fixed = TRUE)[[1L]]
      data.table(cell_type = cell_type_value, compartment = compartment_value,
                 program_id = parts[[1L]], score_layer = parts[[2L]],
                 n_genes_scored = length(gene_sets[[id]]), n_units_scored = ncol(logcpm))
    }))
    rm(pb, y, logcpm, score_mat, score_z)
    gc(verbose = FALSE)
  }
  rm(counts)
  gc(verbose = FALSE)

  scores <- rbindlist(score_rows, fill = TRUE)
  coverage <- rbindlist(coverage_rows, fill = TRUE)
  sample_sheet_all <- rbindlist(sample_rows, fill = TRUE)
  count_manifest <- rbindlist(count_manifest_rows, fill = TRUE)
  fwrite(scores, file.path(outdir, "melanoma_pseudobulk_program_scores.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  write_tsv(coverage, file.path(outdir, "melanoma_program_score_coverage.tsv"))
  write_tsv(sample_sheet_all, file.path(outdir, "melanoma_pseudobulk_sample_sheet.tsv"))
  write_tsv(count_manifest, file.path(outdir, "melanoma_pseudobulk_counts_manifest.tsv"))

  concordance <- dcast(scores[, .(unit_id, cell_type, program_id, score_layer, score_raw)],
                       unit_id + cell_type + program_id ~ score_layer, value.var = "score_raw")
  concordance <- concordance[, .(
    full_nondefinition_spearman = cor(FULL, NONDEFINITION, method = "spearman", use = "pairwise.complete.obs"),
    max_abs_difference = max(abs(FULL - NONDEFINITION), na.rm = TRUE)
  ), by = .(cell_type, program_id)]
  concordance <- merge(catalog[, .(program_id, program_order, compartment, definition_overlap_fraction)],
                       concordance, by = "program_id", all.y = TRUE)
  setorder(concordance, program_order, cell_type)
  write_tsv(concordance, file.path(outdir, "melanoma_full_nondefinition_score_concordance.tsv"))

  expected_units <- sum(vapply(names(CELLTYPE_TO_COMPARTMENT), function(ct) {
    catalog[compartment == CELLTYPE_TO_COMPARTMENT[[ct]], .N]
  }, integer(1L)))
  checks <- data.table(
    check_id = c(
      "stage3_authorized", "six_predefined_celltypes", "all_source_cells_match",
      "all_sample_labels_joined", "all_expected_celltype_program_units_two_layers",
      "all_scores_finite", "all_count_bundles_nonempty", "qc_program_restricted"
    ),
    passed = c(
      TRUE, uniqueN(sample_sheet_all$cell_type) == 6L, TRUE,
      !anyNA(sample_sheet_all$patient) && !anyNA(sample_sheet_all$sample_C_subtype),
      uniqueN(coverage[, paste(cell_type, program_id)]) == expected_units &&
        all(coverage[, .N, by = .(cell_type, program_id)]$N == 2L),
      all(is.finite(scores$score_raw)) && all(is.finite(scores$score_z_global)),
      nrow(count_manifest) == 6L && all(file.exists(count_manifest$file)) && all(file.info(count_manifest$file)$size > 0),
      nrow(catalog[program_id == "MALIGNANT_CELL_CYCLE_QC" & evidence_use == "QC_ONLY"]) == 1L
    )
  )
  gate_pass <- all(checks$passed)
  gate <- data.table(
    gate = "STAGE4_MELANOMA_PSEUDOBULK_AND_SCORE_FREEZE",
    status = if (gate_pass) "PASS" else "FAIL",
    contract_version = CONTRACT_VERSION,
    n_celltypes = uniqueN(sample_sheet_all$cell_type),
    n_sample_celltype_units = nrow(sample_sheet_all),
    n_units_ge20 = sum(sample_sheet_all$eligible_ge20),
    n_celltype_program_units = expected_units,
    n_score_rows = nrow(scores),
    program_modeling_authorized = gate_pass,
    boundary = "Label-independent melanoma pseudobulk scores only; no C3 comparison has been tested at this step."
  )
  write_tsv(checks, file.path(outdir, "stage4_score_checks.tsv"))
  write_tsv(gate, file.path(outdir, "stage4_score_gate.tsv"))

  output_files <- list.files(outdir, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[file.info(output_files)$isdir == FALSE]
  output_manifest <- rbindlist(lapply(output_files, function(path) file_record(basename(path), path, TRUE)))
  write_tsv(output_manifest, file.path(outdir, "output_manifest.tsv"))
  log_event("DONE", "melanoma pseudobulk and frozen program scores completed", list(
    status = gate$status, n_score_rows = nrow(scores),
    elapsed_minutes = round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2)
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
}

# Purpose: Fit melanoma patient-level pseudobulk comparisons and programme enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/15_run_melanoma_patient_deg_gsea.R
stage_single_cell_15_run_melanoma_patient_deg_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet", "counts_manifest", "program_catalog", "program_membership", "fig6_pathways", "hallmark_gmt", "reactome_gmt", "gobp_gmt")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir")])
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`relocate_input` <- bind_stage_helper(shared_relocate_input_032, environment())

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
`write_tsv` <- bind_stage_helper(shared_write_tsv_047, environment())

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
}
