# Purpose: pan pseudobulk analysis stages.
# Configuration: config/single_cell.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/06_build_pseudobulk_counts.R
stage_single_cell_06_build_pseudobulk_counts <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "cell_metadata", "sample_mapping", "exclusion_genes", "g1_summary")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character()
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


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(SummarizedExperiment)
  library(S4Vectors)
})

options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02a_build_pseudobulk_counts.R"
`log_event` <- bind_stage_helper(shared_log_event_008, environment())
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  counts_dir <- file.path(outdir, "counts")
  dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(seurat_rds = unname(paths["seurat_object"]), cell_metadata = unname(paths["cell_metadata"]),
                sample_mapping = unname(paths["sample_mapping"]), exclusion_genes = unname(paths["exclusion_genes"]),
                g1_summary = unname(paths["g1_summary"]))
  missing <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]
  if (length(missing) > 0L) stop("Missing inputs: ", paste(missing, collapse = ", "))

  log_event("START", "build sample-by-compartment pseudobulk raw counts",
            input = paths$seurat_rds, output = outdir)

  g1 <- fread(paths$g1_summary)
  if (!identical(g1$status[[1]], "PASS")) stop("G1 status is not PASS")
  cell_meta <- fread(paths$cell_metadata)
  mapping <- fread(paths$sample_mapping)
  exclusion <- fread(paths$exclusion_genes)

  log_event("CHECKPOINT", "load Seurat object and extract RNA counts")
  obj <- readRDS(paths$seurat_rds)
  if (!"RNA" %in% names(obj@assays) || !"counts" %in% Layers(obj[["RNA"]])) {
    stop("RNA/counts layer is unavailable")
  }
  counts <- LayerData(obj[["RNA"]], layer = "counts")
  genes <- rownames(counts)
  cells <- colnames(counts)
  if (!inherits(counts, "sparseMatrix")) stop("RNA/counts is not sparse")
  if (anyDuplicated(genes)) stop("Duplicated feature names in RNA/counts")
  ord <- match(cells, cell_meta$cell_id)
  if (anyNA(ord) || length(ord) != ncol(counts)) stop("Cell metadata does not match counts columns")
  cell_meta <- cell_meta[ord]
  if (!identical(cell_meta$cell_id, cells)) stop("Cell order verification failed")
  cell_library <- as.numeric(Matrix::colSums(counts))
  rm(obj)
  invisible(gc())

  non_other <- cell_meta$project_broad_compartment != "Other"
  if (!any(non_other)) stop("No non-Other cells")
  cell_meta[, cell_library_size := cell_library]
  cell_meta[, unit_id := paste(sample_unit_key, project_broad_compartment, sep = "||")]

  units <- cell_meta[non_other, .(
    n_cells = .N,
    expected_library_size = sum(cell_library_size),
    study_from_cells = unique(study)[1],
    cancer_type_from_cells = unique(cancer_type)[1]
  ), by = .(unit_id, sample_unit_key, project_broad_compartment)]
  setorder(units, project_broad_compartment, sample_unit_key)
  units[, unit_index := .I]

  eligible_idx <- which(non_other)
  group_index <- match(cell_meta$unit_id[eligible_idx], units$unit_id)
  if (anyNA(group_index)) stop("Failed to map cells to pseudobulk units")

  log_event("CHECKPOINT", "aggregate sparse counts", n_cells = length(eligible_idx), n_units = nrow(units))
  aggregation_matrix <- sparseMatrix(
    i = eligible_idx,
    j = group_index,
    x = 1,
    dims = c(ncol(counts), nrow(units))
  )
  pb <- counts %*% aggregation_matrix
  rownames(pb) <- genes
  colnames(pb) <- units$unit_id
  rm(aggregation_matrix)
  invisible(gc())

  observed_library <- as.numeric(Matrix::colSums(pb))
  detected_genes <- as.integer(Matrix::colSums(pb > 0))
  units[, `:=`(
    library_size = observed_library,
    detected_genes = detected_genes,
    library_size_match = observed_library == expected_library_size,
    eligible_ge20 = n_cells >= 20L,
    eligible_ge30 = n_cells >= 30L,
    eligible_ge50 = n_cells >= 50L
  )]
  units[, log10_library_size := log10(library_size + 1)]
  units[, log10_n_cells := log10(n_cells + 1)]

  mapping_keep <- setdiff(names(mapping), intersect(names(mapping), c("patient_key")))
  sample_sheet <- merge(units, mapping[, c("sample_unit_key", "patient_key", setdiff(mapping_keep, "sample_unit_key")), with = FALSE],
                        by = "sample_unit_key", all.x = TRUE, sort = FALSE)
  sample_sheet[, cohort_block := paste(study, cancer_type, sep = "||")]
  setorder(sample_sheet, unit_index)
  if (anyNA(sample_sheet$patient_key)) stop("Missing patient mapping after aggregation")
  if (!identical(sample_sheet$unit_id, colnames(pb))) stop("Pseudobulk sample order mismatch")

  gene_meta <- data.table(
    gene = genes,
    is_mitochondrial = grepl("^MT-", genes),
    is_ribosomal = grepl("^RP[SL][0-9]", genes),
    is_definition_exclusion = genes %in% exclusion$gene
  )
  fwrite(gene_meta, file.path(outdir, "pseudobulk_gene_metadata.tsv"), sep = "\t", na = "NA")

  log_event("CHECKPOINT", "save compartment-level SummarizedExperiment objects")
  manifest_rows <- list()
  compartments <- unique(sample_sheet$project_broad_compartment)
  for (compartment in compartments) {
    idx <- which(sample_sheet$project_broad_compartment == compartment)
    comp_counts <- pb[, idx, drop = FALSE]
    comp_meta <- copy(sample_sheet[idx])
    rownames(comp_meta) <- comp_meta$unit_id
    se <- SummarizedExperiment(
      assays = list(counts = comp_counts),
      rowData = DataFrame(gene_meta, row.names = gene_meta$gene),
      colData = DataFrame(comp_meta, row.names = comp_meta$unit_id)
    )
    file_path <- file.path(counts_dir, paste0("pseudobulk_", safe_name(compartment), ".rds"))
    saveRDS(se, file_path, compress = "gzip")
    vals <- comp_counts@x
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
      compartment = compartment,
      file = file_path,
      n_genes = nrow(comp_counts),
      n_units = ncol(comp_counts),
      n_nonzero = length(vals),
      total_library_size = sum(comp_counts),
      nonnegative = all(vals >= 0),
      integer_like = all(abs(vals - round(vals)) < 1e-8),
      file_bytes = as.numeric(file.info(file_path)$size),
      md5 = unname(tools::md5sum(file_path))
    )
    rm(se, comp_counts)
    invisible(gc())
  }
  manifest <- rbindlist(manifest_rows)
  fwrite(manifest, file.path(outdir, "pseudobulk_counts_manifest.tsv"), sep = "\t", na = "NA")
  fwrite(sample_sheet, file.path(outdir, "pseudobulk_sample_sheet.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")

  qc <- sample_sheet[, .(
    unit_id, sample_unit_key, patient_key, study, cancer_type,
    project_broad_compartment, n_cells, expected_library_size,
    observed_library_size = library_size,
    library_size_difference = library_size - expected_library_size,
    library_size_match, detected_genes, eligible_ge20, eligible_ge30, eligible_ge50
  )]
  fwrite(qc, file.path(outdir, "pseudobulk_aggregation_qc.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")

  global_expected <- sum(cell_library[non_other])
  global_observed <- sum(pb)
  counts_pass <- all(qc$library_size_match) && global_expected == global_observed &&
    all(manifest$nonnegative) && all(manifest$integer_like)
  gate <- data.table(
    gate = "G2A_pseudobulk_counts",
    status = if (counts_pass) "PASS" else "FAIL",
    n_cells_aggregated = sum(non_other),
    n_units = nrow(sample_sheet),
    n_compartments = uniqueN(sample_sheet$project_broad_compartment),
    n_genes = nrow(pb),
    global_expected_library = global_expected,
    global_observed_library = global_observed,
    n_unit_mismatches = sum(!qc$library_size_match),
    note = "Technical count aggregation gate only; no differential-expression inference was performed."
  )
  fwrite(gate, file.path(outdir, "g2a_counts_gate_summary.tsv"), sep = "\t", na = "NA")

  run_manifest <- data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version",
             "Matrix_version", "SummarizedExperiment_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string,
              as.character(packageVersion("Matrix")), as.character(packageVersion("SummarizedExperiment")), outdir)
  )
  fwrite(run_manifest, file.path(outdir, "02a_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pseudobulk count aggregation completed",
            runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
            n_units = nrow(sample_sheet), status = gate$status[[1]], output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
}

# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/07_build_pseudobulk_model_contract.R
stage_single_cell_07_build_pseudobulk_model_contract <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet", "counts_manifest", "counts_gate")]
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
`relocate_input` <- bind_stage_helper(shared_relocate_input_031, environment())


Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02b_build_pseudobulk_model_contract.R"
`log_event` <- bind_stage_helper(shared_log_event_008, environment())
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())

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
}

# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/08_run_primary_pseudobulk_models.R
stage_single_cell_08_run_primary_pseudobulk_models <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("technical_gate", "model_contract", "counts_manifest", "exclusion_genes")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir", "sample_sets_dir")])
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
`relocate_input` <- bind_stage_helper(shared_relocate_input_031, environment())


Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02c_run_primary_pseudobulk_models.R"
`log_event` <- bind_stage_helper(shared_log_event_008, environment())
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())
`select_one_per_patient` <- bind_stage_helper(shared_select_one_per_patient_005, environment())

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
}

# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/09_run_pseudobulk_stability.R
stage_single_cell_09_run_pseudobulk_stability <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("primary_gate", "model_contract", "counts_manifest")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("counts_dir", "sample_sets_dir", "primary_gene_results_dir")])
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
`relocate_input` <- bind_stage_helper(shared_relocate_input_031, environment())



Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)
script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02d_run_pseudobulk_stability.R"
`log_event` <- bind_stage_helper(shared_log_event_008, environment())
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())
`select_one_per_patient` <- bind_stage_helper(shared_select_one_per_patient_005, environment())

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
}

# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/10_summarize_pseudobulk_evidence.R
stage_single_cell_10_summarize_pseudobulk_evidence <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("primary_gate", "stability_gate", "candidates", "sensitivity", "cohort_meta", "loso", "target_stability_summary", "legacy_deg")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("primary_gene_results_dir")])
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
`relocate_input` <- bind_stage_helper(shared_relocate_input_031, environment())



suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02e_summarize_pseudobulk_evidence.R"

log_event <- function(level, msg, ...) {
  extras <- list(...)
  extra_text <- if (length(extras) == 0L) "" else paste(
    vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)),
    collapse = " "
  )
  cat(sprintf("[%s] %s script=%s msg=%s%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg),
              if (nzchar(extra_text)) " " else "", extra_text))
}
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())

main <- function() {
  started <- Sys.time()
  outdir <- out_dir
  primary_gene_results_dir <- unname(paths["primary_gene_results_dir"])
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (fread(unname(paths["primary_gate"]))$status[[1]] != "PASS") stop("G2B gate is not PASS")
  if (fread(unname(paths["stability_gate"]))$status[[1]] != "PASS") stop("G2C gate is not PASS")
  fwrite(data.table(gate = "G2D_gene_evidence_synthesis", status = "RUNNING",
                    note = "RUNNING prevents stale evidence status after interruption."),
         file.path(outdir, "g2d_evidence_gate.tsv"), sep = "\t", na = "NA")
  log_event("START", "summarize pseudobulk stability and legacy survival evidence")

  candidates <- fread(unname(paths["candidates"]))
  sensitivity <- fread(unname(paths["sensitivity"]))
  cohort_meta <- fread(unname(paths["cohort_meta"]))
  loso <- fread(unname(paths["loso"]))

  sens_sum <- sensitivity[, .(
    n_sensitivity = sum(!is.na(sensitivity_logFC)),
    sensitivity_sign_fraction = mean(direction_concordant, na.rm = TRUE),
    sensitivity_fdr010_fraction = mean(retained_fdr_010, na.rm = TRUE),
    sensitivity_median_logFC = median(sensitivity_logFC, na.rm = TRUE)
  ), by = .(comparison, compartment, gene)]

  evidence <- merge(candidates, sens_sum, by = c("comparison", "compartment", "gene"), all.x = TRUE)
  evidence <- merge(evidence, cohort_meta, by = c("comparison", "compartment", "gene"), all.x = TRUE)
  evidence <- merge(evidence, loso, by = c("comparison", "compartment", "gene"), all.x = TRUE)
  evidence[, random_direction_concordant := sign(main_logFC) == sign(beta_random)]
  evidence[, evidence_tier := fcase(
    main_fdr < 0.10 & sensitivity_sign_fraction >= 0.75 & loso_sign_fraction >= 0.80 &
      n_cohorts >= 2 & cohort_sign_fraction >= 0.67 & random_direction_concordant & p_random < 0.10,
    "A_cross_cohort_stable",
    main_fdr < 0.10 & sensitivity_sign_fraction >= 0.75 & loso_sign_fraction >= 0.80,
    "B_model_stable_context_dependent",
    main_fdr < 0.10,
    "C_overall_only",
    default = "D_leading_nonsignificant"
  )]
  evidence[, evidence_tier := factor(evidence_tier, levels = c(
    "A_cross_cohort_stable", "B_model_stable_context_dependent",
    "C_overall_only", "D_leading_nonsignificant"
  ))]
  evidence[, abs_main_logFC_sort := abs(main_logFC)]
  setorder(evidence, comparison, compartment, evidence_tier, main_fdr, -abs_main_logFC_sort)
  evidence[, abs_main_logFC_sort := NULL]
  evidence[, evidence_tier := as.character(evidence_tier)]
  fwrite(evidence, file.path(outdir, "integrated_gene_evidence.tsv.gz"),
         sep = "\t", na = "NA", compress = "gzip")

  target_summary <- evidence[, .(
    n_candidates = .N,
    n_main_fdr010_in_candidates = sum(main_fdr < 0.10, na.rm = TRUE),
    n_tier_a = sum(evidence_tier == "A_cross_cohort_stable"),
    n_tier_b = sum(evidence_tier == "B_model_stable_context_dependent"),
    n_tier_c = sum(evidence_tier == "C_overall_only"),
    n_tier_d = sum(evidence_tier == "D_leading_nonsignificant"),
    median_sensitivity_sign_fraction_main_sig = median(sensitivity_sign_fraction[main_fdr < 0.10], na.rm = TRUE),
    median_loso_sign_fraction_main_sig = median(loso_sign_fraction[main_fdr < 0.10], na.rm = TRUE),
    median_cohort_sign_fraction_main_sig = median(cohort_sign_fraction[main_fdr < 0.10], na.rm = TRUE)
  ), by = .(comparison, compartment)]
  total_signal <- fread(unname(paths["target_stability_summary"]))[,
    .(comparison, compartment, n_main_fdr010_total_nondefinition = n_main_fdr_010_nondefinition)]
  target_summary <- merge(target_summary, total_signal, by = c("comparison", "compartment"), all.x = TRUE)
  target_summary[, target_evidence_status := fcase(
    n_tier_a >= 10L, "CROSS_COHORT_SUPPORTED",
    n_tier_a > 0L | n_tier_b >= 10L, "MODEL_STABLE_CONTEXT_DEPENDENT",
    n_main_fdr010_total_nondefinition > 0L, "OVERALL_ONLY_NOT_STABLE",
    default = "NO_MAIN_FDR_SIGNAL"
  )]
  fwrite(target_summary, file.path(outdir, "target_gene_evidence_summary.tsv"), sep = "\t", na = "NA")

  # 旧 pooled-cell 证据只有 C1_vs_C4 与当前主端点完全对应；原子比较不强行映射。
  legacy_file <- unname(paths["legacy_deg"])
  legacy <- fread(legacy_file)[comparison == "C1_vs_C4" &
                                project_broad_compartment %in% c("Malignant", "T/NK", "Myeloid", "Stromal/CAF") &
                                p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25]
  legacy_rows <- list()
  for (comp in unique(legacy$project_broad_compartment)) {
    prefix <- paste0("C1_vs_C4__", safe_name(comp))
    new <- fread(relocate_input(paste0(prefix, "__main.tsv.gz"), primary_gene_results_dir))[,
      .(gene, new_logFC = logFC, new_p_value = p_value, new_fdr = fdr, definition_gene)]
    old <- legacy[project_broad_compartment == comp, .(
      comparison, compartment = project_broad_compartment, gene,
      legacy_logFC = avg_log2FC, legacy_fdr = p_val_adj
    )]
    audit <- merge(old, new, by = "gene", all.x = TRUE)
    audit[, direction_concordant := sign(legacy_logFC) == sign(new_logFC)]
    audit[, survival_class := fcase(
      is.na(new_logFC), "lost_not_tested",
      direction_concordant & new_fdr < 0.10, "retained",
      direction_concordant, "attenuated",
      !direction_concordant & new_fdr < 0.10, "reversed",
      default = "lost"
    )]
    legacy_rows[[length(legacy_rows) + 1L]] <- audit
  }
  legacy_audit <- rbindlist(legacy_rows, fill = TRUE)
  fwrite(legacy_audit, file.path(outdir, "legacy_c1c4_pooled_deg_survival.tsv.gz"),
         sep = "\t", na = "NA", compress = "gzip")
  legacy_summary <- legacy_audit[, .(
    n_legacy_significant = .N,
    n_retained = sum(survival_class == "retained"),
    n_attenuated = sum(survival_class == "attenuated"),
    n_reversed = sum(survival_class == "reversed"),
    n_lost = sum(survival_class %in% c("lost", "lost_not_tested")),
    retained_fraction = mean(survival_class == "retained"),
    direction_concordant_fraction = mean(direction_concordant, na.rm = TRUE),
    n_definition_genes = sum(definition_gene == TRUE, na.rm = TRUE)
  ), by = .(comparison, compartment)]
  fwrite(legacy_summary, file.path(outdir, "legacy_c1c4_survival_summary.tsv"), sep = "\t", na = "NA")

  out_gate <- data.table(
    gate = "G2D_gene_evidence_synthesis", status = "PASS",
    n_targets = nrow(target_summary), n_targets_cross_cohort_supported = sum(target_summary$target_evidence_status == "CROSS_COHORT_SUPPORTED"),
    n_targets_model_stable_context_dependent = sum(target_summary$target_evidence_status == "MODEL_STABLE_CONTEXT_DEPENDENT"),
    note = "Gene-level synthesis completed; pathway interpretation and final G2 decision remain pending."
  )
  fwrite(out_gate, file.path(outdir, "g2d_evidence_gate.tsv"), sep = "\t", na = "NA")
  fwrite(data.table(
    item = c("script", "started_at", "finished_at", "elapsed_minutes", "R_version", "output_directory"),
    value = c(script_name, format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              sprintf("%.2f", as.numeric(difftime(Sys.time(), started, units = "mins"))), R.version.string, outdir)
  ), file.path(outdir, "02e_run_manifest.tsv"), sep = "\t", na = "NA")
  log_event("DONE", "pseudobulk gene evidence synthesis completed", status = "PASS", output = outdir)
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  traceback(2)
  quit(save = "no", status = 1L)
})
}

# Purpose: Construct compartment pseudobulks and fit the reported patient-level comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/single_cell/11_run_pseudobulk_pathways.R
stage_single_cell_11_run_pseudobulk_pathways <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("evidence_gate", "integrated_gene_evidence", "target_gene_evidence_summary", "hallmark_gmt", "reactome_gmt", "gobp_gmt", "c7_gmt", "legacy_c7")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("primary_gene_results_dir")])
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
`relocate_input` <- bind_stage_helper(shared_relocate_input_031, environment())



Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

options(stringsAsFactors = FALSE)
set.seed(20260716)
script_name <- "LM_pan_cancer/7.1scRNA_lactate_ecosystem_rebuild/scripts/02f_run_pseudobulk_pathways_and_finalize_g2.R"
`log_event` <- bind_stage_helper(shared_log_event_008, environment())
`safe_name` <- bind_stage_helper(shared_safe_name_036, environment())

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
}
