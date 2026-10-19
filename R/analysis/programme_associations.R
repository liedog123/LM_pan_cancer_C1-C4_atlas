# Purpose: programme associations analysis stages.
# Configuration: config/drug.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Prepare and test state-adjusted programme-feature associations within cancers.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/drug_response/05_build_programme_numeric_inputs.R
stage_drug_response_05_build_programme_numeric_inputs <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient", "mrna", "program_genes", "program_clusters", "rppa", "rppa_dict", "tcpa", "tmb", "mutation", "cnv", "estimate", "cibersort", "mcpcounter", "drug")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

script_name <- "03_build_stage1_shared_numeric_inputs_v0.1_20260810.R"
started <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_015, environment())

norm_id <- function(x) substr(toupper(trimws(as.character(x))), 1L, 15L)
`norm_gene` <- bind_stage_helper(shared_norm_symbol_039, environment())

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
  # 仅构建六类关联检验使用的数值输入，保留原有标准化、患者与特征范围。
  out_dir <- unname(paths["output_dir"])
  bundle_dir <- file.path(out_dir, "bundles")
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  p <- as.list(input_paths)
  bad <- names(p)[!file.exists(unlist(p)) | file.info(unlist(p))$size <= 0]
  if (length(bad)) stop("输入缺失或为空: ", paste(bad, collapse = ","))

  log_event("START", "构建Stage 1共用数值输入", list(
    project = "LM_pan_cancer", output = out_dir,
    formal_model = "programme_feature_associations"
  ))

  patient <- fread(p$patient)
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

    bundle <- list(
      cancer = cancer,
      patient = as.data.frame(meta),
      rna_program = prog,
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
      N_C4 = sum(meta$C_subtype == "C4"),
      N_RNA_Programs = ncol(prog), N_RPPA_All = ncol(rppa),
      N_TCPA_Pathways = ncol(rppa_path), N_CNV = ncol(cnv_mat),
      N_Genomic = ncol(genome), N_Immune = ncol(immune), N_PredDrug = ncol(drug_mat),
      Bundle = out_file, Bundle_Bytes = file.info(out_file)$size
    )
    log_event("CHECKPOINT", "癌种数值bundle完成", list(
      cancer = cancer, n_samples = nrow(meta),
      n_programs = ncol(prog), n_rppa_pathways = ncol(rppa_path)
    ))
  }

  manifest <- rbindlist(manifest_rows)
  program_registry <- rbindlist(program_registry_rows)
  fwrite(manifest, file.path(out_dir, "stage1_numeric_bundle_manifest.tsv"), sep = "\t")
  fwrite(program_registry, file.path(out_dir, "stage1_rna_program_mapping_by_cancer.tsv"), sep = "\t")
  fwrite(unmapped_tcpa, file.path(out_dir, "stage1_tcpa_unmapped_component.tsv"), sep = "\t")

  if (nrow(manifest) != 22L || any(manifest$N_RNA_Programs != 9L) ||
      any(manifest$Bundle_Bytes <= 0)) stop("数值bundle QA失败")

  writeLines(c("Programme-feature numeric inputs prepared.",
               paste0("cancers=", nrow(manifest)),
               "rna_programmes_per_cancer=9",
               "consumer=06_run_programme_feature_associations.R"),
             file.path(out_dir, "programme_numeric_input_summary.txt"))

  log_event("DONE", "Stage 1共用数值输入完成", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    cancers = nrow(manifest),
    output = out_dir
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1L)
})
}

# Purpose: Prepare and test state-adjusted programme-feature associations within cancers.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/drug_response/06_run_programme_feature_associations.R
stage_drug_response_06_run_programme_feature_associations <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("numeric_manifest", "program_clusters")]
required_files <- unname(input_paths)

# Directory metadata only: do not read manifests or enumerate bundle contents.
bundle_dir <- unname(paths["bundle_dir"])
if (!dir.exists(bundle_dir) || file.access(bundle_dir, 5L) != 0L) stop("Bundle directory must be readable and searchable.")

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
package_root <- context$root
protected <- c(dirname(required_files),bundle_dir,config_path,package_root)
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

# Validate only the selected bundle at its original read position.
resolve_bundle_path <- function(value) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value) ||
      value %in% c(".", "..")) stop("Invalid manifest bundle reference.")
  name <- basename(value)
  if (!nzchar(name) || name %in% c(".", "..")) stop("Invalid bundle basename.")
  path <- resolve_path(file.path(bundle_dir, name), bundle_dir)
  if (!is_within(path, bundle_dir) || is_within(path, out_dir) ||
      !file_test("-f", path) || file.access(path, 4L) != 0L ||
      is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Bundle must be a contained readable nonempty regular file.")
  }
  path
}

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

script_name <- "04_run_stage1_trackA_directed_integration_v0.1_20260810.R"
started <- Sys.time()
`log_event` <- bind_stage_helper(shared_log_event_015, environment())

safe_cor <- function(x, y) {
  sx <- sd(x)
  sy <- sd(y)
  if (length(x) < 3L || !is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  z <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(z$estimate), p = z$p.value)
}

test_pair <- function(x, y, meta) {
  d <- data.table(
    x = as.numeric(x), y = as.numeric(y),
    state = factor(meta$C_subtype, levels = c("C4", "C1")),
    age = suppressWarnings(as.numeric(meta$Age)),
    sex = factor(meta$Sex)
  )
  base_ok <- is.finite(d$x) & is.finite(d$y) & !is.na(d$state)
  if (sum(base_ok) < 30L) return(list(Status = "NOT_EVALUABLE_N_LT30"))

  covars <- "state"
  if (mean(is.finite(d$age[base_ok])) >= 0.80 && length(unique(d$age[base_ok & is.finite(d$age)])) >= 3L) {
    covars <- c(covars, "age")
  }
  sex_ok <- base_ok & !is.na(d$sex)
  if (mean(!is.na(d$sex[base_ok])) >= 0.80 && length(unique(d$sex[sex_ok])) == 2L) {
    covars <- c(covars, "sex")
  }

  keep_cols <- c("x", "y", covars)
  cc <- complete.cases(d[, ..keep_cols])
  dd <- d[cc]
  n1 <- sum(dd$state == "C1")
  n4 <- sum(dd$state == "C4")
  if (nrow(dd) < 30L || min(n1, n4) < 10L) return(list(Status = "NOT_EVALUABLE_GROUP_GATE"))
  if (!is.finite(sd(dd$x)) || sd(dd$x) == 0) return(list(Status = "NOT_EVALUABLE_CONSTANT_PROGRAM"))
  if (!is.finite(sd(dd$y)) || sd(dd$y) == 0) return(list(Status = "NOT_EVALUABLE_CONSTANT_TARGET"))

  if ("age" %in% covars) dd[, age := as.numeric(scale(age))]
  form <- as.formula(paste("~", paste(covars, collapse = "+")))
  design <- model.matrix(form, data = dd)
  if (qr(design)$rank < ncol(design)) {
    covars <- "state"
    design <- model.matrix(~ state, data = dd)
  }
  rx <- lm.fit(design, dd$x)$residuals
  ry <- lm.fit(design, dd$y)$residuals
  partial <- safe_cor(rx, ry)
  raw <- safe_cor(dd$x, dd$y)
  if (any(!is.finite(partial)) || any(!is.finite(raw))) {
    return(list(Status = "NOT_EVALUABLE_SINGULAR_RESIDUAL"))
  }

  zx <- as.numeric(scale(dd$x))
  zy <- as.numeric(scale(dd$y))
  list(
    Status = "ESTIMATED", N = nrow(dd), N_C1 = n1, N_C4 = n4,
    Covariates = paste(covars, collapse = "+"),
    Partial_Rho = partial[["rho"]], Partial_P = partial[["p"]],
    Raw_Rho = raw[["rho"]], Raw_P = raw[["p"]],
    Program_C1_minus_C4_Z = mean(zx[dd$state == "C1"]) - mean(zx[dd$state == "C4"]),
    Target_C1_minus_C4_Z = mean(zy[dd$state == "C1"]) - mean(zy[dd$state == "C4"])
  )
}

matrix_family <- function(bundle) {
  list(
    RPPA_PATHWAY = bundle$rppa_tcpa_pathway,
    GENOMIC = cbind(bundle$cnv, bundle$genomic),
    IMMUNE_ESTIMATE = bundle$immune[, grepl("^ESTIMATE_", colnames(bundle$immune)), drop = FALSE],
    IMMUNE_CIBERSORT = bundle$immune[, grepl("^CIBERSORT_", colnames(bundle$immune)), drop = FALSE],
    IMMUNE_MCPCOUNTER = bundle$immune[, grepl("^MCPCOUNTER_", colnames(bundle$immune)), drop = FALSE],
    PREDICTED_DRUG = bundle$predicted_drug
  )
}

main <- function() {
  out_dir <- unname(paths["output_dir"])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- unname(paths["numeric_manifest"])
  cluster_path <- unname(paths["program_clusters"])
  if (!file.exists(manifest_path) || !file.exists(cluster_path)) stop("缺少共用输入或cluster registry")
  manifest <- fread(manifest_path)
  clusters <- fread(cluster_path)[State %in% c("C1", "C4") & Cluster_ID != "OTHER",
                                  .(Program = Cluster_ID, Program_Label = Cluster_Label,
                                    Program_State = State, Semantic_Boundary)]
  if (nrow(manifest) != 22L || nrow(clusters) != 9L) stop("输入Gate不一致")

  log_event("START", "线路A逐癌定向整合开始", list(
    project = "LM_pan_cancer", cancers = nrow(manifest), programs = nrow(clusters),
    method = "partial_spearman_residualized_state_age_sex", output = out_dir
  ))

  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    cancer <- manifest$Cancer[i]
    # Historical manifest paths provide basenames; keep bundles inside the configured directory.
    b <- readRDS(resolve_bundle_path(manifest$Bundle[i]))
    meta <- as.data.table(b$patient)
    ids <- meta$Patient_ID
    if (!identical(rownames(b$rna_program), ids)) stop(cancer, " RNA程序患者顺序不一致")
    families <- matrix_family(b)
    for (family in names(families)) {
      target <- families[[family]]
      if (is.null(dim(target)) || nrow(target) != length(ids)) stop(cancer, " target矩阵维度异常: ", family)
      for (program in colnames(b$rna_program)) {
        for (feature in colnames(target)) {
          ans <- test_pair(b$rna_program[, program], target[, feature], meta)
          rows[[length(rows) + 1L]] <- data.table(
            Cancer = cancer, Program = program, Evidence_Family = family,
            Target_Feature = feature,
            Status = ans$Status,
            N = if (!is.null(ans$N)) ans$N else NA_integer_,
            N_C1 = if (!is.null(ans$N_C1)) ans$N_C1 else NA_integer_,
            N_C4 = if (!is.null(ans$N_C4)) ans$N_C4 else NA_integer_,
            Covariates = if (!is.null(ans$Covariates)) ans$Covariates else NA_character_,
            Partial_Rho = if (!is.null(ans$Partial_Rho)) ans$Partial_Rho else NA_real_,
            Partial_P = if (!is.null(ans$Partial_P)) ans$Partial_P else NA_real_,
            Raw_Rho = if (!is.null(ans$Raw_Rho)) ans$Raw_Rho else NA_real_,
            Raw_P = if (!is.null(ans$Raw_P)) ans$Raw_P else NA_real_,
            Program_C1_minus_C4_Z = if (!is.null(ans$Program_C1_minus_C4_Z)) ans$Program_C1_minus_C4_Z else NA_real_,
            Target_C1_minus_C4_Z = if (!is.null(ans$Target_C1_minus_C4_Z)) ans$Target_C1_minus_C4_Z else NA_real_
          )
        }
      }
    }
    log_event("CHECKPOINT", "线路A癌种完成", list(cancer = cancer, n_samples = nrow(meta)))
  }

  res <- rbindlist(rows, fill = TRUE)
  res <- merge(res, clusters, by = "Program", all.x = TRUE)
  res[, Partial_FDR := NA_real_]
  res[Status == "ESTIMATED", Partial_FDR := p.adjust(Partial_P, method = "BH"),
      by = .(Cancer, Evidence_Family)]
  setcolorder(res, c("Cancer", "Program", "Program_Label", "Program_State", "Semantic_Boundary",
                     "Evidence_Family", "Target_Feature", "Status", "N", "N_C1", "N_C4",
                     "Covariates", "Partial_Rho", "Partial_P", "Partial_FDR", "Raw_Rho", "Raw_P",
                     "Program_C1_minus_C4_Z", "Target_C1_minus_C4_Z"))

  summary <- res[Status == "ESTIMATED", {
    sig <- is.finite(Partial_FDR) & Partial_FDR < 0.05
    n_sig <- sum(sig)
    n_pos_sig <- sum(sig & Partial_Rho > 0)
    n_neg_sig <- sum(sig & Partial_Rho < 0)
    dom_n <- max(n_pos_sig, n_neg_sig)
    .(
      N_Evaluable_Cancers = .N,
      N_Positive_All = sum(Partial_Rho > 0, na.rm = TRUE),
      N_Negative_All = sum(Partial_Rho < 0, na.rm = TRUE),
      Median_Partial_Rho = median(Partial_Rho, na.rm = TRUE),
      N_FDR_Significant = n_sig,
      N_Positive_Significant = n_pos_sig,
      N_Negative_Significant = n_neg_sig,
      Dominant_Significant_Direction = if (n_sig == 0L) "NONE" else if (n_pos_sig >= n_neg_sig) "POSITIVE" else "NEGATIVE",
      Significant_Direction_Consistency = if (n_sig == 0L) NA_real_ else dom_n / n_sig,
      Significant_Cancers = paste(sort(Cancer[sig]), collapse = ";")
    )
  }, by = .(Program, Program_Label, Program_State, Evidence_Family, Target_Feature)]

  orthogonal <- summary[
    Evidence_Family %in% c("RPPA_PATHWAY", "GENOMIC") &
      N_FDR_Significant >= 3L & Significant_Direction_Consistency >= 0.70
  ]
  orthogonal[, Abs_Median_Partial_Rho := abs(Median_Partial_Rho)]
  setorder(orthogonal, -N_FDR_Significant, -Significant_Direction_Consistency,
           -Abs_Median_Partial_Rho, Program, Evidence_Family, Target_Feature)
  orthogonal[, Abs_Median_Partial_Rho := NULL]

  res_path <- file.path(out_dir, "stage1_trackA_associations_by_cancer.tsv.gz")
  summary_path <- file.path(out_dir, "stage1_trackA_cross_cancer_summary.tsv")
  candidate_path <- file.path(out_dir, "stage1_trackA_orthogonal_candidates.tsv")
  fwrite(res, res_path, sep = "\t", compress = "gzip")
  fwrite(summary, summary_path, sep = "\t")
  fwrite(orthogonal, candidate_path, sep = "\t")

  qa <- data.table(
    Check = c(
      "22 cancers processed", "9 programs present", "no pooled patient model",
      "all estimated p finite", "within-family BH reproducible",
      "candidate restricted to orthogonal families", "output files nonempty"
    ),
    Pass = c(
      uniqueN(res$Cancer) == 22L,
      uniqueN(res$Program) == 9L,
      TRUE,
      all(is.finite(res[Status == "ESTIMATED", Partial_P])),
      res[Status == "ESTIMATED", all(abs(Partial_FDR - p.adjust(Partial_P, "BH")) < 1e-12),
          by = .(Cancer, Evidence_Family)][, all(V1)],
      all(orthogonal$Evidence_Family %in% c("RPPA_PATHWAY", "GENOMIC")),
      all(file.exists(c(res_path, summary_path, candidate_path))) &&
        all(file.info(c(res_path, summary_path, candidate_path))$size > 0)
    )
  )
  fwrite(qa, file.path(out_dir, "stage1_trackA_qa.tsv"), sep = "\t")
  if (!all(qa$Pass)) stop("线路A机器QA未全部通过")

  top_text <- if (nrow(orthogonal)) {
    paste(sprintf("%s—%s（%s，%s癌FDR通过）",
                  head(orthogonal$Program, 5), head(orthogonal$Target_Feature, 5),
                  head(orthogonal$Dominant_Significant_Direction, 5),
                  head(orthogonal$N_FDR_Significant, 5)), collapse = "；")
  } else "无候选达到预设的至少3癌种、方向一致性至少70%门槛"

  report <- c(
    "# Stage 1线路A定向整合报告 V0.1", "",
    "日期：2026-08-10  ",
    "状态：`TRACK_A_COMPLETE / MACHINE_QA_PASS / BIOLOGICAL_REVIEW_PENDING`", "",
    sprintf("共完成%s条逐癌程序—目标关联，其中%s条可评价。", nrow(res), res[Status == "ESTIMATED", .N]),
    sprintf("按单癌种、单证据家族BH后有%s条关联FDR<0.05。", res[Status == "ESTIMATED" & Partial_FDR < 0.05, .N]),
    sprintf("相对正交的RPPA/基因组候选单元为%s条。", nrow(orthogonal)),
    paste0("排序靠前的候选：", top_text, "。"), "",
    "主统计已调整C1/C4、年龄和性别，因此检验的是组别之外的跨层协同；它仍是观察性关联，不构成因果或药物机制。",
    "TCGA预测药敏和免疫结果保留为RNA耦联注释，未用于多组学候选的单独晋级。"
  )
  writeLines(report, file.path(out_dir, "stage1_trackA_report_v0.1.md"), useBytes = TRUE)

  log_event("DONE", "线路A定向整合完成", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    n_tests = nrow(res), n_estimated = res[Status == "ESTIMATED", .N],
    n_fdr = res[Status == "ESTIMATED" & Partial_FDR < 0.05, .N],
    n_orthogonal_candidates = nrow(orthogonal), output = out_dir
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1L)
})
}
