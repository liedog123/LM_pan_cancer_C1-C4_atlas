#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 06_prepare_liu_exploratory.R --config PATH [--check-config]\n       Rscript 06_prepare_liu_exploratory.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("source_rda","official_s2_tpm","clinical_workbook","gene_alias","lm12_members","program_membership","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("source_rda","official_s2_tpm","clinical_workbook","gene_alias","lm12_members","program_membership")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),input_dirs,config_path,package_root)
if(cfg_out_dir=="/" || cfg_out_dir==normalizePath(path.expand("~")) ||
   any(vapply(protected,function(p) is_within(p,cfg_out_dir)||is_within(cfg_out_dir,p),logical(1L)))) stop("Output must be independent of inputs/config.")
if(file.exists(cfg_out_dir) && !dir.exists(cfg_out_dir)) stop("Output is not a directory.")
if(dir.exists(cfg_out_dir) && (file.access(cfg_out_dir,7L)!=0L || length(list.files(cfg_out_dir,all.files=TRUE,no..=TRUE)))) stop("Output must be writable and empty.")
parent <- cfg_out_dir
while(!dir.exists(parent)) {
 if(file.exists(parent)) stop("Output has a non-directory ancestor.")
 parent <- dirname(parent)
}
if(file.access(parent,3L)!=0L) stop("No writable searchable output parent.")
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

cfg_paths <- paths
# 普通运行才创建独立空输出；检查模式已退出。
dir.create(unname(cfg_paths["output_dir"]),recursive=TRUE,showWarnings=FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(readxl)
  library(SummarizedExperiment)
})

# Private project root replaced by explicit DCF inputs.
# Module paths are supplied by the checked DCF.
# Prepared source files are independently configured.
out_dir <- unname(cfg_paths["output_dir"])
source_rda <- unname(cfg_paths["source_rda"])
source_tpm <- unname(cfg_paths["official_s2_tpm"])
source_clin <- unname(cfg_paths["clinical_workbook"])

log_event <- function(level, msg) {
  cat(sprintf("[%s] %s script=06_prepare_liu_exploratory.R msg=\"%s\"\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
  flush.console()
}
log_event("START", "tigeR full-gene FPKM conversion and outcome-blind assignment; primary 103, timing sensitivity 119")
required <- c(source_rda, source_tpm, source_clin,
              unname(cfg_paths["gene_alias"]),
              unname(cfg_paths["lm12_members"]),
              unname(cfg_paths["program_membership"]))
if (any(!file.exists(required))) stop("Missing required input: ", paste(required[!file.exists(required)], collapse = " | "))
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L)
  stop("Output directory is non-empty; refusing to overwrite: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_col <- function(rng) readxl::read_xlsx(source_clin, sheet = "Supplemental Table 1", range = rng, col_names = FALSE)
patient_id <- trimws(as.character(read_col("A4:A161")[[1]]))
months_pre <- suppressWarnings(as.numeric(read_col("K4:K161")[[1]]))
tx <- trimws(as.character(read_col("R4:R161")[[1]]))
prior_ctla4 <- suppressWarnings(as.integer(read_col("AJ4:AJ161")[[1]]))
ctx_days <- read_col("AP4:AQ161")
biopsy_context <- suppressWarnings(as.integer(ctx_days[[1]]))
days_to_pd1 <- suppressWarnings(as.numeric(ctx_days[[2]]))
site <- trimws(as.character(read_col("J4:J161")[[1]]))
meta_all <- data.table(patient_id = patient_id, monthsBiopsyPreTx = months_pre, Tx = tx,
                       priorCTLA4 = prior_ctla4, biopsyContext = biopsy_context,
                       daysBiopsyToPD1 = days_to_pd1, biopsy_site = site)
meta_all <- meta_all[!is.na(patient_id) & nzchar(patient_id)]
if (anyDuplicated(meta_all$patient_id)) stop("Table S1 patient IDs are not unique")
meta_all[, timepoint_label := fcase(biopsyContext == 1L, "Pre-Ipi", biopsyContext == 2L, "On-Ipi",
                                    biopsyContext == 3L, "Pre-PD1", biopsyContext == 4L, "On-PD1",
                                    default = "Unknown")]

# RDA 中只抽取样本身份两列；结局列不参与本脚本的状态构建。
e <- new.env(parent = globalenv())
loaded <- load(source_rda, envir = e)
if (!"MEL_phs000452" %in% loaded) stop("Author RDA did not contain MEL_phs000452")
obj <- e$MEL_phs000452
fpkm0 <- as.matrix(SummarizedExperiment::assay(obj))
storage.mode(fpkm0) <- "double"
cd <- as.data.frame(SummarizedExperiment::colData(obj)[, c("sample_id", "patient_name"), drop = FALSE])
if (!all(c("sample_id", "patient_name") %in% names(cd))) stop("tigeR colData is missing sample_id/patient_name")
source_map <- data.table(sample_id = as.character(cd$sample_id), patient_name = as.character(cd$patient_name))
if (!setequal(colnames(fpkm0), source_map$sample_id)) stop("FPKM assay columns do not match tigeR sample_id")
source_map[, patient_id := sub("_.*$", "", patient_name)]
source_map <- source_map[grepl("^Patient[0-9]+$", patient_id)]
if (anyDuplicated(source_map$sample_id) || anyDuplicated(source_map$patient_id))
  stop("tigeR Patient### keys are not one-to-one with samples")
if (nrow(fpkm0) != 56269L || nrow(source_map) != 121L) stop("Unexpected tigeR object dimensions or Patient### count")
if (any(!is.finite(fpkm0)) || any(fpkm0 < 0)) stop("tigeR FPKM contains non-finite or negative values")
if (!setequal(rownames(fpkm0), rownames(obj))) stop("Assay gene rownames mismatch")

# 官方 S2 只用于 Patient ID 交集和表达量尺度对照，不作为状态评分矩阵。
official <- fread(source_tpm, check.names = FALSE, showProgress = FALSE)
if (nrow(official) != 121L || ncol(official) != 20849L) stop("Unexpected official S2 dimensions")
official_ids <- as.character(official[[1]])
official_gene <- toupper(trimws(names(official)[-1]))
if (anyDuplicated(official_ids) || anyDuplicated(official_gene)) stop("Duplicate official S2 sample or gene IDs")
official_tpm <- as.matrix(official[, -1, with = FALSE])
storage.mode(official_tpm) <- "double"
rownames(official_tpm) <- official_ids
colnames(official_tpm) <- official_gene
if (any(!is.finite(official_tpm)) || any(official_tpm < 0)) stop("Official S2 TPM contains invalid values")

# 完整 tigeR FPKM 按每样本全部保留基因换算 TPM；不按 S2 过滤后的 20,848 行再次缩放。
gene_raw <- toupper(trimws(rownames(fpkm0)))
if (anyNA(gene_raw) || any(!nzchar(gene_raw)) || anyDuplicated(gene_raw)) stop("tigeR gene symbols are missing or duplicated")
rownames(fpkm0) <- gene_raw
fpkm_sum_raw <- colSums(fpkm0)
if (any(fpkm_sum_raw <= 0)) stop("tigeR has a zero-sum FPKM sample")
tpm_for_source_check <- sweep(fpkm0, 2, fpkm_sum_raw, "/") * 1e6
common_gene <- intersect(rownames(tpm_for_source_check), colnames(official_tpm))
common_patient <- sort(intersect(source_map$patient_id, official_ids))
if (length(common_patient) != 120L || length(common_gene) < 10000L)
  stop(sprintf("Unexpected source concordance dimensions: patients=%d genes=%d", length(common_patient), length(common_gene)))
cor_row <- function(local_id, official_id, comparison) {
  a <- log2(tpm_for_source_check[common_gene, source_map$sample_id[match(local_id, source_map$patient_id)]] + 1)
  b <- log2(official_tpm[official_id, common_gene] + 1)
  data.table(comparison = comparison, local_patient_id = local_id, official_S2_patient_id = official_id,
             shared_gene_n = length(common_gene), pearson_log2TPM = cor(a, b, method = "pearson"),
             spearman_log2TPM = cor(a, b, method = "spearman"))
}
matched_cor <- rbindlist(lapply(common_patient, function(id) cor_row(id, id, "exact_patient_key")))
shift_ids <- c(common_patient[-1L], common_patient[1L])
shift_cor <- rbindlist(lapply(seq_along(common_patient), function(i)
  cor_row(common_patient[i], shift_ids[i], "fixed_one_position_shift_control")))
patient192_test <- if ("Patient192" %in% source_map$patient_id && "Patient20" %in% official_ids)
  cor_row("Patient192", "Patient20", "unmapped_candidate_Patient192_vs_Patient20") else
  data.table(comparison = "unmapped_candidate_Patient192_vs_Patient20", local_patient_id = "Patient192",
             official_S2_patient_id = "Patient20", shared_gene_n = length(common_gene),
             pearson_log2TPM = NA_real_, spearman_log2TPM = NA_real_)
concordance <- rbindlist(list(matched_cor, shift_cor, patient192_test), fill = TRUE)
source_crosswalk <- data.table(official_S2_patient_id = sort(official_ids))
source_crosswalk[, exact_tigeR_sample_match := source_map$sample_id[match(official_S2_patient_id, source_map$patient_id)]]
source_crosswalk[, exact_tigeR_patient_key_match := !is.na(exact_tigeR_sample_match)]
source_crosswalk[, mapping_status := fifelse(exact_tigeR_patient_key_match, "exact Patient### prefix match", "no exact tigeR Patient### key")]
source_extra <- source_map[!patient_id %in% official_ids,
                           .(official_S2_patient_id = NA_character_, exact_tigeR_sample_match = sample_id,
                             exact_tigeR_patient_key_match = FALSE, mapping_status = paste0("tigeR-only exact key: ", patient_id))]
source_crosswalk <- rbindlist(list(source_crosswalk, source_extra), fill = TRUE)
if (sum(source_crosswalk$exact_tigeR_patient_key_match) != 120L) stop("Expected 120 exact tigeR/S2 Patient### matches")
if (any(source_crosswalk$official_S2_patient_id == "Patient20" & source_crosswalk$exact_tigeR_patient_key_match)) stop("Patient20 unexpectedly mapped")

# 临床主集合与宽时间窗先按 Table S1 定义，再要求 tigeR 有精确患者键。
meta_s2 <- meta_all[patient_id %in% official_ids]
if (any(meta_s2$biopsyContext == 3L & meta_s2$daysBiopsyToPD1 > 0, na.rm = TRUE))
  stop("Positive day offset conflicts with Table S1 Pre-PD1 context")
if (any(meta_s2$biopsyContext == 4L & meta_s2$daysBiopsyToPD1 <= 0, na.rm = TRUE))
  stop("On-PD1 context conflicts with nonpositive day offset")
sample_map <- source_map[, .(patient_id, sample_id, patient_name)]
manifest <- merge(meta_all, sample_map, by = "patient_id", all.x = TRUE, sort = FALSE)
manifest[, has_exact_tigeR_patient_key := !is.na(sample_id)]
manifest[, has_official_S2_TPM_patient := patient_id %in% official_ids]
manifest[, official_primary_timepoint := has_official_S2_TPM_patient & !is.na(biopsyContext) & biopsyContext == 3L & !is.na(daysBiopsyToPD1) & daysBiopsyToPD1 <= 0]
manifest[, official_broad_timepoint := has_official_S2_TPM_patient & !is.na(biopsyContext) & biopsyContext %in% c(1L, 2L, 3L) & !is.na(daysBiopsyToPD1) & daysBiopsyToPD1 <= 0]
manifest[, include_primary103 := official_primary_timepoint & has_exact_tigeR_patient_key]
manifest[, include_broad119 := official_broad_timepoint & has_exact_tigeR_patient_key]
manifest[, primary_reason := fcase(!has_official_S2_TPM_patient, "excluded_not_in_official_S2_expression_cohort",
                                    !official_primary_timepoint, "not_Pre-PD1_context3_day_le_0",
                                    !has_exact_tigeR_patient_key, "excluded_no_exact_tigeR_patient_key",
                                    default = "included_primary_Pre-PD1_context3")]
manifest[, broad_reason := fcase(!has_official_S2_TPM_patient, "excluded_not_in_official_S2_expression_cohort",
                                  !official_broad_timepoint, "not_context1to3_day_le_0",
                                  !has_exact_tigeR_patient_key, "excluded_no_exact_tigeR_patient_key",
                                  default = "included_broad_pre-PD1_context1to3")]
manifest[, day0_prePD1_context := biopsyContext == 3L & daysBiopsyToPD1 == 0]
primary_official <- manifest[official_primary_timepoint == TRUE]
broad_official <- manifest[official_broad_timepoint == TRUE]
primary_meta <- manifest[include_primary103 == TRUE]
broad_meta <- manifest[include_broad119 == TRUE]
if (nrow(primary_official) != 104L || nrow(broad_official) != 120L || nrow(primary_meta) != 103L || nrow(broad_meta) != 119L)
  stop(sprintf("Unexpected eligible/mapped set sizes: official primary=%d broad=%d; mapped primary=%d broad=%d",
               nrow(primary_official), nrow(broad_official), nrow(primary_meta), nrow(broad_meta)))
if (!setequal(primary_official$patient_id[!primary_official$has_exact_tigeR_patient_key], "Patient20"))
  stop("The sole unmapped primary patient is not Patient20")
if (sum(primary_meta$daysBiopsyToPD1 == 0, na.rm = TRUE) != 2L) stop("Unexpected mapped day-zero primary count")

alias <- fread(unname(cfg_paths["gene_alias"]))
alias_map <- setNames(toupper(alias$canonical_symbol), toupper(alias$alias_symbol))
normalize_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- ""
  hit <- x %in% names(alias_map)
  x[hit] <- unname(alias_map[x[hit]])
  x
}
gene <- normalize_symbol(rownames(fpkm0))
keep_gene <- nzchar(gene)
fpkm <- fpkm0[keep_gene, , drop = FALSE]
rownames(fpkm) <- gene[keep_gene]
duplicate_gene_rows <- sum(duplicated(rownames(fpkm)))
if (duplicate_gene_rows > 0L) {
  n_gene_rows <- rowsum(matrix(1, nrow(fpkm), ncol(fpkm)), group = rownames(fpkm), reorder = FALSE)
  fpkm <- rowsum(fpkm, group = rownames(fpkm), reorder = FALSE) / n_gene_rows
}
fpkm_sum <- colSums(fpkm)
if (any(fpkm_sum <= 0)) stop("Zero-sum FPKM sample after gene-symbol normalization")
tpm <- sweep(fpkm, 2, fpkm_sum, "/") * 1e6
logexpr <- log2(tpm + 1)

lm_long <- fread(unname(cfg_paths["lm12_members"]))
lm_sets <- lapply(split(normalize_symbol(lm_long$gene), lm_long$feature_id), unique)
lm_sets <- lm_sets[sort(names(lm_sets))]
if (length(lm_sets) != 12L || !"LM12_Antimicrobial_Peptides" %in% names(lm_sets)) stop("Unexpected frozen LM12 signature catalog")
acid <- c("ATP6AP1", "GAPDH", "HIF1A", "HK2", "PFKP", "PGK1", "PKM", "SLC2A1", "SLC9A7")
pd <- fread(unname(cfg_paths["program_membership"]))
prog_ids <- c(Hypoxia = "MALIGNANT_HYPOXIA", ECM = "STROMAL_ECM_ORGANIZATION", Cytotoxic = "TNK_CYTOTOXICITY_EFFECTOR")
programs <- lapply(prog_ids, function(id) sort(unique(normalize_symbol(pd[program_id == id & include_in_definition_independent_score == TRUE]$gene))))
if (length(programs$Cytotoxic) != 6L) stop("Frozen cytotoxic program does not contain the expected six members")

assign_one_set <- function(set_name, md) {
  setorder(md, patient_id)
  ids <- md$patient_id
  sample_ids <- md$sample_id
  expr <- logexpr[, sample_ids, drop = FALSE]
  colnames(expr) <- ids
  sets <- c(list(Acid_Score = acid), lm_sets)
  coverage <- list()
  present_sets <- list()
  for (nm in names(sets)) {
    gs <- unique(sets[[nm]])
    seen <- intersect(gs, rownames(expr))
    vv <- expr[seen, , drop = FALSE]
    sdv <- if (nrow(vv) > 0L) apply(vv, 1, sd) else numeric()
    variable <- names(sdv)[is.finite(sdv) & sdv > 0]
    pass <- if (nm == "Acid_Score") length(seen) == 9L && length(variable) == 9L else length(seen) / length(gs) >= 0.70 && length(variable) > 0L
    coverage[[length(coverage) + 1L]] <- data.table(sample_set = set_name, score = nm, frozen_n = length(gs),
      present_n = length(seen), variable_n = length(variable), coverage_pass = pass,
      missing = paste(setdiff(gs, seen), collapse = "|"), zero_variance = paste(setdiff(seen, variable), collapse = "|"))
    present_sets[[nm]] <- seen
  }
  cov <- rbindlist(coverage)
  if (!all(cov$coverage_pass)) stop(set_name, ": fixed Acid/LM12 gene coverage failed")
  ss <- GSVA::gsva(GSVA::ssgseaParam(exprData = expr, geneSets = present_sets, normalize = TRUE), verbose = FALSE)
  lm <- ss[names(lm_sets), , drop = FALSE]
  axis <- colMeans(lm[setdiff(rownames(lm), "LM12_Antimicrobial_Peptides"), , drop = FALSE])
  acid_score <- as.numeric(ss["Acid_Score", ]); names(acid_score) <- colnames(ss)

  # 使用已冻结的九基因 z-score + ward.D/k=3 和 LM12 ssGSEA + ward.D/k=2。
  az <- expr[acid, , drop = FALSE]
  az <- sweep(sweep(az, 1, rowMeans(az), "-"), 1, apply(az, 1, sd), "/")
  ca <- cutree(hclust(dist(t(az)), method = "ward.D"), k = 3)
  means_a <- tapply(acid_score, ca, mean)
  if (anyDuplicated(as.numeric(means_a))) stop(set_name, ": Acid cluster means tie; AG direction unresolved")
  ag_map <- setNames(paste0("AG", 1:3), names(sort(means_a, decreasing = TRUE)))
  ci <- cutree(hclust(dist(t(lm)), method = "ward.D"), k = 2)
  means_i <- tapply(axis, ci, mean)
  if (anyDuplicated(as.numeric(means_i))) stop(set_name, ": LM12 cluster means tie; IG direction unresolved")
  ig_map <- setNames(paste0("IG", 1:2), names(sort(means_i, decreasing = TRUE)))
  ans <- data.table(patient_id = ids, sample_id = sample_ids, AG = unname(ag_map[as.character(ca)]),
    IG = unname(ig_map[as.character(ci)]), Acid_Score = acid_score[ids], LM12_inflammatory_axis = axis[ids])
  ans[, C_state := fcase(AG == "AG1" & IG == "IG1", "C1", AG == "AG3" & IG == "IG2", "C4",
                         AG == "AG1" & IG == "IG2", "C2", AG %in% c("AG2", "AG3") & IG == "IG1", "C3",
                         AG == "AG2" & IG == "IG2", "C2")]
  if (anyNA(ans$C_state)) stop(set_name, ": C1-C4 map incomplete")
  for (nm in names(programs)) {
    gs <- unique(programs[[nm]])
    seen <- intersect(gs, rownames(expr))
    vv <- expr[seen, , drop = FALSE]
    sdv <- if (nrow(vv) > 0L) apply(vv, 1, sd) else numeric()
    variable <- names(sdv)[is.finite(sdv) & sdv > 0]
    pass <- if (nm == "Cytotoxic") length(seen) == 6L && length(variable) >= 5L else length(seen) / length(gs) >= 0.80 && length(variable) >= 10L
    cov <- rbind(cov, data.table(sample_set = set_name, score = nm, frozen_n = length(gs), present_n = length(seen),
      variable_n = length(variable), coverage_pass = pass, missing = paste(setdiff(gs, seen), collapse = "|"),
      zero_variance = paste(setdiff(seen, variable), collapse = "|")))
    val <- rep(NA_real_, length(ids)); names(val) <- ids
    if (pass) {
      use <- vv[variable, , drop = FALSE]
      val <- colMeans(sweep(sweep(use, 1, rowMeans(use), "-"), 1, sdv[variable], "/"))
      if (!is.finite(sd(val)) || sd(val) == 0) pass <- FALSE
      else {
        val <- as.numeric(scale(val))
        names(val) <- ids
      }
    }
    if (!pass) val[] <- NA_real_
    ans[, (nm) := val[ids]]
  }
  scored <- merge(ans, md, by = c("patient_id", "sample_id"), all.x = TRUE, sort = FALSE)
  setorder(scored, patient_id)
  scored[, sample_set := set_name]
  list(data = scored, coverage = cov)
}

res103 <- assign_one_set("primary_Pre_PD1_context3_n103", primary_meta)
res119 <- assign_one_set("sensitivity_all_pre_PD1_context1to3_n119", broad_meta)
state_stability <- merge(res103$data[, .(patient_id, primary_state = C_state)],
                         res119$data[, .(patient_id, broad_state = C_state)], by = "patient_id")
state_table <- as.data.table(as.table(table(factor(state_stability$primary_state, levels = c("C1", "C2", "C3", "C4")),
                                             factor(state_stability$broad_state, levels = c("C1", "C2", "C3", "C4")))))
setnames(state_table, c("primary_state", "broad_state", "n"))
stability_summary <- data.table(n_shared = nrow(state_stability), n_same_state = sum(state_stability$primary_state == state_stability$broad_state),
  exact_agreement = mean(state_stability$primary_state == state_stability$broad_state),
  n_primary_day0 = sum(primary_meta$daysBiopsyToPD1 == 0, na.rm = TRUE),
  n_primary_day_negative = sum(primary_meta$daysBiopsyToPD1 < 0, na.rm = TRUE),
  n_broad_Pre_Ipi = sum(broad_meta$biopsyContext == 1L), n_broad_On_Ipi = sum(broad_meta$biopsyContext == 2L),
  n_broad_Pre_PD1 = sum(broad_meta$biopsyContext == 3L))

audit <- rbindlist(list(
  data.table(metric = c("tigeR_source_object", "source_expression_scale", "analysis_expression_scale",
                        "official_S2_use", "TPM_scaling_universe", "official_S2_exact_patient_matches",
                        "official_S2_patients_without_exact_tigeR_match", "tigeR_only_patient_key",
                        "shared_exact_gene_symbols_for_concordance", "matched_log2TPM_Pearson_median",
                        "matched_log2TPM_Pearson_P10", "matched_log2TPM_Spearman_median",
                        "shift_control_log2TPM_Pearson_median", "Patient192_vs_Patient20_log2TPM_Pearson",
                        "Table_S1_primary_context3_day_le0", "primary_103_actual", "Table_S1_broad_context1to3_day_le0",
                        "broad_119_actual", "primary_day0_official", "primary_day0_actual", "duplicate_gene_rows_after_alias",
                        "GSVA_version", "ssGSEA_normalize"),
             value = c("YuLab-SMU/tigeR.data MEL_phs000452.rda", "FPKM", "FPKM/sum(FPKM over full 56,269-gene object)*1e6 then log2(TPM+1)",
                       "patient ID and expression-scale concordance only; not used for state scores", "all tigeR gene rows after frozen symbol alias aggregation by mean",
                       sum(source_crosswalk$exact_tigeR_patient_key_match), "Patient20", "Patient192", length(common_gene),
                       median(matched_cor$pearson_log2TPM), unname(quantile(matched_cor$pearson_log2TPM, .10)),
                       median(matched_cor$spearman_log2TPM), median(shift_cor$pearson_log2TPM), patient192_test$pearson_log2TPM[1],
                       nrow(primary_official), nrow(primary_meta), nrow(broad_official), nrow(broad_meta),
                       sum(primary_official$daysBiopsyToPD1 == 0, na.rm = TRUE), sum(primary_meta$daysBiopsyToPD1 == 0, na.rm = TRUE),
                       duplicate_gene_rows, as.character(packageVersion("GSVA")), "TRUE"))), fill = TRUE)
tx_summary <- primary_meta[, .N, by = Tx][order(Tx)]
tx_summary[, analysis_set := "primary_103"]
prior_summary <- primary_meta[, .N, by = priorCTLA4][order(priorCTLA4)]
prior_summary[, analysis_set := "primary_103"]

fwrite(source_crosswalk, file.path(out_dir, "Liu_tigeR_S2_exact_patient_crosswalk.tsv"), sep = "\t", na = "NA")
fwrite(concordance, file.path(out_dir, "Liu_tigeR_S2_expression_concordance.tsv"), sep = "\t", na = "NA")
fwrite(manifest, file.path(out_dir, "Liu_103_119_patient_inclusion_manifest.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(list(res103$data, res119$data), fill = TRUE), file.path(out_dir, "Liu_103_119_outcome_blind_state_assignments.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(list(res103$coverage, res119$coverage)), file.path(out_dir, "Liu_103_119_member_coverage.tsv"), sep = "\t", na = "NA")
fwrite(state_table, file.path(out_dir, "Liu_103_vs_119_assignment_crosstab.tsv"), sep = "\t", na = "NA")
fwrite(stability_summary, file.path(out_dir, "Liu_103_vs_119_assignment_stability.tsv"), sep = "\t", na = "NA")
fwrite(audit, file.path(out_dir, "Liu_tigeR_input_and_source_audit.tsv"), sep = "\t", na = "NA")
fwrite(tx_summary, file.path(out_dir, "Liu_103_Tx_levels.tsv"), sep = "\t", na = "NA")
fwrite(prior_summary, file.path(out_dir, "Liu_103_priorCTLA4_levels.tsv"), sep = "\t", na = "NA")

log_event("CHECKPOINT", sprintf("source mapping: 120/121 exact Patient### IDs; Patient20 has no exact expression key; Patient192 remains unmapped"))
log_event("DONE", sprintf("outcome-blind assignment complete: primary=%d/%d; broad=%d/%d; shared state agreement=%.3f; output=%s",
                           nrow(res103$data), nrow(primary_official), nrow(res119$data), nrow(broad_official),
                           stability_summary$exact_agreement, out_dir))
