# Purpose: cd8 response analysis stages.
# Configuration: config/cd8.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Prepare CD8 sample metadata, filtered expression and PCA summaries.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/01_qc_and_pca.R
stage_cd8_perturbation_01_qc_and_pca <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["count_rdata"]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# ---- Retained complete source body; CLI/path substitutions only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

input_file <- unname(paths["count_rdata"])
out_dir <- unname(paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
`write_tsv` <- bind_stage_helper(shared_write_tsv_053, environment())
`stop_if` <- bind_stage_helper(shared_stop_if_004, environment())

stop_if(!requireNamespace("edgeR", quietly = TRUE), "R 包 edgeR 未安装。")

e <- new.env(parent = emptyenv())
stop_if(!file.exists(input_file), paste("输入文件不存在:", input_file))
load(input_file, envir = e)
objects <- ls(e, all.names = TRUE)
stop_if(!all(c("count", "sample_anno") %in% objects),
        "RData 必须同时包含 count 和 sample_anno。")
count <- get("count", envir = e)
sample_anno <- get("sample_anno", envir = e)
stop_if(!identical(colnames(count), rownames(sample_anno)),
        "count 列名与 sample_anno 行名不一致。")

anno <- sample_anno[colnames(count), c("Donor", "Treatment", "Time"), drop = FALSE]
anno$SampleID <- rownames(anno)
factorial <- as.character(anno$Treatment) %in% c("HCl", "LA", "NC")
fa <- anno[factorial, , drop = FALSE]
fa$Donor <- factor(as.character(fa$Donor))
fa$Treatment <- factor(as.character(fa$Treatment), levels = c("NC", "HCl", "LA"))
time_levels <- c("1h", "6h", "24h", "48h")
stop_if(!all(as.character(fa$Time) %in% time_levels),
        "发现未登记的 factorial Time 水平。")
fa$Time <- factor(as.character(fa$Time), levels = time_levels)
stop_if(anyNA(fa$Donor) || anyNA(fa$Treatment) || anyNA(fa$Time),
        "factorial 设计因子包含 NA。")

design <- stats::model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = fa)
design_rank <- qr(design)$rank
design_full_rank <- design_rank == ncol(design)
stop_if(!design_full_rank, "factorial 设计矩阵不满秩，停止 QC/PCA。")

# 先进行 count-aware 过滤和 TMM 标准化，再用于样本结构 QC。
dge <- edgeR::DGEList(counts = count[, fa$SampleID, drop = FALSE], samples = fa)
keep <- edgeR::filterByExpr(dge, design = design)
stop_if(sum(keep) < 1000L, "过滤后基因数过少，停止 QC/PCA。")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

library_size <- as.numeric(dge$samples$lib.size)
norm_factor <- as.numeric(dge$samples$norm.factors)
detected_genes <- colSums(dge$counts > 0)
qc <- data.frame(
  SampleID = fa$SampleID,
  Donor = as.character(fa$Donor),
  Treatment = as.character(fa$Treatment),
  Time = as.character(fa$Time),
  library_size = library_size,
  TMM_norm_factor = norm_factor,
  filtered_detected_genes = as.numeric(detected_genes),
  stringsAsFactors = FALSE
)

gene_variance <- apply(log_cpm, 1L, stats::var)
n_top <- min(2000L, nrow(log_cpm))
top_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(n_top)]
pca <- stats::prcomp(t(log_cpm[top_genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_var <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
pca_coordinates <- data.frame(
  SampleID = rownames(pca$x),
  Donor = as.character(fa$Donor),
  Treatment = as.character(fa$Treatment),
  Time = as.character(fa$Time),
  PC1 = pca$x[, 1L],
  PC2 = pca$x[, 2L],
  PC3 = pca$x[, 3L],
  stringsAsFactors = FALSE
)
pca_variance <- data.frame(
  component = paste0("PC", seq_along(pca_var)),
  variance_fraction = pca_var,
  variance_percent = 100 * pca_var,
  stringsAsFactors = FALSE
)

sample_cor <- stats::cor(log_cpm, method = "pearson")
sample_cor_df <- data.frame(SampleID = rownames(sample_cor), sample_cor,
                            check.names = FALSE, stringsAsFactors = FALSE)

qa <- rbind(
  data.frame(check_id = "factorial_samples", status = ifelse(nrow(fa) == 36L, "PASS", "HOLD"),
             detail = nrow(fa), stringsAsFactors = FALSE),
  data.frame(check_id = "design_full_rank", status = ifelse(design_full_rank, "PASS", "HOLD"),
             detail = paste("rank", design_rank, "of", ncol(design)), stringsAsFactors = FALSE),
  data.frame(check_id = "genes_after_filter", status = ifelse(sum(keep) >= 1000L, "PASS", "HOLD"),
             detail = sum(keep), stringsAsFactors = FALSE),
  data.frame(check_id = "finite_log_cpm", status = ifelse(all(is.finite(log_cpm)), "PASS", "HOLD"),
             detail = paste(nrow(log_cpm), "genes x", ncol(log_cpm), "samples"), stringsAsFactors = FALSE),
  data.frame(check_id = "pca_complete", status = ifelse(all(is.finite(pca$x)), "PASS", "HOLD"),
             detail = paste("top variable genes", n_top), stringsAsFactors = FALSE)
)

write_tsv(qc, file.path(out_dir, "sample_qc_after_filter.tsv"))
write_tsv(pca_coordinates, file.path(out_dir, "pca_coordinates.tsv"))
write_tsv(pca_variance, file.path(out_dir, "pca_variance.tsv"))
write_tsv(sample_cor_df, file.path(out_dir, "sample_correlation.tsv"))
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
write_tsv(data.frame(
  input_file = normalizePath(input_file, mustWork = TRUE),
  factorial_samples = nrow(fa), total_genes = nrow(count),
  genes_after_filter = sum(keep), top_variable_genes = n_top,
  design_rank = design_rank, design_columns = ncol(design),
  stringsAsFactors = FALSE
), file.path(out_dir, "run_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
saveRDS(log_cpm, file.path(out_dir, "filtered_tmm_logCPM.rds"))

# QC 图仅用于候选分析审阅，不是正式 Fig. 7 panel。
cond_cols <- c(NC = "#2C7FB8", HCl = "#D95F0E", LA = "#31A354")
png(file.path(out_dir, "pca_condition.png"), width = 1800, height = 1500, res = 220)
plot(pca_coordinates$PC1, pca_coordinates$PC2,
     col = unname(cond_cols[pca_coordinates$Treatment]), pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1L]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2L]),
     main = "Bulk perturbation PCA: treatment")
legend("topright", legend = names(cond_cols), col = cond_cols, pch = 19,
       bty = "n", title = "Treatment")
dev.off()

donor_cols <- c(V1 = "#756BB1", V2 = "#31A354", V3 = "#E6550D")
png(file.path(out_dir, "pca_donor.png"), width = 1800, height = 1500, res = 220)
plot(pca_coordinates$PC1, pca_coordinates$PC2,
     col = unname(donor_cols[pca_coordinates$Donor]), pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1L]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2L]),
     main = "Bulk perturbation PCA: donor")
legend("topright", legend = names(donor_cols), col = donor_cols, pch = 19,
       bty = "n", title = "Donor")
dev.off()

cor_cols <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101L)
png(file.path(out_dir, "sample_correlation_heatmap.png"), width = 1900, height = 1700, res = 220)
op <- par(mar = c(8, 8, 3, 2))
image(seq_len(ncol(sample_cor)), seq_len(nrow(sample_cor)), t(sample_cor[nrow(sample_cor):1L, , drop = FALSE]),
      col = cor_cols, zlim = c(0.7, 1), axes = FALSE, xlab = "", ylab = "", main = "Sample correlation (filtered TMM logCPM)")
axis(1, at = seq_len(ncol(sample_cor)), labels = colnames(sample_cor), las = 2, cex.axis = 0.45)
axis(2, at = seq_len(nrow(sample_cor)), labels = rev(rownames(sample_cor)), las = 2, cex.axis = 0.45)
box()
par(op)
dev.off()

cat("QC_STATUS\t", ifelse(all(qa$status == "PASS"), "PASS", "HOLD"), "\n", sep = "")
cat("GENES_AFTER_FILTER\t", sum(keep), "\n", sep = "")
cat("PC1_PC2_PERCENT\t", sprintf("%.2f,%.2f", 100 * pca_var[1L], 100 * pca_var[2L]), "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
}

# Purpose: Fit donor-adjusted CD8 treatment-time responses and fixed programme enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/02_differential_response.R
stage_cd8_perturbation_02_differential_response <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["count_rdata"]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

# ---- Retained complete source body; CLI/path substitutions only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

input_file <- unname(paths["count_rdata"])
out_dir <- unname(paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
`write_tsv` <- bind_stage_helper(shared_write_tsv_053, environment())
`stop_if` <- bind_stage_helper(shared_stop_if_004, environment())

stop_if(!requireNamespace("edgeR", quietly = TRUE), "R 包 edgeR 未安装。")

e <- new.env(parent = emptyenv())
stop_if(!file.exists(input_file), paste("输入文件不存在:", input_file))
load(input_file, envir = e)
stop_if(!all(c("count", "sample_anno") %in% ls(e, all.names = TRUE)),
        "RData 必须同时包含 count 和 sample_anno。")
count <- get("count", envir = e)
sample_anno <- get("sample_anno", envir = e)
stop_if(!identical(colnames(count), rownames(sample_anno)),
        "count 列名与 sample_anno 行名不一致。")

anno <- sample_anno[colnames(count), c("Donor", "Treatment", "Time"), drop = FALSE]
anno$SampleID <- rownames(anno)
factorial <- as.character(anno$Treatment) %in% c("HCl", "LA", "NC")
fa <- anno[factorial, , drop = FALSE]
fa$Donor <- factor(as.character(fa$Donor))
fa$Treatment <- factor(as.character(fa$Treatment), levels = c("NC", "HCl", "LA"))
time_levels <- c("1h", "6h", "24h", "48h")
stop_if(!all(as.character(fa$Time) %in% time_levels),
        "发现未登记的 factorial Time 水平。")
fa$Time <- factor(as.character(fa$Time), levels = time_levels)
stop_if(anyNA(fa$Donor) || anyNA(fa$Treatment) || anyNA(fa$Time),
        "factorial 设计因子包含 NA。")

design <- stats::model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = fa)
colnames(design) <- make.unique(colnames(design))
design_rank <- qr(design)$rank
stop_if(design_rank != ncol(design), "差异模型设计矩阵不满秩。")

dge <- edgeR::DGEList(counts = count[, fa$SampleID, drop = FALSE], samples = fa)
keep <- edgeR::filterByExpr(dge, design = design)
stop_if(sum(keep) < 1000L, "过滤后基因数过少。")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")
dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
fit <- edgeR::glmQLFit(dge, design, robust = TRUE)

coef_names <- colnames(design)
find_coef <- function(term) {
  hit <- which(coef_names == term)
  if (length(hit) != 1L) {
    hit <- which(coef_names %in% c(term, sub(":", ":", term, fixed = TRUE)))
  }
  stop_if(length(hit) != 1L, paste("找不到唯一模型系数:", term))
  hit
}

contrast_vector <- function(left, right, time) {
  v <- numeric(ncol(design))
  add_mean <- function(treatment, sign) {
    base_term <- paste0("Treatment", treatment)
    if (treatment != "NC") v[find_coef(base_term)] <<- v[find_coef(base_term)] + sign
    if (time != "1h" && treatment != "NC") {
      interaction_terms <- c(
        paste0("Time", time, ":Treatment", treatment),
        paste0("Treatment", treatment, ":Time", time)
      )
      hit <- which(coef_names %in% interaction_terms)
      stop_if(length(hit) != 1L, paste("找不到时间交互项:", paste(interaction_terms, collapse = "/")))
      v[hit] <<- v[hit] + sign
    }
  }
  add_mean(left, 1)
  add_mean(right, -1)
  v
}

contrast_specs <- expand.grid(
  contrast = c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"),
  time = time_levels,
  stringsAsFactors = FALSE
)
contrast_specs$left <- c("LA", "LA", "HCl")[match(contrast_specs$contrast,
                                                    c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))]
contrast_specs$right <- c("NC", "HCl", "NC")[match(contrast_specs$contrast,
                                                       c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))]
contrast_specs$contrast_vector <- lapply(seq_len(nrow(contrast_specs)), function(i) {
  contrast_vector(contrast_specs$left[[i]], contrast_specs$right[[i]], contrast_specs$time[[i]])
})

format_de <- function(test, spec) {
  tab <- as.data.frame(edgeR::topTags(test, n = Inf, sort.by = "none")$table,
                       stringsAsFactors = FALSE)
  tab$gene_id <- rownames(tab)
  tab <- tab[, c("gene_id", setdiff(colnames(tab), "gene_id")), drop = FALSE]
  tab$contrast <- spec$contrast
  tab$time <- spec$time
  tab$left <- spec$left
  tab$right <- spec$right
  tab
}

contrast_manifest <- do.call(rbind, lapply(seq_len(nrow(contrast_specs)), function(i) {
  v <- contrast_specs$contrast_vector[[i]]
  data.frame(
    contrast_id = paste(contrast_specs$contrast[[i]], contrast_specs$time[[i]], sep = "__"),
    contrast = contrast_specs$contrast[[i]], time = contrast_specs$time[[i]],
    left = contrast_specs$left[[i]], right = contrast_specs$right[[i]],
    n_left = sum(fa$Treatment == contrast_specs$left[[i]] & fa$Time == contrast_specs$time[[i]]),
    n_right = sum(fa$Treatment == contrast_specs$right[[i]] & fa$Time == contrast_specs$time[[i]]),
    contrast_nonzero = paste(paste(coef_names[which(v != 0)], v[v != 0], sep = "="), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

for (i in seq_len(nrow(contrast_specs))) {
  spec <- contrast_specs[i, , drop = FALSE]
  test <- edgeR::glmQLFTest(fit, contrast = contrast_specs$contrast_vector[[i]])
  result <- format_de(test, spec)
  out_name <- paste0("de_", spec$contrast, "_", spec$time, ".tsv")
  write_tsv(result, file.path(out_dir, out_name))
}

interaction_terms <- grep("^Time.*:Treatment|^Treatment.*:Time", coef_names)
stop_if(length(interaction_terms) != 6L, "时间×处理交互项数量不是预期的 6。")
interaction_test <- edgeR::glmQLFTest(fit, coef = interaction_terms)
interaction_result <- as.data.frame(edgeR::topTags(interaction_test, n = Inf, sort.by = "none")$table,
                                    stringsAsFactors = FALSE)
interaction_result$gene_id <- rownames(interaction_result)
interaction_result <- interaction_result[, c("gene_id", setdiff(colnames(interaction_result), "gene_id")), drop = FALSE]
write_tsv(interaction_result, file.path(out_dir, "time_treatment_interaction.tsv"))

sample_design <- data.frame(fa, stringsAsFactors = FALSE)
write_tsv(sample_design, file.path(out_dir, "factorial_sample_design.tsv"))
write_tsv(data.frame(coefficient = coef_names, index = seq_along(coef_names), stringsAsFactors = FALSE),
          file.path(out_dir, "model_coefficients.tsv"))
write_tsv(contrast_manifest, file.path(out_dir, "contrast_manifest.tsv"))
write_tsv(data.frame(
  input_file = normalizePath(input_file, mustWork = TRUE),
  factorial_samples = nrow(fa), genes_before_filter = nrow(count),
  genes_after_filter = sum(keep), design_rank = design_rank,
  design_columns = ncol(design), n_interaction_terms = length(interaction_terms),
  stringsAsFactors = FALSE
), file.path(out_dir, "run_summary.tsv"))
write_tsv(data.frame(
  check_id = c("factorial_samples", "design_full_rank", "genes_after_filter",
               "all_time_contrasts", "interaction_test"),
  status = c(
    ifelse(nrow(fa) == 36L, "PASS", "HOLD"),
    ifelse(design_rank == ncol(design), "PASS", "HOLD"),
    ifelse(sum(keep) >= 1000L, "PASS", "HOLD"),
    ifelse(nrow(contrast_specs) == 12L, "PASS", "HOLD"),
    ifelse(length(interaction_terms) == 6L, "PASS", "HOLD")
  ),
  detail = c(nrow(fa), paste(design_rank, "of", ncol(design)), sum(keep),
             nrow(contrast_specs), length(interaction_terms)),
  stringsAsFactors = FALSE
), file.path(out_dir, "qa_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

cat("DE_STATUS\tPASS\n")
cat("GENES_AFTER_FILTER\t", sum(keep), "\n", sep = "")
cat("CONTRAST_FILES\t", nrow(contrast_specs), "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
}

# Purpose: Fit donor-adjusted CD8 treatment-time responses and fixed programme enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/03_fixed_programme_gsea.R
stage_cd8_perturbation_03_fixed_programme_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("count_rdata","lm12_membership","melanoma_counts","melanoma_deg")]
gmt_basenames <- c("h.all.v2025.1.Hs.symbols.gmt","c2.cp.reactome.v2025.1.Hs.symbols.gmt","c5.go.bp.v2024.1.Hs.symbols.gmt","c7.all.v2025.1.Hs.symbols.gmt")
de_contrasts <- c("LA_minus_NC")
de_basenames <- as.vector(outer(de_contrasts, c("1h","6h","24h","48h"), function(x,y) paste0("de_",x,"_",y,".tsv")))
required_files <- c(unname(input_paths), file.path(unname(paths["msigdb_dir"]),gmt_basenames), file.path(unname(paths["bulk_de_dir"]),de_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("msigdb_dir","bulk_de_dir")])
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

cfg_paths <- paths

# ---- Retained complete source body; documented path adaptations only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(GSVA)
  library(fgsea)
})

# v0.2 专门补足：基因集 ID 映射、donor 一致性/区间、leading-edge 拆解和 GSEA 可复现性。
analysis_seed <- 20260911L
boundary_seed <- 20260912L
set.seed(analysis_seed)

out_dir <- unname(cfg_paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  bulk_input = unname(cfg_paths["count_rdata"]),
  lm12_long = unname(cfg_paths["lm12_membership"]),
  msigdb_root = unname(cfg_paths["msigdb_dir"]),
  melanoma_counts = unname(cfg_paths["melanoma_counts"]),
  melanoma_deg = unname(cfg_paths["melanoma_deg"]),
  bulk_de_root = unname(cfg_paths["bulk_de_dir"])
)

msigdb_files <- c(
  hallmark = file.path(paths$msigdb_root, "h.all.v2025.1.Hs.symbols.gmt"),
  reactome = file.path(paths$msigdb_root, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"),
  go_bp = file.path(paths$msigdb_root, "c5.go.bp.v2024.1.Hs.symbols.gmt"),
  c7 = file.path(paths$msigdb_root, "c7.all.v2025.1.Hs.symbols.gmt")
)

fixed_msigdb <- data.table(
  feature_id = c(
    "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
    "REACTOME_TCR_SIGNALING",
    "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN",
    "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
  ),
  feature_label = c(
    "CD8 activation",
    "TCR signaling",
    "T-cell mediated cytotoxicity",
    "IFN-gamma response",
    "Memory-enriched (Goldrath effector vs memory, DN)",
    "Effector-enriched (Goldrath effector vs memory, UP)",
    "Glycolysis",
    "Hypoxia",
    "Unfolded protein response"
  ),
  source_collection = c("C5_GO_BP", "C2_REACTOME", "C5_GO_BP", "HALLMARK", "C7_IMMUNO", "C7_IMMUNO", "HALLMARK", "HALLMARK", "HALLMARK"),
  source_file_key = c("go_bp", "reactome", "go_bp", "hallmark", "c7", "c7", "hallmark", "hallmark", "hallmark"),
  direction = c("activity", "activity", "activity", "activity", "memory_enriched", "effector_enriched", "activity", "activity", "activity")
)
`write_tsv` <- bind_stage_helper(shared_write_tsv_046, environment())
`stop_if` <- bind_stage_helper(shared_stop_if_003, environment())

for (f in c(unlist(paths[c("bulk_input", "lm12_long", "melanoma_counts", "melanoma_deg", "bulk_de_root")]), unname(msigdb_files))) {
  stop_if(!file.exists(f), paste("missing input:", f))
}

input_paths <- c(unlist(paths[c("bulk_input", "lm12_long", "melanoma_counts", "melanoma_deg")], use.names = TRUE), msigdb_files)
file_manifest <- data.table(
  role = c(c("bulk_input", "lm12_long", "melanoma_counts", "melanoma_deg"), paste0("msigdb_", names(msigdb_files))),
  path = unname(input_paths)
)
file_manifest[, `:=`(exists = file.exists(path), bytes = as.numeric(file.info(path)$size), sha256 = unname(tools::sha256sum(path)))]
write_tsv(file_manifest, file.path(out_dir, "input_file_manifest.tsv"))

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  toks <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(toks, function(x) x[[1L]], character(1))
  sets <- lapply(toks, function(x) unique(x[-seq_len(min(2L, length(x)))]))
  names(sets) <- ids
  sets <- sets[!duplicated(names(sets))]
  list(sets = sets, n_source = lengths(sets))
}

gmt_cache <- lapply(msigdb_files, read_gmt)
names(gmt_cache) <- names(msigdb_files)
msigdb_sets <- lapply(seq_len(nrow(fixed_msigdb)), function(i) {
  key <- fixed_msigdb$source_file_key[i]
  id <- fixed_msigdb$feature_id[i]
  stop_if(!id %in% names(gmt_cache[[key]]$sets), paste("fixed MSigDB set missing:", id, key))
  gmt_cache[[key]]$sets[[id]]
})
names(msigdb_sets) <- fixed_msigdb$feature_id

lm12_tab <- fread(paths$lm12_long, sep = "\t", quote = "\"", data.table = TRUE)
stop_if(!all(c("feature_id", "gene", "pathway") %in% names(lm12_tab)), "LM12 membership columns incomplete")
lm12_ids <- unique(as.character(lm12_tab$feature_id))
lm12_sets <- setNames(lapply(lm12_ids, function(id) unique(as.character(lm12_tab[feature_id == id, gene]))), lm12_ids)
lm12_meta <- data.table(feature_id = lm12_ids, feature_label = vapply(lm12_ids, function(id) unique(as.character(lm12_tab[feature_id == id, pathway]))[1L], character(1)))
lm12_meta[, `:=`(feature_layer = "LM12_FIG2_FROZEN", source_collection = "FIG2_FROZEN_LM12", source_file_key = "lm12_long", direction = "activity")]

msigdb_meta <- copy(fixed_msigdb)
msigdb_meta[, `:=`(feature_layer = "MSIGDB_FULL_SET", source_path = unname(msigdb_files[source_file_key]))]
msigdb_meta[, n_genes_source := lengths(msigdb_sets[feature_id])]
msigdb_meta[, set_sha256 := vapply(msigdb_sets[feature_id], function(g) digest::digest(sort(unique(g)), algo = "sha256"), character(1))]

all_sets <- c(setNames(lm12_sets, paste0("LM12::", names(lm12_sets))), setNames(msigdb_sets, paste0("MSIGDB::", names(msigdb_sets))))
feature_meta <- rbind(
  lm12_meta[, .(feature_id = paste0("LM12::", feature_id), feature_label, feature_layer, source_collection, direction)],
  msigdb_meta[, .(feature_id = paste0("MSIGDB::", feature_id), feature_label, feature_layer, source_collection, direction)]
)
feature_meta[, feature_id_raw := sub("^[^:]+::", "", feature_id)]
set_label <- setNames(feature_meta$feature_label, feature_meta$feature_id)
set_layer <- setNames(feature_meta$feature_layer, feature_meta$feature_id)
source_size <- data.table(feature_id = names(all_sets), n_genes_source = lengths(all_sets))
source_path <- c(setNames(rep(paths$lm12_long, length(lm12_sets)), names(all_sets)[seq_along(lm12_sets)]),
                 setNames(unname(msigdb_files[msigdb_meta$source_file_key]), names(all_sets)[(length(lm12_sets) + 1L):length(all_sets)]))

write_tsv(rbind(
  lm12_meta[, .(feature_id, feature_label, feature_layer, source_collection, source_path = paths$lm12_long, direction, n_genes_source = unname(lengths(lm12_sets)[feature_id]), set_sha256 = vapply(lm12_sets[feature_id], function(g) digest::digest(sort(unique(g)), algo = "sha256"), character(1)))],
  msigdb_meta[, .(feature_id, feature_label, feature_layer, source_collection, source_path, direction, n_genes_source, set_sha256)]
), file.path(out_dir, "fixed_gene_set_manifest.tsv"))

load_count <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  stop_if(!all(c("count", "sample_anno") %in% ls(e, all.names = TRUE)), "bulk RData lacks count/sample_anno")
  list(count = get("count", e), sample_anno = get("sample_anno", e))
}

score_ssgsea <- function(logcpm, sets, dataset) {
  present <- lapply(sets, function(g) intersect(unique(g), rownames(logcpm)))
  names(present) <- names(sets)
  status <- data.table(dataset = dataset, feature_id = names(present), n_genes_input = lengths(sets), n_genes_present = lengths(present), coverage = lengths(present) / lengths(sets), status = ifelse(lengths(present) >= 3L, "PASS", "HOLD"))
  stop_if(any(status$status == "HOLD"), paste("set coverage below 3 genes:", paste(status[status == "HOLD", feature_id], collapse = ",")))
  set.seed(analysis_seed)
  score <- GSVA::gsva(GSVA::ssgseaParam(exprData = logcpm, geneSets = present, minSize = 1L, alpha = 0.25, normalize = TRUE), verbose = FALSE)
  score <- as.matrix(score)
  returned_ids <- rownames(score)
  stop_if(is.null(returned_ids) || !setequal(returned_ids, names(present)), paste(dataset, "GSVA returned IDs do not match fixed gene-set IDs"))
  returned_index <- match(names(present), returned_ids)
  stop_if(anyNA(returned_index), paste(dataset, "GSVA ID matching failed"))
  score <- score[returned_index, , drop = FALSE]
  rownames(score) <- names(present)
  mapping <- data.table(dataset = dataset, expected_id = names(present), returned_id = returned_ids[returned_index], returned_index = returned_index, mapping_ok = names(present) == returned_ids[returned_index])
  stop_if(!all(mapping$mapping_ok), paste(dataset, "GSVA ID mapping is not exact"))
  list(score = score, present = present, coverage = status, mapping = mapping)
}

make_long <- function(score, sample_sheet, score_layer) {
  long <- melt(as.data.table(t(score), keep.rownames = "sample_id"), id.vars = "sample_id", variable.name = "feature_id", value.name = "ssgsea_score")
  long[, feature_id := as.character(feature_id)]
  long <- merge(long, as.data.table(sample_sheet), by = "sample_id", all.x = TRUE, sort = FALSE)
  long[, `:=`(score_layer = score_layer, feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]))]
  stop_if(anyNA(long$feature_label) || anyNA(long$feature_layer), paste(score_layer, "feature label mapping failed"))
  setcolorder(long, c("sample_id", "feature_id", "feature_label", "feature_layer", "score_layer"))
  long[]
}

model_score_contrasts <- function(score, design, contrast_specs, feature_meta_local) {
  fit <- limma::lmFit(score, design)
  coef_names <- colnames(design)
  find_coef <- function(term) {
    hit <- which(coef_names == term)
    stop_if(length(hit) != 1L, paste("coefficient not unique:", term))
    hit
  }
  out <- lapply(seq_len(nrow(contrast_specs)), function(i) {
    cf <- limma::contrasts.fit(fit, contrast_specs$vector[[i]])
    cf <- limma::eBayes(cf, robust = TRUE)
    tt <- as.data.table(limma::topTable(cf, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature_id")
    stop_if(!setequal(tt$feature_id, rownames(score)), "score contrast lost feature IDs")
    tt <- tt[match(rownames(score), feature_id)]
    se <- as.numeric(cf$stdev.unscaled[, 1L] * cf$sigma)
    names(se) <- rownames(cf$stdev.unscaled)
    se <- se[match(tt$feature_id, names(se))]
    df_resid <- max(cf$df.residual, na.rm = TRUE)
    tt[, `:=`(contrast = contrast_specs$contrast[i], time = contrast_specs$time[i], SE = se, df = df_resid)]
    tt[, `:=`(CI_low = logFC - qt(0.975, df_resid) * SE, CI_high = logFC + qt(0.975, df_resid) * SE)]
    tt[, .(feature_id, contrast, time, logFC, AveExpr, t, P.Value, adj.P.Val, B, SE, df, CI_low, CI_high)]
  })
  merge(rbindlist(out, fill = TRUE), feature_meta_local, by = "feature_id", all.x = TRUE, sort = FALSE)
}

bulk_de_stats <- function(time) {
  f <- file.path(paths$bulk_de_root, paste0("de_LA_minus_NC_", time, ".tsv"))
  stop_if(!file.exists(f), paste("missing existing bulk DE:", f))
  d <- fread(f, sep = "\t")
  d[, rank_score := sign(logFC) * -log10(pmax(PValue, .Machine$double.xmin))]
  d[!duplicated(gene_id), .(gene = gene_id, rank_score, logFC, PValue, FDR)]
}

run_fgsea <- function(stats, sets, context, time = NA_character_, donor = NA_character_, seed = analysis_seed) {
  stats <- stats[is.finite(stats)]
  stats <- stats[!duplicated(names(stats))]
  set.seed(seed)
  res <- fgsea::fgseaMultilevel(pathways = sets, stats = stats, minSize = 3L, maxSize = 5000L, eps = 0)
  res <- as.data.table(res)
  res[, pathway := as.character(pathway)]
  res[, `:=`(context = context, time = time, donor = donor, feature_id = pathway, seed = seed)]
  res[, pathway := NULL]
  res[]
}

build_bulk <- function() {
  obj <- load_count(paths$bulk_input)
  count <- obj$count
  ann <- as.data.table(obj$sample_anno, keep.rownames = "sample_id")
  stop_if(!identical(colnames(count), ann$sample_id), "bulk count/sample annotation order mismatch")
  ann <- ann[Treatment %in% c("NC", "HCl", "LA")]
  ann[, `:=`(Donor = factor(Donor), Treatment = factor(Treatment, levels = c("NC", "HCl", "LA")), Time = factor(Time, levels = c("1h", "6h", "24h", "48h")))]
  stop_if(anyNA(ann$Donor) || anyNA(ann$Treatment) || anyNA(ann$Time), "bulk factorial metadata contains NA")
  design <- model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = ann)
  stop_if(qr(design)$rank != ncol(design), "bulk design not full rank")
  dge <- edgeR::DGEList(counts = count[, ann$sample_id, drop = FALSE], samples = ann)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  scored <- score_ssgsea(logcpm, all_sets, "bulk")
  sample_sheet <- as.data.frame(ann)
  score_long <- make_long(scored$score, sample_sheet, "bulk_ssGSEA")
  coef_names <- colnames(design)
  find_coef <- function(term) {
    hit <- which(coef_names == term)
    stop_if(length(hit) != 1L, paste("coefficient not unique:", term))
    hit
  }
  contrast_vector <- function(left, right, time) {
    v <- numeric(ncol(design))
    add <- function(treatment, sign) {
      if (treatment != "NC") {
        v[find_coef(paste0("Treatment", treatment))] <<- v[find_coef(paste0("Treatment", treatment))] + sign
        if (time != "1h") {
          hit <- which(coef_names %in% c(paste0("Time", time, ":Treatment", treatment), paste0("Treatment", treatment, ":Time", time)))
          stop_if(length(hit) != 1L, "bulk interaction coefficient missing")
          v[hit] <<- v[hit] + sign
        }
      }
    }
    add(left, 1); add(right, -1); v
  }
  specs <- CJ(contrast = c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"), time = c("1h", "6h", "24h", "48h"), unique = TRUE)
  specs[, `:=`(left = c("LA", "LA", "HCl")[match(contrast, c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))], right = c("NC", "HCl", "NC")[match(contrast, c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))])]
  specs[, vector := lapply(seq_len(.N), function(i) contrast_vector(left[i], right[i], time[i]))]
  contrasts <- model_score_contrasts(scored$score, design, specs, feature_meta)
  donor_rows <- list()
  for (i in seq_len(nrow(specs))) {
    sp <- specs[i]
    sub <- score_long[Time == sp$time & Treatment %in% c(sp$left, sp$right)]
    left <- sub[Treatment == sp$left, .(Donor = as.character(Donor), feature_id, left_score = ssgsea_score)]
    right <- sub[Treatment == sp$right, .(Donor = as.character(Donor), feature_id, right_score = ssgsea_score)]
    dd <- merge(left, right, by = c("Donor", "feature_id"), all = FALSE)
    dd[, `:=`(contrast = sp$contrast, time = sp$time, donor_delta = left_score - right_score)]
    donor_rows[[length(donor_rows) + 1L]] <- dd[, .(contrast, time, feature_id, donor = Donor, donor_delta)]
  }
  donor_delta <- rbindlist(donor_rows, fill = TRUE)
  donor_summary <- donor_delta[, .(n_donors = uniqueN(donor), mean_delta = mean(donor_delta), median_delta = median(donor_delta), sd_delta = sd(donor_delta), n_positive = sum(donor_delta > 0), n_negative = sum(donor_delta < 0), n_zero = sum(donor_delta == 0)), by = .(contrast, time, feature_id)]
  donor_summary[, se_delta := sd_delta / sqrt(n_donors)]
  donor_summary[, `:=`(CI_low = mean_delta - qt(0.975, pmax(n_donors - 1, 1)) * se_delta, CI_high = mean_delta + qt(0.975, pmax(n_donors - 1, 1)) * se_delta, donor_all_same_direction = n_positive == n_donors | n_negative == n_donors)]
  list(score_long = score_long, score_contrasts = contrasts, donor_delta = donor_delta, donor_summary = donor_summary, sample_sheet = sample_sheet, ann = ann, logcpm = logcpm, present = scored$present, coverage = scored$coverage, mapping = scored$mapping, design = design, design_rank = qr(design)$rank)
}

build_melanoma <- function() {
  obj <- readRDS(paths$melanoma_counts)
  stop_if(!all(c("counts", "sample_sheet") %in% names(obj)), "melanoma CD8 pseudobulk bundle incomplete")
  count <- obj$counts; ss_all <- as.data.table(obj$sample_sheet)
  req <- c("unit_id", "sample_ID_long", "patient", "sample_C_subtype", "pre_post", "eligible_ge20", "n_cells", "library_size")
  stop_if(!all(req %in% names(ss_all)), "melanoma sample sheet incomplete")
  ss_all <- ss_all[eligible_ge20 == TRUE & sample_C_subtype %in% c("C1", "C4")]
  stop_if(anyDuplicated(ss_all$unit_id) > 0L, "duplicate melanoma CD8 unit IDs")
  ss <- copy(ss_all); selection_mode <- "all_patients_unique"
  if (anyDuplicated(ss$patient) > 0L) {
    ss[, pre_priority := fifelse(pre_post == "Pre", 0L, 1L)]
    setorder(ss, patient, pre_priority, -n_cells, -library_size, sample_ID_long)
    ss <- ss[, .SD[1L], by = patient]; ss[, pre_priority := NULL]
    selection_mode <- "independent_patient_fallback_pre_priority"
  }
  stop_if(anyDuplicated(ss$patient) > 0L, "duplicate melanoma patients after selection")
  selection_manifest <- copy(ss_all); selection_manifest[, retained := unit_id %in% ss$unit_id]; selection_manifest[, selection_reason := fifelse(retained, selection_mode, "dropped_duplicate_patient")]
  count <- count[, ss$unit_id, drop = FALSE]
  group <- factor(ss$sample_C_subtype, levels = c("C4", "C1")); ss[, group := group]; ss[, pre_post := factor(pre_post)]
  design <- model.matrix(~ 0 + group + pre_post, data = ss)
  stop_if(qr(design)$rank != ncol(design), "melanoma CD8 design not full rank")
  dge <- edgeR::DGEList(counts = count, samples = ss); dge <- edgeR::calcNormFactors(dge, method = "TMM"); logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  scored <- score_ssgsea(logcpm, all_sets, "melanoma_CD8")
  sample_sheet <- data.frame(sample_id = ss$unit_id, sample_ID_long = ss$sample_ID_long, patient = ss$patient, state = ss$sample_C_subtype, pre_post = as.character(ss$pre_post), cell_type = "CD8 T", n_cells = ss$n_cells, stringsAsFactors = FALSE)
  score_long <- make_long(scored$score, sample_sheet, "melanoma_CD8_ssGSEA")
  fit <- limma::lmFit(scored$score, design); cf <- limma::contrasts.fit(fit, c(-1, 1, rep(0, ncol(design) - 2L))); cf <- limma::eBayes(cf, robust = TRUE)
  tt <- as.data.table(limma::topTable(cf, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature_id")
  stop_if(!setequal(tt$feature_id, rownames(scored$score)), "melanoma score contrast lost feature IDs")
  tt <- tt[match(rownames(scored$score), feature_id)]
  se <- as.numeric(cf$stdev.unscaled[, 1L] * cf$sigma); names(se) <- rownames(cf$stdev.unscaled); se <- se[match(tt$feature_id, names(se))]
  df_resid <- max(cf$df.residual, na.rm = TRUE)
  tt[, `:=`(contrast = "C1_minus_C4", SE = se, df = df_resid)]
  tt[, `:=`(CI_low = logFC - qt(0.975, df_resid) * SE, CI_high = logFC + qt(0.975, df_resid) * SE)]
  tt <- tt[, .(feature_id, contrast, logFC, AveExpr, t, P.Value, adj.P.Val, B, SE, df, CI_low, CI_high)]
  tt <- merge(tt, feature_meta, by = "feature_id", all.x = TRUE, sort = FALSE)
  group_summary <- score_long[, .(mean_score = mean(ssgsea_score), sd_score = sd(ssgsea_score), n_samples = .N, n_patients = uniqueN(patient)), by = .(feature_id, feature_label, feature_layer, state)]
  deg <- fread(cmd = paste("gzip -cd", shQuote(paths$melanoma_deg)), sep = "\t", quote = "\"", data.table = TRUE); deg <- deg[cell_type == "CD8 T" & comparison == "C1_vs_C4"]
  union_genes <- unique(unlist(all_sets)); deg_fixed <- deg[gene %in% union_genes]
  list(score_long = score_long, score_contrasts = tt, group_summary = group_summary, sample_sheet = sample_sheet, present = scored$present, coverage = scored$coverage, mapping = scored$mapping, design = design, design_rank = qr(design)$rank, n_samples = nrow(ss), n_patients = uniqueN(ss$patient), deg_fixed = deg_fixed, selection_manifest = selection_manifest, selection_mode = selection_mode, logcpm = logcpm)
}

if (FALSE) {
# v0.1 tail retained below only as an auditable historical code path; it is not executed.
bulk <- build_bulk(); mel <- build_melanoma()
membership <- merge(feature_meta, merge(bulk$coverage[, .(feature_id, n_genes_bulk = n_genes_present, bulk_coverage = coverage)], mel$coverage[, .(feature_id, n_genes_melanoma = n_genes_present, melanoma_coverage = coverage)], by = "feature_id", all = TRUE), by = "feature_id", all = TRUE, sort = FALSE)
membership[, `:=`(n_genes_source = ifelse(grepl("^LM12::", feature_id), lengths(lm12_sets)[sub("^LM12::", "", feature_id)], msigdb_meta$n_genes_source[match(sub("^MSIGDB::", "", feature_id), msigdb_meta$feature_id)]), source_path = ifelse(grepl("^LM12::", feature_id), paths$lm12_long, unname(msigdb_files[msigdb_meta$source_file_key[match(sub("^MSIGDB::", "", feature_id), msigdb_meta$feature_id)] ]))) ]
membership[, `:=`(bulk_score_status = ifelse(n_genes_bulk >= 3L, "PASS", "HOLD"), melanoma_score_status = ifelse(n_genes_melanoma >= 3L, "PASS", "HOLD"))]

write_tsv(membership, file.path(out_dir, "fixed_gene_set_coverage_manifest.tsv"))
write_tsv(rbind(bulk$mapping, mel$mapping, fill = TRUE), file.path(out_dir, "score_id_mapping_qa.tsv"))
write_tsv(bulk$score_long, file.path(out_dir, "bulk_fixed_readout_ssgsea_long.tsv"))
write_tsv(bulk$score_contrasts, file.path(out_dir, "bulk_fixed_readout_ssgsea_contrasts.tsv"))
write_tsv(bulk$donor_delta, file.path(out_dir, "bulk_fixed_readout_donor_deltas.tsv"))
write_tsv(bulk$donor_summary, file.path(out_dir, "bulk_fixed_readout_donor_summary.tsv"))
write_tsv(mel$score_long, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_long.tsv"))
write_tsv(mel$score_contrasts, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_c1_vs_c4.tsv"))
write_tsv(mel$group_summary, file.path(out_dir, "melanoma_cd8_fixed_readout_group_summary.tsv"))
write_tsv(mel$selection_manifest, file.path(out_dir, "melanoma_cd8_patient_selection_manifest.tsv"))
write_tsv(mel$deg_fixed, file.path(out_dir, "melanoma_cd8_fixed_readout_gene_deg.tsv"))
write_tsv(bridge, file.path(out_dir, "melanoma_cd8_bulk_bridge_by_time.tsv"))
write_tsv(bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_bridge_summary.tsv"))
write_tsv(gene_bridge, file.path(out_dir, "melanoma_cd8_bulk_gene_bridge.tsv"))
write_tsv(flatten_gsea(bulk_fgsea_raw), file.path(out_dir, "bulk_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(mel_fgsea_raw), file.path(out_dir, "melanoma_cd8_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(bulk_donor_gsea_raw), file.path(out_dir, "bulk_donor_directional_preranked_gsea.tsv"))
write_tsv(gsea_bridge, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_by_time.tsv"))
write_tsv(gsea_bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_summary.tsv"))
write_tsv(gsea_repro, file.path(out_dir, "gsea_reproducibility_boundary.tsv"))
write_tsv(method_direction_audit, file.path(out_dir, "ssgsea_gsea_direction_audit.tsv"))
write_tsv(memory_direction_audit, file.path(out_dir, "memory_direction_audit.tsv"))
write_tsv(diag_specs, file.path(out_dir, "diagnostic_gene_set_manifest.tsv"))
write_tsv(edge_long, file.path(out_dir, "gsea_leading_edge_long.tsv"))
write_tsv(edge_cat, file.path(out_dir, "gsea_leading_edge_category_long.tsv"))
write_tsv(edge_summary, file.path(out_dir, "gsea_leading_edge_category_summary.tsv"))
write_tsv(gsea_repro[boundary_a | boundary_b | boundary_changed], file.path(out_dir, "gsea_boundary_focus.tsv"))
write_tsv(data.table(
  object = c("bulk_count", "melanoma_cd8_pseudobulk"),
  analysis_unit = c("donor x condition x time sample", "sample x native CD8 T pseudobulk; patient biological unit"),
  primary_comparison = c("time-matched LA_minus_NC", "C1_minus_C4"),
  n_independent_units = c(length(unique(bulk$sample_sheet$Donor)), mel$n_patients),
  interpretation = c("CD8 bulk aggregate response", "melanoma CD8 state association")
), file.path(out_dir, "analysis_unit_manifest.tsv"))

mapping_ok <- all(bulk$mapping$mapping_ok) && all(mel$mapping$mapping_ok)
gsea_rows_ok <- nrow(bulk_fgsea_raw) == length(all_sets) * 4L && nrow(mel_fgsea_raw) == length(all_sets)
donor_rows_ok <- nrow(bulk$donor_delta[contrast == "LA_minus_NC"]) == length(all_sets) * 4L * 3L
repro_ok <- nrow(gsea_repro) > 0L && all(is.finite(gsea_repro$abs_NES_delta)) && all(is.finite(gsea_repro$abs_FDR_delta))
qa <- data.table(
  check_id = c("bulk_factorial_samples", "bulk_design_rank", "melanoma_c1c4_samples", "melanoma_c1c4_patients", "fixed_readout_count", "score_id_mapping_exact", "bulk_coverage", "melanoma_coverage", "all_ssgsea_finite", "bulk_primary_donor_rows", "gsea_rows_all_fixed", "donor_gsea_rows", "gsea_seed_reproducible_run", "leading_edge_table"),
  status = c(
    ifelse(nrow(bulk$sample_sheet) == 36L, "PASS", "HOLD"),
    ifelse(bulk$design_rank == ncol(bulk$design), "PASS", "HOLD"),
    ifelse(mel$n_samples >= 6L, "PASS", "HOLD"),
    ifelse(mel$n_patients >= 6L, "PASS", "HOLD"),
    ifelse(length(all_sets) == 21L, "PASS", "HOLD"),
    ifelse(mapping_ok, "PASS", "HOLD"),
    ifelse(all(membership$bulk_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(membership$melanoma_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(is.finite(bulk$score_long$ssgsea_score)) && all(is.finite(mel$score_long$ssgsea_score)), "PASS", "HOLD"),
    ifelse(donor_rows_ok, "PASS", "HOLD"),
    ifelse(gsea_rows_ok, "PASS", "HOLD"),
    ifelse(nrow(bulk_donor_gsea_raw) == length(all_sets) * 4L * 3L, "PASS", "HOLD"),
    ifelse(repro_ok, "PASS", "HOLD"),
    ifelse(nrow(edge_long) > 0L, "PASS", "HOLD")
  ),
  detail = c(nrow(bulk$sample_sheet), paste(bulk$design_rank, "of", ncol(bulk$design)), mel$n_samples, mel$n_patients, length(all_sets), paste(all(bulk$mapping$mapping_ok), all(mel$mapping$mapping_ok), sep = ";"), sprintf("min %.3f", min(membership$bulk_coverage)), sprintf("min %.3f", min(membership$melanoma_coverage)), paste(nrow(bulk$score_long), nrow(mel$score_long), sep = " + "), nrow(bulk$donor_delta[contrast == "LA_minus_NC"]), paste(nrow(bulk_fgsea_raw), nrow(mel_fgsea_raw), sep = " + "), nrow(bulk_donor_gsea_raw), sprintf("max |delta NES| %.4g; max |delta FDR| %.4g", max(gsea_repro$abs_NES_delta), max(gsea_repro$abs_FDR_delta)), nrow(edge_long))
)
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
overall_status <- if (all(qa$status == "PASS")) "PASS_CANDIDATE" else "HOLD_CANDIDATE"
writeLines(c(paste0("status=", overall_status), paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)), paste0("analysis_seed=", analysis_seed), paste0("boundary_seed=", boundary_seed), paste0("fixed_readouts=", length(all_sets)), paste0("bulk_samples=", nrow(bulk$sample_sheet), ";bulk_donors=", length(unique(bulk$sample_sheet$Donor))), paste0("melanoma_cd8_samples=", mel$n_samples, ";melanoma_patients=", mel$n_patients, ";selection=", mel$selection_mode), paste0("gsea_rows=", nrow(bulk_fgsea_raw), "+", nrow(mel_fgsea_raw)), paste0("leading_edge_rows=", nrow(edge_long))), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("STATUS\t", overall_status, "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")

bulk_primary <- bulk$score_contrasts[contrast == "LA_minus_NC"]
mel_primary <- mel$score_contrasts
bridge <- merge(mel_primary[, .(feature_id, melanoma_logFC = logFC, melanoma_p = P.Value, melanoma_fdr = adj.P.Val)], bulk_primary[, .(feature_id, time, bulk_logFC = logFC, bulk_p = P.Value, bulk_fdr = adj.P.Val)], by = "feature_id", allow.cartesian = TRUE)
donor_support <- merge(bulk$donor_delta[contrast == "LA_minus_NC"], mel_primary[, .(feature_id, melanoma_logFC = logFC)], by = "feature_id", all.x = TRUE)
donor_support[, donor_same_direction := sign(donor_delta) == sign(melanoma_logFC)]
support_summary <- donor_support[, .(n_donors = uniqueN(donor), n_donors_same_direction = sum(donor_same_direction, na.rm = TRUE), donor_same_fraction = mean(donor_same_direction, na.rm = TRUE)), by = .(feature_id, time)]
bridge <- merge(bridge, support_summary, by = c("feature_id", "time"), all.x = TRUE)
bridge[, `:=`(feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]), melanoma_direction = fifelse(melanoma_logFC > 0, "positive", fifelse(melanoma_logFC < 0, "negative", "zero")), bulk_direction = fifelse(bulk_logFC > 0, "positive", fifelse(bulk_logFC < 0, "negative", "zero")))]
bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
bridge[, evidence_class := ifelse(melanoma_fdr < 0.05 & bulk_fdr < 0.05 & donor_same_fraction >= 2/3, ifelse(direction_relation == "same_direction", "supported_concordance", "supported_divergence"), ifelse(direction_relation == "same_direction" & (donor_same_fraction >= 1 | melanoma_p < 0.05 | bulk_p < 0.05), "suggestive_concordance", "insufficient_or_mixed"))]
bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(bridge, layer_order, feature_id, time)
bridge[, layer_order := NULL]
bridge_summary <- bridge[, .(bridge_class = if (all(evidence_class == "supported_concordance")) "supported_concordance" else if (any(evidence_class == "supported_divergence")) "supported_divergence_or_mixed" else if (any(evidence_class == "suggestive_concordance")) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), mean_abs_bulk_logFC = mean(abs(bulk_logFC), na.rm = TRUE), melanoma_logFC = first(melanoma_logFC), melanoma_fdr = first(melanoma_fdr)), by = .(feature_id, feature_label, feature_layer)]

gene_bridge <- merge(mel$deg_fixed[, .(gene, melanoma_logFC = log2fc_c1_minus_ref, melanoma_p = p_value, melanoma_fdr = fdr)], rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) { d <- bulk_de_stats(tm); d[, time := tm]; d[gene %in% unique(unlist(all_sets))] })), by = "gene", allow.cartesian = TRUE)
gene_bridge[, direction_relation := fifelse(sign(melanoma_logFC) == sign(logFC), "same_direction", "opposite_direction")]
gene_bridge[, set_membership := vapply(gene, function(g) paste(names(all_sets)[vapply(all_sets, function(s) g %in% s, logical(1))], collapse = ";"), character(1))]

bulk_fgsea <- rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) {
  d <- bulk_de_stats(tm); stats <- d$rank_score; names(stats) <- d$gene
  run_fgsea(stats, msigdb_sets, "bulk_LA_minus_NC", tm)
}), fill = TRUE)
mel_deg <- mel$deg_fixed
mel_all_deg <- fread(cmd = paste("gzip -cd", shQuote(paths$melanoma_deg)), sep = "\t", quote = "\"", data.table = TRUE)[cell_type == "CD8 T" & comparison == "C1_vs_C4"]
mel_stats <- mel_all_deg$moderated_t; names(mel_stats) <- mel_all_deg$gene
mel_fgsea <- run_fgsea(mel_stats, msigdb_sets, "melanoma_CD8_C1_minus_C4")
gsea_bridge <- merge(
  mel_fgsea[, .(feature_id, melanoma_NES = NES, melanoma_p = pval, melanoma_fdr = padj, melanoma_size = size)],
  bulk_fgsea[, .(feature_id, time, bulk_NES = NES, bulk_p = pval, bulk_fdr = padj, bulk_size = size)],
  by = "feature_id", allow.cartesian = TRUE
)
gsea_bridge[, `:=`(
  feature_label = unname(set_label[feature_id]),
  feature_layer = unname(set_layer[feature_id]),
  melanoma_direction = fifelse(melanoma_NES > 0, "positive", fifelse(melanoma_NES < 0, "negative", "zero")),
  bulk_direction = fifelse(bulk_NES > 0, "positive", fifelse(bulk_NES < 0, "negative", "zero"))
)]
gsea_bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
gsea_bridge[, evidence_class := ifelse(melanoma_fdr < 0.05 & bulk_fdr < 0.05, ifelse(direction_relation == "same_direction", "supported_concordance", "supported_divergence"), ifelse(direction_relation == "same_direction" & (melanoma_p < 0.05 | bulk_p < 0.05), "suggestive_concordance", "insufficient_or_mixed"))]
gsea_bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(gsea_bridge, layer_order, feature_id, time)
gsea_bridge[, layer_order := NULL]
gsea_bridge_summary <- gsea_bridge[, .(bridge_class = if (all(evidence_class == "supported_concordance")) "supported_concordance" else if (any(evidence_class == "supported_divergence")) "supported_divergence_or_mixed" else if (any(evidence_class == "suggestive_concordance")) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), melanoma_NES = first(melanoma_NES), melanoma_fdr = first(melanoma_fdr), mean_bulk_abs_NES = mean(abs(bulk_NES), na.rm = TRUE)), by = .(feature_id, feature_label, feature_layer)]

write_tsv(membership, file.path(out_dir, "fixed_gene_set_coverage_manifest.tsv"))
write_tsv(bulk$score_long, file.path(out_dir, "bulk_fixed_readout_ssgsea_long.tsv"))
write_tsv(bulk$score_contrasts, file.path(out_dir, "bulk_fixed_readout_ssgsea_contrasts.tsv"))
write_tsv(bulk$donor_delta, file.path(out_dir, "bulk_fixed_readout_donor_deltas.tsv"))
write_tsv(mel$score_long, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_long.tsv"))
write_tsv(mel$score_contrasts, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_c1_vs_c4.tsv"))
write_tsv(mel$group_summary, file.path(out_dir, "melanoma_cd8_fixed_readout_group_summary.tsv"))
write_tsv(mel$selection_manifest, file.path(out_dir, "melanoma_cd8_patient_selection_manifest.tsv"))
write_tsv(mel$deg_fixed, file.path(out_dir, "melanoma_cd8_fixed_readout_gene_deg.tsv"))
write_tsv(bridge, file.path(out_dir, "melanoma_cd8_bulk_bridge_by_time.tsv"))
write_tsv(bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_bridge_summary.tsv"))
write_tsv(gene_bridge, file.path(out_dir, "melanoma_cd8_bulk_gene_bridge.tsv"))
write_tsv(bulk_fgsea, file.path(out_dir, "bulk_msigdb_preranked_gsea.tsv"))
write_tsv(mel_fgsea, file.path(out_dir, "melanoma_cd8_msigdb_preranked_gsea.tsv"))
write_tsv(gsea_bridge, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_by_time.tsv"))
write_tsv(gsea_bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_summary.tsv"))
write_tsv(rbind(bulk$coverage[, `:=`(dataset = "bulk")], mel$coverage[, `:=`(dataset = "melanoma_CD8")]), file.path(out_dir, "score_coverage_qa.tsv"))
write_tsv(data.table(
  check_id = c("bulk_factorial_samples", "bulk_design_rank", "melanoma_c1c4_samples", "melanoma_c1c4_patients", "melanoma_unique_patient_units", "fixed_readout_count", "bulk_coverage", "melanoma_coverage", "all_ssgsea_finite", "bridge_rows", "bulk_fgsea_rows", "melanoma_fgsea_rows"),
  status = c(ifelse(nrow(bulk$sample_sheet) == 36L, "PASS", "HOLD"), ifelse(bulk$design_rank == ncol(bulk$design), "PASS", "HOLD"), ifelse(mel$n_samples >= 6L, "PASS", "HOLD"), ifelse(mel$n_patients >= 6L, "PASS", "HOLD"), ifelse(!anyDuplicated(mel$sample_sheet$patient), "PASS", "HOLD"), ifelse(length(all_sets) == 21L, "PASS", "HOLD"), ifelse(all(membership$bulk_score_status == "PASS"), "PASS", "HOLD"), ifelse(all(membership$melanoma_score_status == "PASS"), "PASS", "HOLD"), ifelse(all(is.finite(bulk$score_long$ssgsea_score)) && all(is.finite(mel$score_long$ssgsea_score)), "PASS", "HOLD"), ifelse(nrow(bridge) == length(all_sets) * 4L, "PASS", "HOLD"), ifelse(nrow(bulk_fgsea) == length(msigdb_sets) * 4L, "PASS", "HOLD"), ifelse(nrow(mel_fgsea) == length(msigdb_sets), "PASS", "HOLD")),
  detail = c(nrow(bulk$sample_sheet), paste(bulk$design_rank, "of", ncol(bulk$design)), mel$n_samples, mel$n_patients, uniqueN(mel$sample_sheet$patient), length(all_sets), sprintf("min %.3f", min(membership$bulk_coverage)), sprintf("min %.3f", min(membership$melanoma_coverage)), paste(nrow(bulk$score_long), nrow(mel$score_long), sep = " + "), nrow(bridge), nrow(bulk_fgsea), nrow(mel_fgsea))
), file.path(out_dir, "qa_summary.tsv"))
write_tsv(data.table(object = c("bulk_count", "melanoma_cd8_pseudobulk"), analysis_unit = c("donor x condition x time sample", "sample x native CD8 T pseudobulk; patient biological unit"), primary_comparison = c("time-matched LA_minus_NC", "C1_minus_C4"), n_independent_units = c(length(unique(bulk$sample_sheet$Donor)), mel$n_patients), interpretation = c("CD8 bulk aggregate response", "melanoma CD8 state association")), file.path(out_dir, "analysis_unit_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
writeLines(c("status=PASS_CANDIDATE", paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)), paste0("fixed_readouts=", length(all_sets)), paste0("bulk_samples=", nrow(bulk$sample_sheet), ";bulk_donors=", length(unique(bulk$sample_sheet$Donor))), paste0("melanoma_cd8_samples=", mel$n_samples, ";melanoma_patients=", mel$n_patients, ";selection=", mel$selection_mode), paste0("bridge_rows=", nrow(bridge))), file.path(out_dir, "run_summary.txt"))
cat("STATUS\tPASS_CANDIDATE\n")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
}

# ----------------------------- v0.2 execution -----------------------------
bulk <- build_bulk()
mel <- build_melanoma()

# donor-specific rank enrichment is descriptive and supplies a direction gate;
# it does not replace the formal limma rank used for the bulk GSEA.
build_bulk_donor_gsea <- function(bulk_obj) {
  out <- list(); k <- 0L
  for (donor in unique(as.character(bulk_obj$ann$Donor))) for (tm in c("1h", "6h", "24h", "48h")) {
    la <- bulk_obj$ann$sample_id[as.character(bulk_obj$ann$Donor) == donor & as.character(bulk_obj$ann$Time) == tm & as.character(bulk_obj$ann$Treatment) == "LA"]
    nc <- bulk_obj$ann$sample_id[as.character(bulk_obj$ann$Donor) == donor & as.character(bulk_obj$ann$Time) == tm & as.character(bulk_obj$ann$Treatment) == "NC"]
    stop_if(length(la) != 1L || length(nc) != 1L, paste("donor rank pair missing", donor, tm))
    stats <- bulk_obj$logcpm[, la] - bulk_obj$logcpm[, nc]
    names(stats) <- rownames(bulk_obj$logcpm)
    k <- k + 1L
    out[[k]] <- run_fgsea(stats, all_sets, "bulk_donor_LA_minus_NC", tm, donor, analysis_seed)
  }
  rbindlist(out, fill = TRUE)
}

bulk_donor_gsea_raw <- build_bulk_donor_gsea(bulk)
bulk_fgsea_raw <- rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) {
  d <- bulk_de_stats(tm); stats <- d$rank_score; names(stats) <- d$gene
  run_fgsea(stats, all_sets, "bulk_LA_minus_NC", tm, seed = analysis_seed)
}), fill = TRUE)
mel_all_deg <- fread(cmd = paste("gzip -cd", shQuote(paths$melanoma_deg)), sep = "\t", quote = "\"", data.table = TRUE)
mel_all_deg <- mel_all_deg[cell_type == "CD8 T" & comparison == "C1_vs_C4"]
mel_stats <- mel_all_deg$moderated_t
names(mel_stats) <- mel_all_deg$gene
mel_fgsea_raw <- run_fgsea(mel_stats, all_sets, "melanoma_CD8_C1_minus_C4", seed = analysis_seed)

flatten_gsea <- function(x) {
  x <- copy(as.data.table(x))
  x[, leadingEdge_genes := vapply(seq_len(.N), function(i) paste(as.character(leadingEdge[[i]]), collapse = ";"), character(1))]
  x[, leadingEdge := NULL]
  x[]
}

gsea_repro_one <- function(stats, context, time = NA_character_) {
  a <- run_fgsea(stats, all_sets, context, time, seed = analysis_seed)
  b <- run_fgsea(stats, all_sets, context, time, seed = boundary_seed)
  out <- merge(a[, .(feature_id, NES_a = NES, FDR_a = padj, p_a = pval)], b[, .(feature_id, NES_b = NES, FDR_b = padj, p_b = pval)], by = "feature_id", all = TRUE)
  out[, `:=`(context = context, time = time, seed_a = analysis_seed, seed_b = boundary_seed, abs_NES_delta = abs(NES_a - NES_b), abs_FDR_delta = abs(FDR_a - FDR_b), boundary_a = FDR_a >= 0.02 & FDR_a <= 0.10, boundary_b = FDR_b >= 0.02 & FDR_b <= 0.10, boundary_changed = (FDR_a < 0.05) != (FDR_b < 0.05))]
  out[]
}

gsea_repro <- rbindlist(c(
  list(gsea_repro_one(mel_stats, "melanoma_CD8_C1_minus_C4")),
  lapply(c("1h", "6h", "24h", "48h"), function(tm) { d <- bulk_de_stats(tm); stats <- d$rank_score; names(stats) <- d$gene; gsea_repro_one(stats, "bulk_LA_minus_NC", tm) })
), fill = TRUE)

# Program-level ssGSEA bridge with explicit donor summaries and t-based intervals.
bulk_primary <- bulk$score_contrasts[contrast == "LA_minus_NC"]
mel_primary <- mel$score_contrasts
bridge <- merge(
  mel_primary[, .(feature_id, melanoma_logFC = logFC, melanoma_p = P.Value, melanoma_fdr = adj.P.Val, melanoma_CI_low = CI_low, melanoma_CI_high = CI_high)],
  bulk_primary[, .(feature_id, time, bulk_logFC = logFC, bulk_p = P.Value, bulk_fdr = adj.P.Val, bulk_CI_low = CI_low, bulk_CI_high = CI_high)],
  by = "feature_id", allow.cartesian = TRUE
)
donor_support <- merge(bulk$donor_delta[contrast == "LA_minus_NC"], mel_primary[, .(feature_id, melanoma_logFC = logFC)], by = "feature_id", all.x = TRUE)
donor_support[, donor_relation := fifelse(donor_delta == 0 | melanoma_logFC == 0, "zero", fifelse(sign(donor_delta) == sign(melanoma_logFC), "same", "opposite"))]
support_summary <- donor_support[, .(n_donors = uniqueN(donor), n_donors_same_direction = sum(donor_relation == "same"), n_donors_opposite_direction = sum(donor_relation == "opposite"), n_donors_zero = sum(donor_relation == "zero")), by = .(feature_id, time)]
support_summary[, `:=`(donor_same_fraction = n_donors_same_direction / n_donors, donor_relation_consistency = fifelse(n_donors_same_direction == n_donors, "same_3of3", fifelse(n_donors_opposite_direction == n_donors, "opposite_3of3", fifelse(pmax(n_donors_same_direction, n_donors_opposite_direction) >= 2, "2of3_consistent", "mixed"))))]
bridge <- merge(bridge, support_summary, by = c("feature_id", "time"), all.x = TRUE)
bridge[, `:=`(feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]), melanoma_direction = fifelse(melanoma_logFC > 0, "positive", fifelse(melanoma_logFC < 0, "negative", "zero")), bulk_direction = fifelse(bulk_logFC > 0, "positive", fifelse(bulk_logFC < 0, "negative", "zero")))]
bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
classify_bridge <- function(mel_fdr, bulk_fdr, mel_p, bulk_p, relation, n_same, n_opp) {
  if (is.na(relation) || relation == "insufficient_evidence") return("insufficient_evidence")
  if (mel_fdr < 0.05 && bulk_fdr < 0.05 && relation == "same_direction" && n_same == 3) return("supported_concordance")
  if (mel_fdr < 0.05 && bulk_fdr < 0.05 && relation == "opposite_direction" && n_opp == 3) return("supported_divergence")
  if (relation == "same_direction" && (mel_p < 0.05 || bulk_p < 0.05) && n_same >= 2) return("suggestive_concordance")
  if (relation == "opposite_direction" && (mel_p < 0.05 || bulk_p < 0.05) && n_opp >= 2) return("suggestive_divergence")
  if (max(n_same, n_opp) >= 2) return("donor_consistent_direction_only")
  "time_dependent_or_mixed"
}
bridge[, evidence_class := mapply(classify_bridge, melanoma_fdr, bulk_fdr, melanoma_p, bulk_p, direction_relation, n_donors_same_direction, n_donors_opposite_direction)]
bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(bridge, layer_order, feature_id, time)
bridge[, layer_order := NULL]
bridge_summary <- bridge[, .(bridge_class = if (any(evidence_class == "supported_concordance")) "supported_concordance_present" else if (any(evidence_class == "supported_divergence")) "supported_divergence_present" else if (any(evidence_class %in% c("suggestive_concordance", "suggestive_divergence"))) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), mean_abs_bulk_logFC = mean(abs(bulk_logFC), na.rm = TRUE), melanoma_logFC = first(melanoma_logFC), melanoma_fdr = first(melanoma_fdr)), by = .(feature_id, feature_label, feature_layer)]

# GSEA bridge uses formal bulk/melanoma ranks plus donor-specific directional consistency.
gsea_bridge <- merge(
  mel_fgsea_raw[, .(feature_id, melanoma_NES = NES, melanoma_p = pval, melanoma_fdr = padj, melanoma_size = size)],
  bulk_fgsea_raw[, .(feature_id, time, bulk_NES = NES, bulk_p = pval, bulk_fdr = padj, bulk_size = size)],
  by = "feature_id", allow.cartesian = TRUE
)
gsea_donor <- merge(bulk_donor_gsea_raw[, .(feature_id, time, donor, donor_NES = NES)], mel_fgsea_raw[, .(feature_id, melanoma_NES = NES)], by = "feature_id", all.x = TRUE)
gsea_donor[, donor_relation := fifelse(donor_NES == 0 | melanoma_NES == 0, "zero", fifelse(sign(donor_NES) == sign(melanoma_NES), "same", "opposite"))]
gsea_donor_summary <- gsea_donor[, .(n_donors = uniqueN(donor), n_donors_same_direction = sum(donor_relation == "same"), n_donors_opposite_direction = sum(donor_relation == "opposite"), n_donors_zero = sum(donor_relation == "zero")), by = .(feature_id, time)]
gsea_donor_summary[, `:=`(donor_same_fraction = n_donors_same_direction / n_donors, donor_relation_consistency = fifelse(n_donors_same_direction == n_donors, "same_3of3", fifelse(n_donors_opposite_direction == n_donors, "opposite_3of3", fifelse(pmax(n_donors_same_direction, n_donors_opposite_direction) >= 2, "2of3_consistent", "mixed"))))]
gsea_bridge <- merge(gsea_bridge, gsea_donor_summary, by = c("feature_id", "time"), all.x = TRUE)
gsea_bridge[, `:=`(feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]), melanoma_direction = fifelse(melanoma_NES > 0, "positive", fifelse(melanoma_NES < 0, "negative", "zero")), bulk_direction = fifelse(bulk_NES > 0, "positive", fifelse(bulk_NES < 0, "negative", "zero")))]
gsea_bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
gsea_bridge[, evidence_class := mapply(classify_bridge, melanoma_fdr, bulk_fdr, melanoma_p, bulk_p, direction_relation, n_donors_same_direction, n_donors_opposite_direction)]
gsea_bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(gsea_bridge, layer_order, feature_id, time)
gsea_bridge[, layer_order := NULL]
gsea_bridge_summary <- gsea_bridge[, .(bridge_class = if (any(evidence_class == "supported_concordance")) "supported_concordance_present" else if (any(evidence_class == "supported_divergence")) "supported_divergence_present" else if (any(evidence_class %in% c("suggestive_concordance", "suggestive_divergence"))) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), melanoma_NES = first(melanoma_NES), melanoma_fdr = first(melanoma_fdr), mean_bulk_abs_NES = mean(abs(bulk_NES), na.rm = TRUE)), by = .(feature_id, feature_label, feature_layer)]

# 方法方向审计：同一 feature 的 ssGSEA 与 GSEA 不强行视为同一结果。
method_direction_audit <- merge(
  mel_primary[, .(feature_id, context = "melanoma_CD8", time = NA_character_, ssgsea_effect = logFC, ssgsea_direction = sign(logFC))],
  mel_fgsea_raw[, .(feature_id, gsea_effect = NES, gsea_direction = sign(NES))], by = "feature_id", all = TRUE
)
method_direction_audit[, `:=`(method_direction_discordant = !is.na(ssgsea_direction) & !is.na(gsea_direction) & ssgsea_direction != 0 & gsea_direction != 0 & ssgsea_direction != gsea_direction, interpretation = "ssGSEA is sample-level score; GSEA is rank enrichment; discordance requires review")]
bulk_method <- merge(bulk_primary[, .(feature_id, time, ssgsea_effect = logFC, ssgsea_direction = sign(logFC))], bulk_fgsea_raw[, .(feature_id, time, gsea_effect = NES, gsea_direction = sign(NES))], by = c("feature_id", "time"), all = TRUE)
bulk_method[, `:=`(context = "bulk_LA_minus_NC", method_direction_discordant = !is.na(ssgsea_direction) & !is.na(gsea_direction) & ssgsea_direction != 0 & gsea_direction != 0 & ssgsea_direction != gsea_direction, interpretation = "ssGSEA is sample-level score; GSEA is rank enrichment; discordance requires review")]
method_direction_audit <- rbind(method_direction_audit, bulk_method, fill = TRUE)
memory_direction_audit <- method_direction_audit[feature_id %in% c("MSIGDB::GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN", "MSIGDB::GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP")]
memory_direction_audit[, semantic_note := fifelse(grepl("_DN$", feature_id), "DN set is lower in effector than memory; interpret as the memory-enriched side of the paired signature", "UP set is higher in effector than memory; interpret as the effector-enriched side of the paired signature")]

# Leading-edge 拆解使用本地完整 MSigDB 诊断集合，不另造 CD8 marker 集合。
diag_specs <- data.table(
  category = c("cell_cycle_E2F", "cell_cycle_G2M", "IFN_gamma", "T_cell_cytotoxicity", "TCR_signaling", "hypoxia", "UPR", "glycolysis"),
  source_file_key = c("hallmark", "hallmark", "hallmark", "go_bp", "reactome", "hallmark", "hallmark", "hallmark"),
  feature_id = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "GOBP_T_CELL_MEDIATED_CYTOTOXICITY", "REACTOME_TCR_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_GLYCOLYSIS")
)
diag_specs[, exists := mapply(function(key, id) id %in% names(gmt_cache[[key]]$sets), source_file_key, feature_id)]
diag_specs[, n_genes_source := mapply(function(key, id) if (id %in% names(gmt_cache[[key]]$sets)) length(gmt_cache[[key]]$sets[[id]]) else NA_integer_, source_file_key, feature_id)]
diag_specs[, source_path := unname(msigdb_files[source_file_key])]
diag_existing <- diag_specs[exists == TRUE]
diag_sets <- setNames(lapply(seq_len(nrow(diag_existing)), function(i) gmt_cache[[diag_existing$source_file_key[i]]]$sets[[diag_existing$feature_id[i]]]), diag_existing$category)
gene_to_cat <- list()
for (cat in names(diag_sets)) for (g in diag_sets[[cat]]) gene_to_cat[[g]] <- unique(c(gene_to_cat[[g]], cat))

gsea_primary_raw <- rbind(mel_fgsea_raw, bulk_fgsea_raw, fill = TRUE)
edge_rows <- list()
k <- 0L
for (i in seq_len(nrow(gsea_primary_raw))) {
  genes <- as.character(gsea_primary_raw$leadingEdge[[i]])
  if (!length(genes)) next
  k <- k + 1L
  edge_rows[[k]] <- data.table(context = gsea_primary_raw$context[i], time = gsea_primary_raw$time[i], feature_id = gsea_primary_raw$feature_id[i], gene = genes, rank_in_leading_edge = seq_along(genes))
}
edge_long <- rbindlist(edge_rows, fill = TRUE)
if (nrow(edge_long)) {
  edge_cat <- rbindlist(lapply(seq_len(nrow(edge_long)), function(i) {
    cats <- gene_to_cat[[edge_long$gene[i]]]; if (is.null(cats) || !length(cats)) cats <- "unclassified"
    data.table(context = edge_long$context[i], time = edge_long$time[i], feature_id = edge_long$feature_id[i], gene = edge_long$gene[i], rank_in_leading_edge = edge_long$rank_in_leading_edge[i], category = cats)
  }), fill = TRUE)
  edge_summary <- edge_cat[, .(n_genes = uniqueN(gene), mean_rank = mean(rank_in_leading_edge), genes = paste(head(unique(gene), 25L), collapse = ";")), by = .(context, time, feature_id, category)]
} else {
  edge_cat <- data.table(); edge_summary <- data.table()
}

# 基因级桥接只用于候选优先级；不把名义方向直接称为验证。
mel_gene <- mel$deg_fixed[order(fdr, p_value)][!duplicated(gene)]
gene_bridge <- merge(mel_gene[, .(gene, melanoma_logFC = log2fc_c1_minus_ref, melanoma_p = p_value, melanoma_fdr = fdr)], rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) { d <- bulk_de_stats(tm); d[, time := tm]; d[gene %in% unique(unlist(all_sets))] })), by = "gene", allow.cartesian = TRUE)
gene_bridge[, direction_relation := fifelse(sign(melanoma_logFC) == sign(logFC), "same_direction", "opposite_direction")]
gene_bridge[, set_membership := vapply(gene, function(g) paste(names(all_sets)[vapply(all_sets, function(s) g %in% s, logical(1))], collapse = ";"), character(1))]
gene_bridge[, candidate_gate := fifelse(melanoma_fdr < 0.05 & FDR < 0.05 & direction_relation == "same_direction", "strict_two_sided_gate", "not_passed")]

membership <- merge(feature_meta, source_size, by = "feature_id", all = TRUE, sort = FALSE)
membership[, source_path := unname(source_path[feature_id])]
membership <- merge(membership, bulk$coverage[, .(feature_id, n_genes_bulk = n_genes_present, bulk_coverage = coverage)], by = "feature_id", all = TRUE)
membership <- merge(membership, mel$coverage[, .(feature_id, n_genes_melanoma = n_genes_present, melanoma_coverage = coverage)], by = "feature_id", all = TRUE)
membership[, `:=`(bulk_score_status = ifelse(n_genes_bulk >= 3L, "PASS", "HOLD"), melanoma_score_status = ifelse(n_genes_melanoma >= 3L, "PASS", "HOLD"))]

write_tsv(membership, file.path(out_dir, "fixed_gene_set_coverage_manifest.tsv"))
write_tsv(rbind(bulk$mapping, mel$mapping, fill = TRUE), file.path(out_dir, "score_id_mapping_qa.tsv"))
write_tsv(bulk$score_long, file.path(out_dir, "bulk_fixed_readout_ssgsea_long.tsv"))
write_tsv(bulk$score_contrasts, file.path(out_dir, "bulk_fixed_readout_ssgsea_contrasts.tsv"))
write_tsv(bulk$donor_delta, file.path(out_dir, "bulk_fixed_readout_donor_deltas.tsv"))
write_tsv(bulk$donor_summary, file.path(out_dir, "bulk_fixed_readout_donor_summary.tsv"))
write_tsv(mel$score_long, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_long.tsv"))
write_tsv(mel$score_contrasts, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_c1_vs_c4.tsv"))
write_tsv(mel$group_summary, file.path(out_dir, "melanoma_cd8_fixed_readout_group_summary.tsv"))
write_tsv(mel$selection_manifest, file.path(out_dir, "melanoma_cd8_patient_selection_manifest.tsv"))
write_tsv(mel$deg_fixed, file.path(out_dir, "melanoma_cd8_fixed_readout_gene_deg.tsv"))
write_tsv(bridge, file.path(out_dir, "melanoma_cd8_bulk_bridge_by_time.tsv"))
write_tsv(bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_bridge_summary.tsv"))
write_tsv(gene_bridge, file.path(out_dir, "melanoma_cd8_bulk_gene_bridge.tsv"))
write_tsv(flatten_gsea(bulk_fgsea_raw), file.path(out_dir, "bulk_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(mel_fgsea_raw), file.path(out_dir, "melanoma_cd8_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(bulk_donor_gsea_raw), file.path(out_dir, "bulk_donor_directional_preranked_gsea.tsv"))
write_tsv(gsea_bridge, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_by_time.tsv"))
write_tsv(gsea_bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_summary.tsv"))
write_tsv(gsea_repro, file.path(out_dir, "gsea_reproducibility_boundary.tsv"))
write_tsv(method_direction_audit, file.path(out_dir, "ssgsea_gsea_direction_audit.tsv"))
write_tsv(memory_direction_audit, file.path(out_dir, "memory_direction_audit.tsv"))
write_tsv(diag_specs, file.path(out_dir, "diagnostic_gene_set_manifest.tsv"))
write_tsv(edge_long, file.path(out_dir, "gsea_leading_edge_long.tsv"))
write_tsv(edge_cat, file.path(out_dir, "gsea_leading_edge_category_long.tsv"))
write_tsv(edge_summary, file.path(out_dir, "gsea_leading_edge_category_summary.tsv"))
write_tsv(gsea_repro[boundary_a | boundary_b | boundary_changed], file.path(out_dir, "gsea_boundary_focus.tsv"))
write_tsv(data.table(
  object = c("bulk_count", "melanoma_cd8_pseudobulk"),
  analysis_unit = c("donor x condition x time sample", "sample x native CD8 T pseudobulk; patient biological unit"),
  primary_comparison = c("time-matched LA_minus_NC", "C1_minus_C4"),
  n_independent_units = c(length(unique(bulk$sample_sheet$Donor)), mel$n_patients),
  interpretation = c("CD8 bulk aggregate response", "melanoma CD8 state association")
), file.path(out_dir, "analysis_unit_manifest.tsv"))

mapping_ok <- all(bulk$mapping$mapping_ok) && all(mel$mapping$mapping_ok)
gsea_rows_ok <- nrow(bulk_fgsea_raw) == length(all_sets) * 4L && nrow(mel_fgsea_raw) == length(all_sets)
donor_rows_ok <- nrow(bulk$donor_delta[contrast == "LA_minus_NC"]) == length(all_sets) * 4L * 3L
repro_ok <- nrow(gsea_repro) > 0L && all(is.finite(gsea_repro$abs_NES_delta)) && all(is.finite(gsea_repro$abs_FDR_delta))
qa <- data.table(
  check_id = c("bulk_factorial_samples", "bulk_design_rank", "melanoma_c1c4_samples", "melanoma_c1c4_patients", "fixed_readout_count", "score_id_mapping_exact", "bulk_coverage", "melanoma_coverage", "all_ssgsea_finite", "bulk_primary_donor_rows", "gsea_rows_all_fixed", "donor_gsea_rows", "gsea_seed_reproducible_run", "leading_edge_table"),
  status = c(
    ifelse(nrow(bulk$sample_sheet) == 36L, "PASS", "HOLD"),
    ifelse(bulk$design_rank == ncol(bulk$design), "PASS", "HOLD"),
    ifelse(mel$n_samples >= 6L, "PASS", "HOLD"),
    ifelse(mel$n_patients >= 6L, "PASS", "HOLD"),
    ifelse(length(all_sets) == 21L, "PASS", "HOLD"),
    ifelse(mapping_ok, "PASS", "HOLD"),
    ifelse(all(membership$bulk_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(membership$melanoma_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(is.finite(bulk$score_long$ssgsea_score)) && all(is.finite(mel$score_long$ssgsea_score)), "PASS", "HOLD"),
    ifelse(donor_rows_ok, "PASS", "HOLD"),
    ifelse(gsea_rows_ok, "PASS", "HOLD"),
    ifelse(nrow(bulk_donor_gsea_raw) == length(all_sets) * 4L * 3L, "PASS", "HOLD"),
    ifelse(repro_ok, "PASS", "HOLD"),
    ifelse(nrow(edge_long) > 0L, "PASS", "HOLD")
  ),
  detail = c(nrow(bulk$sample_sheet), paste(bulk$design_rank, "of", ncol(bulk$design)), mel$n_samples, mel$n_patients, length(all_sets), paste(all(bulk$mapping$mapping_ok), all(mel$mapping$mapping_ok), sep = ";"), sprintf("min %.3f", min(membership$bulk_coverage)), sprintf("min %.3f", min(membership$melanoma_coverage)), paste(nrow(bulk$score_long), nrow(mel$score_long), sep = " + "), nrow(bulk$donor_delta[contrast == "LA_minus_NC"]), paste(nrow(bulk_fgsea_raw), nrow(mel_fgsea_raw), sep = " + "), nrow(bulk_donor_gsea_raw), sprintf("max |delta NES| %.4g; max |delta FDR| %.4g", max(gsea_repro$abs_NES_delta), max(gsea_repro$abs_FDR_delta)), nrow(edge_long))
)
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
overall_status <- if (all(qa$status == "PASS")) "PASS_CANDIDATE" else "HOLD_CANDIDATE"
writeLines(c(paste0("status=", overall_status), paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)), paste0("analysis_seed=", analysis_seed), paste0("boundary_seed=", boundary_seed), paste0("fixed_readouts=", length(all_sets)), paste0("bulk_samples=", nrow(bulk$sample_sheet), ";bulk_donors=", length(unique(bulk$sample_sheet$Donor))), paste0("melanoma_cd8_samples=", mel$n_samples, ";melanoma_patients=", mel$n_patients, ";selection=", mel$selection_mode), paste0("gsea_rows=", nrow(bulk_fgsea_raw), "+", nrow(mel_fgsea_raw)), paste0("leading_edge_rows=", nrow(edge_long))), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("STATUS\t", overall_status, "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
}

# Purpose: Fit donor-adjusted CD8 treatment-time responses and fixed programme enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/04_hcl_control_phenotype.R
stage_cd8_perturbation_04_hcl_control_phenotype <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("count_rdata","lm12_membership")]
gmt_basenames <- c("h.all.v2025.1.Hs.symbols.gmt","c2.cp.reactome.v2025.1.Hs.symbols.gmt","c5.go.bp.v2024.1.Hs.symbols.gmt","c7.all.v2025.1.Hs.symbols.gmt")
de_contrasts <- c("LA_minus_NC","LA_minus_HCl","HCl_minus_NC")
de_basenames <- as.vector(outer(de_contrasts, c("1h","6h","24h","48h"), function(x,y) paste0("de_",x,"_",y,".tsv")))
required_files <- c(unname(input_paths), file.path(unname(paths["msigdb_dir"]),gmt_basenames), file.path(unname(paths["bulk_de_dir"]),de_basenames))
bridge_basenames <- c("fixed_gene_set_manifest.tsv","bulk_all_fixed_preranked_gsea.tsv","bulk_fixed_readout_ssgsea_contrasts.tsv","bulk_fixed_readout_donor_summary.tsv")
required_files <- c(required_files,file.path(unname(paths["bridge_dir"]),bridge_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("msigdb_dir","bulk_de_dir","bridge_dir")])
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

cfg_paths <- paths

# ---- Ordinary-only fixed-resource path routing; no set/content changes ----
relocate_fixed_source <- function(layer, source_path) {
  if (length(layer) != 1L || is.na(layer) || length(source_path) != 1L ||
      is.na(source_path) || !nzchar(source_path)) stop("Invalid fixed-resource identity.")
  if (layer == "LM12_FIG2_FROZEN") return(unname(cfg_paths["lm12_membership"]))
  if (layer != "MSIGDB_FULL_SET") stop("Unknown fixed-resource layer: ", layer)
  file_name <- basename(source_path)
  if (!file_name %in% gmt_basenames) stop("Unregistered fixed GMT basename: ", file_name)
  file.path(unname(cfg_paths["msigdb_dir"]), file_name)
}

# ---- Retained complete source body; documented path adaptations only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(fgsea)
})

analysis_seed <- 20260916L
input_path <- unname(cfg_paths["count_rdata"])
de_root <- unname(cfg_paths["bulk_de_dir"])
bridge_root <- unname(cfg_paths["bridge_dir"])
out_dir <- unname(cfg_paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
`write_tsv` <- bind_stage_helper(shared_write_tsv_046, environment())
`stop_if` <- bind_stage_helper(shared_stop_if_003, environment())
times <- c("1h", "6h", "24h", "48h")
contrasts <- c("LA_minus_NC", "HCl_minus_NC", "LA_minus_HCl")

load_count <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  stop_if(!all(c("count", "sample_anno") %in% ls(e, all.names = TRUE)), "count/sample_anno missing")
  list(count = get("count", e), anno = as.data.table(get("sample_anno", e), keep.rownames = "sample_id"))
}

obj <- load_count(input_path)
count <- obj$count
anno <- obj$anno
stop_if(!identical(colnames(count), anno$sample_id), "count and annotation order mismatch")
anno$sample_id <- as.character(anno$sample_id)
anno$Donor <- as.character(anno$Donor)
anno$Treatment <- as.character(anno$Treatment)
anno$Time <- as.character(anno$Time)
observed_levels <- rbindlist(list(
  data.table(level_type = "Treatment", level = paste(sort(unique(anno$Treatment)), collapse = ";")),
  data.table(level_type = "Time", level = paste(sort(unique(anno$Time)), collapse = ";")),
  data.table(level_type = "Donor", level = paste(sort(unique(anno$Donor)), collapse = ";"))
))
write_tsv(observed_levels, file.path(out_dir, "observed_design_levels.tsv"))

factorial <- anno[Treatment %in% c("NC", "HCl", "LA") & Time %in% times]
factorial$Treatment <- factor(factorial$Treatment, levels = c("NC", "HCl", "LA"))
factorial$Time <- factor(factorial$Time, levels = times)
factorial$Donor <- factor(factorial$Donor)
stop_if(nrow(factorial) != 36L, paste("unexpected factorial sample count:", nrow(factorial)))
stop_if(anyDuplicated(factorial[, .(Donor, Treatment, Time)]) > 0L, "duplicate donor-treatment-time sample")
design <- model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = factorial)
stop_if(qr(design)$rank != ncol(design), "factorial design is rank deficient")
write_tsv(factorial[, .(sample_id, Donor = as.character(Donor), Treatment = as.character(Treatment), Time = as.character(Time))],
          file.path(out_dir, "factorial_sample_manifest.tsv"))

fixed_manifest <- fread(file.path(bridge_root, "fixed_gene_set_manifest.tsv"))
stop_if(!all(c("feature_id", "feature_label", "feature_layer", "source_path") %in% names(fixed_manifest)),
        "fixed gene-set manifest columns missing")
# Ordinary run: map historical manifest resource paths to configured inputs.
fixed_manifest$source_path <- vapply(seq_len(nrow(fixed_manifest)), function(i)
  relocate_fixed_source(as.character(fixed_manifest$feature_layer[i]),
                        as.character(fixed_manifest$source_path[i])), character(1L))
lm12_path <- fixed_manifest[feature_layer == "LM12_FIG2_FROZEN", unique(source_path)][1L]
lm12_tab <- fread(lm12_path)
stop_if(!all(c("feature_id", "gene") %in% names(lm12_tab)), "LM12 membership columns missing")
lm12_ids <- unique(as.character(lm12_tab$feature_id))
lm12_sets <- lapply(lm12_ids, function(id) unique(as.character(lm12_tab[feature_id == id, gene])))
names(lm12_sets) <- paste0("LM12::", lm12_ids)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  toks <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(toks, function(x) x[[1L]], character(1))
  sets <- lapply(toks, function(x) unique(x[-seq_len(min(2L, length(x)))]))
  names(sets) <- ids
  sets[!duplicated(names(sets))]
}
msig_rows <- fixed_manifest[feature_layer == "MSIGDB_FULL_SET"]
gmt_cache <- list()
msig_sets <- list()
for (i in seq_len(nrow(msig_rows))) {
  path <- msig_rows$source_path[i]
  if (is.null(gmt_cache[[path]])) gmt_cache[[path]] <- read_gmt(path)
  raw_id <- msig_rows$feature_id[i]
  stop_if(!raw_id %in% names(gmt_cache[[path]]), paste("MSigDB set missing:", raw_id))
  msig_sets[[paste0("MSIGDB::", raw_id)]] <- gmt_cache[[path]][[raw_id]]
}
all_sets <- c(lm12_sets, msig_sets)
set_meta <- rbind(
  fixed_manifest[feature_layer == "LM12_FIG2_FROZEN", .(
    feature_id = paste0("LM12::", feature_id), feature_label, feature_layer, source_collection,
    n_genes_source, source_path
  )],
  fixed_manifest[feature_layer == "MSIGDB_FULL_SET", .(
    feature_id = paste0("MSIGDB::", feature_id), feature_label, feature_layer, source_collection,
    n_genes_source, source_path
  )]
)
set_meta <- set_meta[match(names(all_sets), feature_id)]
stop_if(anyNA(set_meta$feature_id) || length(all_sets) != nrow(set_meta), "set metadata mapping failed")
write_tsv(set_meta, file.path(out_dir, "fixed_gene_set_manifest_reused.tsv"))

load_de <- function(contrast, tm) {
  path <- file.path(de_root, paste0("de_", contrast, "_", tm, ".tsv"))
  stop_if(!file.exists(path), paste("missing DE input:", path))
  d <- fread(path)
  stop_if(!all(c("gene_id", "logFC", "PValue", "FDR") %in% names(d)), paste("DE columns missing:", path))
  d$rank_score <- sign(d$logFC) * -log10(pmax(d$PValue, .Machine$double.xmin))
  d <- d[is.finite(d$rank_score) & !duplicated(d$gene_id)]
  stats <- d$rank_score
  names(stats) <- d$gene_id
  list(table = d, stats = stats)
}

flatten_gsea <- function(x) {
  x <- as.data.table(x)
  if (nrow(x) > 0L && "leadingEdge" %in% names(x)) {
    x$leadingEdge_genes <- vapply(x$leadingEdge, function(g) paste(as.character(g), collapse = ";"), character(1))
    x$leadingEdge <- NULL
  }
  x
}

run_gsea_one <- function(stats, contrast, tm) {
  present <- lapply(all_sets, function(g) intersect(unique(g), names(stats)))
  n_present <- lengths(present)
  testable <- n_present >= 3L & n_present <= 5000L
  out <- data.table()
  if (any(testable)) {
    set.seed(analysis_seed)
    fg <- fgsea::fgseaMultilevel(pathways = all_sets[testable], stats = stats,
                                  minSize = 3L, maxSize = 5000L, eps = 0)
    out <- flatten_gsea(fg)
    out$feature_id <- as.character(out$pathway)
    out$pathway <- NULL
    out$contrast <- contrast
    out$time <- tm
    out$n_genes_tested <- n_present[match(out$feature_id, names(all_sets))]
    out$status <- "PASS"
  }
  hold_ids <- names(all_sets)[!testable]
  if (length(hold_ids)) {
    hold <- data.table(feature_id = hold_ids, pval = NA_real_, padj = NA_real_, log2err = NA_real_,
                       ES = NA_real_, NES = NA_real_, size = n_present[hold_ids],
                       contrast = contrast, time = tm, n_genes_tested = n_present[hold_ids],
                       status = "HOLD", leadingEdge_genes = NA_character_)
    out <- rbindlist(list(out, hold), fill = TRUE)
  }
  out <- merge(out, set_meta[, .(feature_id, feature_label, feature_layer, source_collection,
                                 n_genes_source)], by = "feature_id", all.x = TRUE, sort = FALSE)
  out[]
}

existing_la <- fread(file.path(bridge_root, "bulk_all_fixed_preranked_gsea.tsv"))
existing_la <- existing_la[context == "bulk_LA_minus_NC" & time %in% times]
existing_la$contrast <- "LA_minus_NC"
existing_la$status <- "REUSED_FROZEN"
existing_la$n_genes_source <- set_meta$n_genes_source[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$feature_label <- set_meta$feature_label[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$feature_layer <- set_meta$feature_layer[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$source_collection <- set_meta$source_collection[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$n_genes_tested <- existing_la$size
existing_la <- existing_la[, .(contrast, time, feature_id, feature_label, feature_layer,
                               source_collection, n_genes_source, n_genes_tested, status,
                               pval, padj, log2err, ES, NES, size, leadingEdge_genes)]

new_gsea <- rbindlist(lapply(contrasts[contrasts != "LA_minus_NC"], function(ct) {
  rbindlist(lapply(times, function(tm) run_gsea_one(load_de(ct, tm)$stats, ct, tm)), fill = TRUE)
}), fill = TRUE)
gsea_three <- rbindlist(list(existing_la, new_gsea), fill = TRUE)
gsea_three <- gsea_three[order(match(contrast, contrasts), match(time, times),
                               match(feature_id, set_meta$feature_id))]
write_tsv(gsea_three, file.path(out_dir, "gsea_three_contrasts.tsv"))
write_tsv(new_gsea, file.path(out_dir, "gsea_new_control_contrasts.tsv"))

reuse_ref <- fread(file.path(bridge_root, "bulk_all_fixed_preranked_gsea.tsv"))[
  context == "bulk_LA_minus_NC" & time %in% times,
  .(feature_id, time, source_NES = NES, source_FDR = padj, source_ES = ES)
]
reuse_check <- merge(existing_la[, .(feature_id, time, reused_NES = NES, reused_FDR = padj, reused_ES = ES)],
                     reuse_ref, by = c("feature_id", "time"), all = TRUE)
reuse_check$NES_delta <- reuse_check$reused_NES - reuse_check$source_NES
reuse_check$FDR_delta <- reuse_check$reused_FDR - reuse_check$source_FDR
reuse_check$ES_delta <- reuse_check$reused_ES - reuse_check$source_ES
reuse_check$exact_reuse <- (is.na(reuse_check$reused_NES) == is.na(reuse_check$source_NES)) &
  (is.na(reuse_check$reused_NES) | abs(reuse_check$NES_delta) < 1e-12) &
  (is.na(reuse_check$reused_FDR) | abs(reuse_check$FDR_delta) < 1e-12) &
  (is.na(reuse_check$reused_ES) | abs(reuse_check$ES_delta) < 1e-12)
write_tsv(reuse_check, file.path(out_dir, "la_nc_reuse_check.tsv"))

score <- fread(file.path(bridge_root, "bulk_fixed_readout_ssgsea_contrasts.tsv"))[
  contrast %in% contrasts & time %in% times]
donor <- fread(file.path(bridge_root, "bulk_fixed_readout_donor_summary.tsv"))[
  contrast %in% contrasts & time %in% times]
stop_if(nrow(score) != 21L * 4L * 3L, paste("unexpected score rows:", nrow(score)))
stop_if(nrow(donor) != 21L * 4L * 3L, paste("unexpected donor score rows:", nrow(donor)))
write_tsv(score, file.path(out_dir, "program_score_three_contrasts.tsv"))
write_tsv(donor, file.path(out_dir, "program_donor_effects_three_contrasts.tsv"))

key_labels <- c("Hypoxia", "Glycolysis", "Unfolded protein response", "TCR signaling",
                "T-cell mediated cytotoxicity", "Effector-enriched (Goldrath effector vs memory, UP)",
                "Memory-enriched (Goldrath effector vs memory, DN)",
                "Antiviral mechanism by IFN-stimulated genes",
                "Class I MHC mediated antigen processing & presentation")
key_score <- score[feature_label %in% key_labels]
key_donor <- donor[feature_id %in% unique(key_score$feature_id)]
write_tsv(key_score, file.path(out_dir, "key_program_score_summary.tsv"))
write_tsv(key_donor, file.path(out_dir, "key_program_donor_summary.tsv"))

marker_groups <- list(
  activation = c("CD69", "IL2RA"),
  cytotoxic_effector = c("GZMB", "PRF1", "NKG7", "GZMA", "IFNG"),
  checkpoint_exhaustion_associated = c("PDCD1", "TIGIT", "LAG3", "CTLA4", "HAVCR2", "TOX"),
  memory_associated = c("CCR7", "IL7R", "TCF7", "SELL", "LEF1")
)
marker_manifest <- rbindlist(lapply(names(marker_groups), function(g) data.table(
  marker_group = g, gene = marker_groups[[g]], prespecified = TRUE
)))
marker_manifest$in_count <- marker_manifest$gene %in% rownames(count)
marker_manifest$n_samples_detected <- vapply(marker_manifest$gene, function(g) {
  if (!g %in% rownames(count)) return(0L)
  sum(count[g, ] > 0)
}, integer(1))
marker_manifest$detection_status <- ifelse(!marker_manifest$in_count, "NOT_IN_COUNT",
                                            ifelse(marker_manifest$n_samples_detected == 0L, "NO_COUNTS", "PRESENT"))
write_tsv(marker_manifest, file.path(out_dir, "cd8_marker_manifest.tsv"))

fac_ids <- factorial$sample_id
dge <- edgeR::DGEList(counts = count[, fac_ids, drop = FALSE], samples = factorial)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
expr <- rbindlist(lapply(marker_manifest$gene, function(g) {
  if (!g %in% rownames(logcpm)) return(data.table(gene = g))
  data.table(gene = g, sample_id = fac_ids, logCPM = as.numeric(logcpm[g, fac_ids]))
}), fill = TRUE)
expr <- merge(expr, factorial[, .(sample_id, Donor = as.character(Donor),
                                  Treatment = as.character(Treatment), Time = as.character(Time))],
              by = "sample_id", all.x = TRUE, sort = FALSE)
expr <- merge(expr, marker_manifest[, .(marker_group, gene, detection_status)],
              by = "gene", all.x = TRUE, sort = FALSE)
write_tsv(expr, file.path(out_dir, "cd8_marker_expression.tsv"))

marker_de <- rbindlist(lapply(contrasts, function(ct) {
  rbindlist(lapply(times, function(tm) {
    d <- load_de(ct, tm)$table
    x <- data.table(gene = marker_manifest$gene)
    x <- merge(x, d[, .(gene = gene_id, logFC, PValue, FDR)], by = "gene", all.x = TRUE, sort = FALSE)
    x$contrast <- ct
    x$time <- tm
    x[]
  }), fill = TRUE)
}), fill = TRUE)
marker_de <- merge(marker_de, marker_manifest[, .(marker_group, gene, detection_status, n_samples_detected)],
                   by = "gene", all.x = TRUE, sort = FALSE)
marker_de$de_status <- ifelse(marker_de$detection_status != "PRESENT", marker_de$detection_status,
                              ifelse(is.na(marker_de$logFC), "FILTERED_FROM_DE", "EVALUABLE"))
marker_de <- marker_de[order(match(marker_group, names(marker_groups)), match(gene, marker_manifest$gene),
                             match(contrast, contrasts), match(time, times))]
write_tsv(marker_de, file.path(out_dir, "cd8_marker_de.tsv"))

marker_donor_rows <- rbindlist(lapply(contrasts, function(ct) {
  pair <- switch(ct,
                 LA_minus_NC = c("LA", "NC"),
                 HCl_minus_NC = c("HCl", "NC"),
                 LA_minus_HCl = c("LA", "HCl"))
  rbindlist(lapply(times, function(tm) {
    left <- expr[Treatment == pair[1] & Time == tm & !is.na(logCPM),
                 .(Donor, gene, left_logCPM = logCPM)]
    right <- expr[Treatment == pair[2] & Time == tm & !is.na(logCPM),
                  .(Donor, gene, right_logCPM = logCPM)]
    x <- merge(left, right, by = c("Donor", "gene"), all = FALSE)
    x$contrast <- ct
    x$time <- tm
    x$donor_delta_logCPM <- x$left_logCPM - x$right_logCPM
    x[, .(contrast, time, Donor, gene, left_logCPM, right_logCPM, donor_delta_logCPM)]
  }), fill = TRUE)
}), fill = TRUE)
marker_donor_rows <- merge(marker_donor_rows,
                            marker_manifest[, .(marker_group, gene, detection_status)],
                            by = "gene", all.x = TRUE, sort = FALSE)
write_tsv(marker_donor_rows, file.path(out_dir, "cd8_marker_donor_effects.tsv"))
marker_donor_summary <- marker_donor_rows[, .(
  n_donors = uniqueN(Donor),
  mean_delta_logCPM = mean(donor_delta_logCPM),
  median_delta_logCPM = median(donor_delta_logCPM),
  sd_delta_logCPM = sd(donor_delta_logCPM),
  n_positive = sum(donor_delta_logCPM > 0),
  n_negative = sum(donor_delta_logCPM < 0),
  n_zero = sum(donor_delta_logCPM == 0)
), by = .(contrast, time, gene, marker_group, detection_status)]
marker_donor_summary$all_donors_same_direction <-
  marker_donor_summary$n_positive == marker_donor_summary$n_donors |
  marker_donor_summary$n_negative == marker_donor_summary$n_donors
write_tsv(marker_donor_summary, file.path(out_dir, "cd8_marker_donor_summary.tsv"))

qa <- data.table(
  check_id = c("observed_treatment_levels", "factorial_samples", "design_full_rank",
               "three_contrasts", "fixed_set_count", "new_gsea_rows", "la_nc_exact_reuse",
               "score_rows", "donor_rows", "marker_rows", "marker_expression_rows",
               "marker_donor_rows", "marker_donor_summary_rows",
               "no_unexpected_fourth_condition"),
  status = c(
    ifelse(identical(sort(unique(anno$Treatment)), sort(c("HCl", "LA", "NC", "Unsti"))), "PASS", "HOLD"),
    ifelse(nrow(factorial) == 36L, "PASS", "HOLD"),
    ifelse(qr(design)$rank == ncol(design), "PASS", "HOLD"),
    ifelse(all(contrasts %in% unique(score$contrast)), "PASS", "HOLD"),
    ifelse(length(all_sets) == 21L, "PASS", "HOLD"),
    ifelse(nrow(new_gsea) == length(all_sets) * 4L * 2L, "PASS", "HOLD"),
    ifelse(all(reuse_check$exact_reuse), "PASS", "HOLD"),
    ifelse(nrow(score) == 252L, "PASS", "HOLD"),
    ifelse(nrow(donor) == 252L && all(donor$n_donors == 3L), "PASS", "HOLD"),
    ifelse(nrow(marker_de) == nrow(marker_manifest) * 12L, "PASS", "HOLD"),
    ifelse(nrow(expr) > 0L, "PASS", "HOLD"),
    ifelse(nrow(marker_donor_rows) == nrow(marker_manifest) * 12L * 3L, "PASS", "HOLD"),
    ifelse(nrow(marker_donor_summary) == nrow(marker_manifest) * 12L, "PASS", "HOLD"),
    ifelse(identical(sort(unique(anno$Treatment)), sort(c("HCl", "LA", "NC", "Unsti"))), "PASS", "HOLD")
  ),
  detail = c(
    paste(sort(unique(anno$Treatment)), collapse = ";"), nrow(factorial),
    paste(qr(design)$rank, "of", ncol(design)), paste(contrasts, collapse = ";"),
    length(all_sets), nrow(new_gsea), sum(reuse_check$exact_reuse), nrow(score),
    nrow(donor), nrow(marker_de), nrow(expr), nrow(marker_donor_rows),
    nrow(marker_donor_summary), paste(sort(unique(anno$Treatment)), collapse = ";")
  )
)
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
overall <- if (all(qa$status == "PASS")) "PASS_CANDIDATE" else "HOLD_CANDIDATE"
writeLines(c(
  paste0("status=", overall),
  paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)),
  paste0("analysis_seed=", analysis_seed),
  paste0("factorial_samples=", nrow(factorial), ";donors=", uniqueN(factorial$Donor)),
  paste0("fixed_sets=", length(all_sets), ";new_control_gsea_rows=", nrow(new_gsea)),
  paste0("marker_genes=", uniqueN(marker_manifest$gene), ";marker_de_rows=", nrow(marker_de))
), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("STATUS\t", overall, "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
}
