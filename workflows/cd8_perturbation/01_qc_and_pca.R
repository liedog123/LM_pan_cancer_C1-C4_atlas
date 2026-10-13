#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 01_qc_and_pca.R --config PATH [--check-config]\n       Rscript 01_qc_and_pca.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("count_rdata","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["count_rdata"]
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

# ---- Retained complete source body; CLI/path substitutions only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

input_file <- unname(paths["count_rdata"])
out_dir <- unname(paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path, row.names = FALSE) {
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE, row.names = row.names,
    na = "NA", col.names = TRUE
  )
}

stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

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

