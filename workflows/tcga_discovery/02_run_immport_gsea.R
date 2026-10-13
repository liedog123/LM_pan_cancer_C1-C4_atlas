#!/usr/bin/env Rscript
# Continue TCGA discovery from separately obtained processed intermediates.

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 02_run_immport_gsea.R --config PATH [--check-config]\n",
      "       Rscript 02_run_immport_gsea.R --help\n\n",
      "DCF fields: deg_results_dir, immport_gmt, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation preprocesses gene sets and runs per-cancer fgsea.\n",
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
required_fields <- c("deg_results_dir", "immport_gmt", "output_dir")
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
deg_results_dir <- unname(paths["deg_results_dir"])
immport_gmt_file_path <- unname(paths["immport_gmt"])
out_dir <- unname(paths["output_dir"])
# Use the source filename pattern unchanged, including the unescaped period.
if (!dir.exists(deg_results_dir) || file.access(deg_results_dir, mode = 5L) != 0L) {
  stop("deg_results_dir must be a readable searchable directory.", call. = FALSE)
}
deg_metadata_files <- list.files(deg_results_dir, pattern = "_DEG_full_results.csv$", full.names = TRUE)
if (length(deg_metadata_files) == 0L || any(!file_test("-f", deg_metadata_files)) ||
    any(file.access(deg_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(deg_metadata_files)$size) | file.info(deg_metadata_files)$size <= 0)) {
  stop("Expected readable nonempty regular full-DEG files.", call. = FALSE)
}
deg_real_paths <- vapply(deg_metadata_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)
required_files <- c(immport_gmt_file_path, deg_real_paths)

# Protect inputs and existing outputs; checks create no directories.
protected_paths <- c(deg_results_dir, dirname(required_files), config_path)
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

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
})
output_dir <- out_dir
for (path in c(output_dir, file.path(output_dir, "individual_full_results"),
               file.path(output_dir, "individual_sig_results"))) {
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(path)) {
    stop("Cannot create output directory.", call. = FALSE)
  }
}

GSEA_PADJ_THRESHOLD <- 0.05          # 显著性阈值

message(paste("... 正在加载特定的 IMMPORT GMT 文件:", immport_gmt_file_path))

if (!file.exists(immport_gmt_file_path)) {
  stop(paste("错误: 无法找到指定的GMT文件:", immport_gmt_file_path))
}

# 读取指定的GMT文件
immport_pathways_df <- read.gmt(immport_gmt_file_path)
pathways_of_interest <- split(immport_pathways_df$gene, immport_pathways_df$term)

message(paste("... 成功从", basename(immport_gmt_file_path), "加载了", length(pathways_of_interest), "个通路。"))

# 限制基因集大小在 [5, 500]
pathway_sizes <- sapply(pathways_of_interest, length)
pathways_filtered_by_size <- pathways_of_interest[pathway_sizes >= 5 & pathway_sizes <= 500]
message(paste("... 按大小 [5, 500] 筛选后剩余:", length(pathways_filtered_by_size), "个通路。"))

# 去冗余筛选: Jaccard相似度 <= 0.5
genesets_to_check <- pathways_filtered_by_size
jaccard_matrix <- matrix(0, nrow=length(genesets_to_check), ncol=length(genesets_to_check),
                         dimnames=list(names(genesets_to_check), names(genesets_to_check)))

# 仅当通路数 > 1 时才计算 Jaccard
if (length(genesets_to_check) > 1) {
  for (i in 1:(length(genesets_to_check) - 1)) {
    for (j in (i + 1):length(genesets_to_check)) {
      intersect_size <- length(intersect(genesets_to_check[[i]], genesets_to_check[[j]]))
      union_size <- length(union(genesets_to_check[[i]], genesets_to_check[[j]]))
      jaccard_matrix[i, j] <- intersect_size / union_size
    }
  }
  highly_correlated <- which(jaccard_matrix > 0.5, arr.ind = TRUE)
  pathways_to_remove <- c()
  if(nrow(highly_correlated) > 0) {
    for(i in 1:nrow(highly_correlated)) {
      p1_name <- rownames(jaccard_matrix)[highly_correlated[i, 1]]
      p2_name <- colnames(jaccard_matrix)[highly_correlated[i, 2]]
      # 优先移除名字更长的（通常是子通路）
      if(nchar(p1_name) > nchar(p2_name)) pathways_to_remove <- c(pathways_to_remove, p1_name)
      else pathways_to_remove <- c(pathways_to_remove, p2_name)
    }
  }
  final_pathways <- genesets_to_check[!names(genesets_to_check) %in% unique(pathways_to_remove)]
} else {
  final_pathways <- genesets_to_check
  message("... 通路数 <= 1，跳过去冗余步骤。")
}

message(paste("... 去冗余后最终用于分析的基因集:", length(final_pathways), "个。"))


# ===============================
# --- 4. GSEA分析循环 ---
# ===============================
deg_files <- list.files(deg_results_dir, pattern = "_DEG_full_results.csv$", full.names = TRUE)
all_gsea_results <- list()
# 用于可视化的癌种列表
cancer_types_processed <- c()

for (deg_file in deg_files) {
  cancer_name <- sub("_DEG_full_results.csv$", "", basename(deg_file))
  message(paste("... 正在处理:", cancer_name))

  deg_results <- read.csv(deg_file)
  if (!"t" %in% colnames(deg_results)) {
    stop(paste("错误: 差异分析结果文件", basename(deg_file), "中缺少 't' 列。请确保差异分析脚本生成了t统计量。"))
  }
  ranks <- deframe(select(deg_results, gene, t))
  ranks <- ranks[!is.na(ranks) & !is.na(names(ranks)) & names(ranks) != ""]
  ranks <- sort(ranks, decreasing = TRUE)

  # 检查是否有通路可供分析
  if (length(final_pathways) == 0) {
    warning(paste("没有可供分析的通路，跳过GSEA步骤:", cancer_name))
    next
  }

  fgsea_res <- fgsea(pathways = final_pathways, stats = ranks, minSize = 10, maxSize = 500, nPermSimple = 100000)

  if (nrow(fgsea_res) > 0) {
    # 记录成功处理的癌种
    cancer_types_processed <- c(cancer_types_processed, cancer_name)

    fgsea_res_for_csv <- fgsea_res %>%
      mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ","))
    # 保存当前癌种的GSEA结果
    write.csv(fgsea_res_for_csv,
              file.path(output_dir, "individual_full_results", paste0(cancer_name, "_GSEA_full_results.csv")),
              row.names = FALSE)

    # 筛选并保存当前癌种的显著结果
    fgsea_res_sig <- fgsea_res %>% filter(padj < GSEA_PADJ_THRESHOLD) %>% arrange(desc(NES))
    if(nrow(fgsea_res_sig) > 0) {
      fgsea_res_sig_for_csv <- fgsea_res_sig %>%
        mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ","))

      write.csv(fgsea_res_sig_for_csv,
                file.path(output_dir, "individual_sig_results", paste0(cancer_name, "_GSEA_sig_results.csv")),
                row.names = FALSE)
    }

    all_gsea_results[[cancer_name]] <- fgsea_res
  }
}
saveRDS(all_gsea_results, file = file.path(output_dir, "all_gsea_results.rds"))

if (length(all_gsea_results) == 0) {
  stop("所有癌种均未产生GSEA结果，停止后续分析。")
}
