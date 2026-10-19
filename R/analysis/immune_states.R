# Purpose: immune states analysis stages.
# Configuration: config/states.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Identify acidification-associated immune pathway enrichment in TCGA.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/tcga_discovery/01_build_high_low_deg.R
stage_tcga_discovery_01_build_high_low_deg <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Identify acidification-associated immune pathway enrichment in TCGA.
# Task: tcga_immune_discovery; step: s01
# Inputs (DCF fields): processed_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task tcga_immune_discovery --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
processed_data <- unname(paths["processed_data"])
out_dir <- unname(paths["output_dir"])
required_files <- processed_data

# Protect inputs and existing outputs; checks create no directories.
protected_paths <- c(dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_026, environment())
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})
pancancer_data_list <- load_named_object(processed_data, "pancancer_data_list")
main_output_dir <- out_dir
full_results_dir <- file.path(main_output_dir, "full_degs_for_gsea")
sig_results_dir <- file.path(main_output_dir, "significant_degs")
for (path in c(out_dir, full_results_dir, sig_results_dir)) {
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(path)) {
    stop("Cannot create output directory.", call. = FALSE)
  }
}

cancer_types <- names(pancancer_data_list)
for (cancer in cancer_types) {
  message(paste("正在处理:", cancer))

  # 提取当前癌种的数据
  current_data <- pancancer_data_list[[cancer]]
  expr_matrix_full <- current_data$expression_data
  clinical_info_full <- current_data$clinical_data

  if (is.null(expr_matrix_full) || is.null(clinical_info_full) || nrow(clinical_info_full) == 0 || ncol(expr_matrix_full) == 0) {
    warning(paste("跳过", cancer, "，因为数据不完整。"))
    next
  }

  # 筛选出仅包含 "High" 和 "Low" 组的样本
  clinical_info <- clinical_info_full %>%
    filter(Group %in% c("High", "Low"))
  if(nrow(clinical_info) < 2 || length(unique(clinical_info$Group)) < 2) {
    warning(paste("跳过", cancer, "，因为过滤后High/Low组样本不足。"))
    next
  }
  expr_matrix <- expr_matrix_full[, clinical_info$sample]

  # 创建分组信息
  group_list <- factor(clinical_info$Group, levels = c("Low", "High"))
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- c("Low", "High")

  message("...正在进行差异分析")
  v <- voom(expr_matrix + 1, design, plot = FALSE)
  fit <- lmFit(v, design)
  contrast_matrix <- makeContrasts(High_vs_Low = High - Low, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  full_results <- topTable(fit2, number = Inf)
  full_results$gene <- rownames(full_results)

  output_path_full <- file.path(full_results_dir,
                                paste0(cancer, "_DEG_full_results.csv"))
  write.csv(full_results, file = output_path_full, row.names = FALSE)

  # 保存显著差异基因
  adj_p_value_threshold <- 0.05
  logfc_threshold <- 1.0

  significant_degs <- full_results %>%
    filter(adj.P.Val < adj_p_value_threshold, abs(logFC) > logfc_threshold)

  if (nrow(significant_degs) > 0) {
    output_path_sig <- file.path(sig_results_dir,
                                 paste0(cancer, "_DEG_significant.csv"))
    write.csv(significant_degs, file = output_path_sig, row.names = FALSE)
    message(paste("...处理完成！找到", nrow(significant_degs), "个显著差异基因。"))
  } else {
    message("...处理完成！未找到符合条件的显著差异基因。")
  }
}
}

# Purpose: Identify acidification-associated immune pathway enrichment in TCGA.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/tcga_discovery/02_run_immport_gsea.R
stage_tcga_discovery_02_run_immport_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Identify acidification-associated immune pathway enrichment in TCGA.
# Task: tcga_immune_discovery; step: s02
# Inputs (DCF fields): deg_results_dir,immport_gmt
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task tcga_immune_discovery --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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

GSEA_PADJ_THRESHOLD <- 0.05

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
}

# Purpose: Construct the adopted twelve immune pathway definitions from cross-cancer evidence.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/lm12_construction/01_build_votes_and_stable_genes.R
stage_lm12_construction_01_build_votes_and_stable_genes <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Construct the adopted twelve immune pathway definitions from cross-cancer evidence.
# Task: lm12_definition; step: s01
# Inputs (DCF fields): gsea_results_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lm12_definition --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
full_result_dir <- unname(paths["gsea_results_dir"])
out_dir <- unname(paths["output_dir"])
# Check only file metadata using the source's filename rule; resolve symlink targets.
if (!dir.exists(full_result_dir) || file.access(full_result_dir, mode = 5L) != 0L) {
  stop("gsea_results_dir must be a readable searchable directory.", call. = FALSE)
}
gsea_metadata_files <- list.files(full_result_dir,
  pattern = "_GSEA_full_results\\.csv$", full.names = TRUE)
if (length(gsea_metadata_files) == 0L ||
    any(!file_test("-f", gsea_metadata_files)) ||
    any(file.access(gsea_metadata_files, mode = 4L) != 0L) ||
    any(is.na(file.info(gsea_metadata_files)$size) | file.info(gsea_metadata_files)$size <= 0)) {
  stop("Expected at least one readable nonempty regular GSEA result CSV.", call. = FALSE)
}
gsea_real_paths <- vapply(gsea_metadata_files, normalizePath, character(1L),
                         winslash = "/", mustWork = TRUE)
required_files <- gsea_real_paths


# Protect inputs, configuration and existing results; check mode creates no directories.
protected_paths <- c(full_result_dir, dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
})
gene_set_dir <- file.path(out_dir, "gene_sets")
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
if (!dir.create(gene_set_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(gene_set_dir)) {
  stop("Cannot create gene_sets output directory.", call. = FALSE)
}

padj_cutoff <- 0.05
core_freq_cutoff <- 0.60
minor_freq_cutoff <- 0.40
leading_edge_freq_cutoff <- 0.70

parse_leading_edge <- function(x) {
  x <- ifelse(is.na(x), "", x)
  genes <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes[genes != ""]
}
`safe_filename` <- bind_stage_helper(shared_safe_filename_041, environment())

gsea_files <- list.files(
  full_result_dir,
  pattern = "_GSEA_full_results\\.csv$",
  full.names = TRUE
)
if (length(gsea_files) == 0) {
  stop("未找到每癌种 GSEA 结果文件。")
}

message("读取每癌种 GSEA 结果...")
gsea_df <- map_dfr(gsea_files, function(f) {
  cancer <- sub("_GSEA_full_results\\.csv$", "", basename(f))
  read.csv(f, check.names = FALSE) %>%
    mutate(cancer = cancer)
})

n_cancers <- length(unique(gsea_df$cancer))
message("GSEA 癌种数: ", n_cancers)

select_pathways_one_direction <- function(direction_label, direction_filter) {
  sig_df <- gsea_df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    filter(direction_filter(NES)) %>%
    mutate(direction = direction_label)

  vote_df <- sig_df %>%
    distinct(direction, pathway, cancer, .keep_all = TRUE) %>%
    group_by(direction, pathway) %>%
    summarise(
      vote_count = n_distinct(cancer),
      vote_freq = vote_count / n_cancers,
      mean_NES = mean(NES, na.rm = TRUE),
      median_NES = median(NES, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      recurrence_class = case_when(
        vote_freq >= core_freq_cutoff ~ "Core",
        vote_freq >= minor_freq_cutoff ~ "Minor",
        TRUE ~ "Sporadic"
      )
    ) %>%
    arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)

  selected_df <- vote_df %>%
    filter(vote_freq >= minor_freq_cutoff)

  list(sig_df = sig_df, vote_df = vote_df, selected_df = selected_df)
}

high_res <- select_pathways_one_direction("High_Acid_Positive", function(x) x > 0)
low_res <- select_pathways_one_direction("Low_Acid_Negative", function(x) x < 0)

all_vote_df <- bind_rows(high_res$vote_df, low_res$vote_df) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)
selected_pathways <- bind_rows(high_res$selected_df, low_res$selected_df) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)

write_tsv(all_vote_df, file.path(out_dir, "all_pathway_votes_original_like.tsv"))
write_tsv(selected_pathways, file.path(out_dir, "selected_core_minor_pathways_original_like.tsv"))

message("构建 pathway-level stable leading-edge gene set...")
sig_all <- bind_rows(high_res$sig_df, low_res$sig_df)
gene_set_rows <- list()
gene_set_summary <- selected_pathways %>%
  mutate(
    stable_gene_count = 0L,
    leading_edge_threshold_count = ceiling(vote_count * leading_edge_freq_cutoff)
  )

for (i in seq_len(nrow(selected_pathways))) {
  row <- selected_pathways[i, ]
  sub_df <- sig_all %>%
    filter(direction == row$direction, pathway == row$pathway) %>%
    distinct(cancer, .keep_all = TRUE)

  leading_genes <- unlist(lapply(sub_df$leadingEdge, parse_leading_edge), use.names = FALSE)
  gene_counts <- sort(table(leading_genes), decreasing = TRUE)
  threshold_count <- ceiling(row$vote_count * leading_edge_freq_cutoff)
  stable_genes <- names(gene_counts)[as.integer(gene_counts) >= threshold_count]

  gene_set_summary$stable_gene_count[i] <- length(stable_genes)

  if (length(stable_genes) > 0) {
    gene_set_rows[[length(gene_set_rows) + 1]] <- data.frame(
      direction = row$direction,
      pathway = row$pathway,
      recurrence_class = row$recurrence_class,
      gene = stable_genes,
      leading_edge_count = as.integer(gene_counts[stable_genes]),
      vote_count = row$vote_count,
      leading_edge_threshold_count = threshold_count,
      stringsAsFactors = FALSE
    )

    prefix <- paste(
      row$direction,
      row$recurrence_class,
      safe_filename(row$pathway),
      sep = "__"
    )
    writeLines(unique(stable_genes), file.path(gene_set_dir, paste0(prefix, ".txt")))
  }
}

gene_set_long <- bind_rows(gene_set_rows) %>%
  arrange(direction, desc(vote_count), pathway, desc(leading_edge_count), gene)

write_tsv(gene_set_summary, file.path(out_dir, "original_like_gene_set_summary.tsv"))
write_tsv(gene_set_long, file.path(out_dir, "original_like_gene_sets_long.tsv"))

summary_lines <- c(
  "Original-like pathway audit",
  paste0("Cancer types: ", n_cancers),
  paste0("padj cutoff: ", padj_cutoff),
  paste0("Recurrence rule: Core >= ", core_freq_cutoff, ", Minor >= ", minor_freq_cutoff),
  paste0("Stable leading-edge rule: freq >= ", leading_edge_freq_cutoff),
  paste0(
    "Selected pathways: ", nrow(selected_pathways),
    " (High_Acid_Positive=", sum(selected_pathways$direction == "High_Acid_Positive"),
    ", Low_Acid_Negative=", sum(selected_pathways$direction == "Low_Acid_Negative"), ")"
  ),
  paste0("Non-empty stable gene sets: ", n_distinct(gene_set_long$pathway)),
  paste0("Output directory: ", out_dir)
)

writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
}

# Purpose: Construct the adopted twelve immune pathway definitions from cross-cancer evidence.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/lm12_construction/02_finalize_lm12_gene_sets.R
stage_lm12_construction_02_finalize_lm12_gene_sets <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Construct the adopted twelve immune pathway definitions from cross-cancer evidence.
# Task: lm12_definition; step: s02
# Inputs (DCF fields): all_votes,gene_summary,gene_long,candidate_terms,family_dictionary
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task lm12_definition --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
all_vote_file <- unname(paths["all_votes"])
gene_summary_file <- unname(paths["gene_summary"])
gene_long_file <- unname(paths["gene_long"])
candidate_terms_file <- unname(paths["candidate_terms"])
family_dict_file <- unname(paths["family_dictionary"])
out_dir <- unname(paths["output_dir"])
required_files <- c(all_vote_file, gene_summary_file, gene_long_file, candidate_terms_file, family_dict_file)

# Protect inputs, configuration and existing results; check mode creates no directories.
protected_paths <- c(dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})
gene_set_dir <- file.path(out_dir, "gene_sets")
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
if (!dir.create(gene_set_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(gene_set_dir)) {
  stop("Cannot create gene_sets output directory.", call. = FALSE)
}

required_files <- c(
  all_vote_file,
  gene_summary_file,
  gene_long_file,
  candidate_terms_file,
  family_dict_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("缺失输入文件:\n", paste(missing_files, collapse = "\n"))
}

strict_recurrent_cutoff <- 0.80
strict_core_cutoff <- 0.90
original_like_core_cutoff <- 0.60
original_like_minor_cutoff <- 0.40
`safe_filename` <- bind_stage_helper(shared_safe_filename_041, environment())

annotation_df <- data.frame(
  pathway = c(
    "Antiviral mechanism by IFN-stimulated genes",
    "MHC class II antigen presentation",
    "Fcgamma receptor (FCGR) dependent phagocytosis",
    "Interleukin-1 family signaling",
    "Class I MHC mediated antigen processing & presentation",
    "Neutrophil degranulation",
    "Interleukin-17 signaling",
    "Cytosolic sensors of pathogen-associated DNA",
    "Toll-like Receptor Cascades",
    "TCR signaling",
    "C-type lectin receptors (CLRs)",
    "Antimicrobial peptides"
  ),
  feature_id = c(
    "LM12_IFN_ISG",
    "LM12_MHC_II",
    "LM12_FCGR_Phagocytosis",
    "LM12_IL1",
    "LM12_MHC_I",
    "LM12_Neutrophil_Degranulation",
    "LM12_IL17",
    "LM12_Cytosolic_DNA_Sensing",
    "LM12_TLR_Cascades",
    "LM12_TCR_Signaling",
    "LM12_CLR",
    "LM12_Antimicrobial_Peptides"
  ),
  short_label = c(
    "IFN_ISG",
    "MHC_II",
    "FCGR_Phagocytosis",
    "IL1",
    "MHC_I",
    "Neutrophil_Degranulation",
    "IL17",
    "Cytosolic_DNA_Sensing",
    "TLR_Cascades",
    "TCR_Signaling",
    "CLR",
    "Antimicrobial_Peptides"
  ),
  functional_group = c(
    "ifn_antigen_presentation",
    "ifn_antigen_presentation",
    "myeloid_inflammation",
    "myeloid_inflammation",
    "ifn_antigen_presentation",
    "myeloid_inflammation",
    "myeloid_inflammation",
    "innate_sensing",
    "innate_sensing",
    "t_cell_signaling",
    "innate_sensing",
    "low_acid_barrier_defense"
  ),
  display_order = seq_len(12),
  stringsAsFactors = FALSE
)

expected_pathways <- annotation_df$pathway

all_votes <- read_tsv(all_vote_file, show_col_types = FALSE)
gene_summary <- read_tsv(gene_summary_file, show_col_types = FALSE) %>%
  select(direction, pathway, stable_gene_count, leading_edge_threshold_count)
gene_long <- read_tsv(gene_long_file, show_col_types = FALSE)
candidate_terms <- read_tsv(candidate_terms_file, show_col_types = FALSE)
family_dict <- read_tsv(family_dict_file, show_col_types = FALSE)

primary_candidates <- candidate_terms %>%
  filter(analysis_tier == "primary") %>%
  distinct(candidate_set, scope_layer, family_id, local_term, match_type, note)

candidate_votes <- all_votes %>%
  inner_join(primary_candidates, by = c("pathway" = "local_term")) %>%
  left_join(gene_summary, by = c("direction", "pathway")) %>%
  mutate(
    stable_gene_count = replace_na(stable_gene_count, 0L),
    original_like_class = case_when(
      vote_freq >= original_like_core_cutoff ~ "Core",
      vote_freq >= original_like_minor_cutoff ~ "Minor",
      TRUE ~ "Sporadic"
    ),
    strict_class = case_when(
      vote_freq >= strict_core_cutoff ~ "Core",
      vote_freq >= strict_recurrent_cutoff ~ "Minor",
      TRUE ~ "Below_Strict"
    )
  ) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)

write_tsv(candidate_votes, file.path(out_dir, "predefined_primary_candidate_all_votes.tsv"))

selected_original_like <- candidate_votes %>%
  filter(vote_freq >= original_like_minor_cutoff) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), pathway)
write_tsv(
  selected_original_like,
  file.path(out_dir, "selected_primary_candidates_original_like_60_40.tsv")
)

selected_strict <- candidate_votes %>%
  filter(vote_freq >= strict_recurrent_cutoff) %>%
  arrange(direction, desc(vote_freq), desc(abs(mean_NES)), desc(stable_gene_count), pathway)
write_tsv(selected_strict, file.path(out_dir, "selected_primary_candidates_strict_80.tsv"))

unexpected_pathways <- setdiff(selected_strict$pathway, expected_pathways)
missing_expected_pathways <- setdiff(expected_pathways, selected_strict$pathway)
if (length(unexpected_pathways) > 0 || length(missing_expected_pathways) > 0) {
  stop(
    paste0(
      "严格 0.80 结果与预期 12 条 pathway 不一致。\nUnexpected: ",
      paste(unexpected_pathways, collapse = "; "),
      "\nMissing: ",
      paste(missing_expected_pathways, collapse = "; ")
    )
  )
}

lm_imm_12_manifest <- selected_strict %>%
  left_join(annotation_df, by = "pathway") %>%
  left_join(
    family_dict %>%
      distinct(
        scope_layer,
        family_id,
        family_label = fig1_or_project_basis,
        family_description = note
      ),
    by = c("scope_layer", "family_id")
  ) %>%
  mutate(
    selected_rule = "strict_primary_candidate_vote_freq_ge_0.80",
    feature_class = strict_class,
    direction_label = direction
  ) %>%
  arrange(display_order) %>%
  select(
    display_order,
    feature_id,
    short_label,
    pathway,
    direction_label,
    feature_class,
    vote_count,
    vote_freq,
    mean_NES,
    median_NES,
    stable_gene_count,
    leading_edge_threshold_count,
    scope_layer,
    family_id,
    family_label,
    functional_group,
    selected_rule
  )

lm_imm_12_gene_sets_long <- gene_long %>%
  semi_join(lm_imm_12_manifest, by = c("pathway")) %>%
  left_join(
    lm_imm_12_manifest %>%
      select(feature_id, short_label, pathway, functional_group, direction_label, feature_class),
    by = "pathway"
  ) %>%
  arrange(match(feature_id, lm_imm_12_manifest$feature_id), desc(leading_edge_count), gene) %>%
  select(
    feature_id,
    short_label,
    pathway,
    functional_group,
    direction = direction_label,
    feature_class,
    gene,
    leading_edge_count,
    vote_count,
    leading_edge_threshold_count
  )

for (feature_name in lm_imm_12_manifest$feature_id) {
  sub_df <- lm_imm_12_gene_sets_long %>%
    filter(feature_id == feature_name)
  pathway_name <- unique(sub_df$pathway)
  filename <- paste0(feature_name, "__", safe_filename(pathway_name), ".txt")
  writeLines(unique(sub_df$gene), file.path(gene_set_dir, filename))
}

function_annotation_df <- lm_imm_12_manifest %>%
  mutate(
    functional_group_label = case_when(
      functional_group == "ifn_antigen_presentation" ~ "IFN_and_antigen_presentation",
      functional_group == "t_cell_signaling" ~ "T_cell_signaling",
      functional_group == "innate_sensing" ~ "Innate_sensing",
      functional_group == "myeloid_inflammation" ~ "Myeloid_inflammation",
      functional_group == "low_acid_barrier_defense" ~ "Low_acid_barrier_defense",
      TRUE ~ functional_group
    )
  ) %>%
  select(
    feature_id,
    pathway,
    functional_group,
    functional_group_label,
    family_id,
    family_label,
    direction_label
  )

write_tsv(lm_imm_12_manifest, file.path(out_dir, "lm_imm_12_manifest.tsv"))
write_tsv(lm_imm_12_manifest, file.path(out_dir, "lm_imm_12_pathway_vote_summary.tsv"))
write_tsv(lm_imm_12_gene_sets_long, file.path(out_dir, "lm_imm_12_gene_sets_long.tsv"))
write_tsv(function_annotation_df, file.path(out_dir, "lm_imm_12_functional_annotation.tsv"))

summary_lines <- c(
  "LM-Imm 12 pathway-level gene set finalization",
  paste0("Primary candidate term count: ", n_distinct(primary_candidates$local_term)),
  paste0("Strict recurrence cutoff: ", strict_recurrent_cutoff),
  paste0("Strict Core cutoff: ", strict_core_cutoff),
  paste0(
    "Strict selected features: ", nrow(lm_imm_12_manifest),
    " (High_Acid_Positive=", sum(lm_imm_12_manifest$direction_label == "High_Acid_Positive"),
    ", Low_Acid_Negative=", sum(lm_imm_12_manifest$direction_label == "Low_Acid_Negative"), ")"
  ),
  paste0(
    "Functional groups: ",
    paste(unique(function_annotation_df$functional_group_label), collapse = "; ")
  ),
  paste0("Output directory: ", out_dir)
)

writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
}

# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/shared_definitions/01_assign_ig_and_c1c4_from_processed_input.R
stage_shared_definitions_01_assign_ig_and_c1c4_from_processed_input <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Task: c1c4_state_assignment; step: s01
# Inputs (DCF fields): processed_data,lm12_manifest,lm12_genes
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_state_assignment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
is_test <- "--test" %in% args
cancers_arg <- grep("^--cancers=", args, value = TRUE)
input_data_file <- unname(paths["processed_data"])
manifest_file <- unname(paths["lm12_manifest"])
gene_long_file <- unname(paths["lm12_genes"])
result_dir <- unname(paths["output_dir"])
required_files <- c(input_data_file, manifest_file, gene_long_file)

# 防止输出覆盖输入、配置或既有结果；检查模式不会创建目录。
protected_paths <- c(dirname(required_files), config_path)
if (result_dir == "/" || result_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, result_dir) || is_within(result_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(result_dir) && !dir.exists(result_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(result_dir)) {
  if (file.access(result_dir, mode = 7L) != 0L ||
      length(list.files(result_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All three inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- result_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(GSVA)
})
if (!dir.create(result_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(result_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

message("=== 6.0 AG3 x IG2 -> C1-C4 subtype build ===")
message("Input data: ", input_data_file)
message("LM12 manifest: ", manifest_file)
message("Output dir: ", result_dir)

load(input_data_file)
if (!exists("pancancer_data_list")) {
  stop("Object `pancancer_data_list` not found in ", input_data_file)
}

manifest <- read_tsv(manifest_file, show_col_types = FALSE) %>%
  arrange(display_order)
gene_long <- read_tsv(gene_long_file, show_col_types = FALSE)

feature_metadata <- manifest %>%
  transmute(
    display_order,
    feature_id,
    short_label,
    pathway,
    functional_group,
    direction_label,
    stable_gene_count
  )

feature_ids <- feature_metadata$feature_id
immune_high_features <- setdiff(feature_ids, "LM12_Antimicrobial_Peptides")
row_label_map <- setNames(feature_metadata$short_label, feature_metadata$feature_id)

gene_sets <- lapply(feature_ids, function(fid) {
  gene_long %>%
    filter(feature_id == fid) %>%
    pull(gene) %>%
    unique() %>%
    sort()
})
names(gene_sets) <- feature_ids

selected_cancers <- names(pancancer_data_list)
if (length(cancers_arg) > 0) {
  selected_cancers <- strsplit(sub("^--cancers=", "", cancers_arg[1]), ",")[[1]]
}
if (is_test && length(cancers_arg) == 0) {
  selected_cancers <- head(selected_cancers, 3)
}
selected_cancers <- intersect(selected_cancers, names(pancancer_data_list))
if (length(selected_cancers) == 0) {
  stop("没有可处理的癌种。")
}

message("Cancers to process: ", paste(selected_cancers, collapse = ", "))

ag_levels <- c("AG1", "AG2", "AG3")
ig_levels <- c("IG1", "IG2")
c_levels <- c("C1", "C2", "C3", "C4")
agig_levels <- as.vector(outer(ag_levels, ig_levels, paste, sep = "_"))

ag_map <- c(High = "AG1", Medium = "AG2", Low = "AG3")
ag_label_map <- c(
  AG1 = "acid_high",
  AG2 = "acid_intermediate",
  AG3 = "acid_low"
)
ig_label_map <- c(
  IG1 = "lm12_inflammatory_remodeling_high",
  IG2 = "lm12_inflammatory_remodeling_low"
)
c_map <- c(
  AG1_IG1 = "C1",
  AG1_IG2 = "C2",
  AG2_IG1 = "C3",
  AG2_IG2 = "C2",
  AG3_IG1 = "C3",
  AG3_IG2 = "C4"
)
c_label_map <- c(
  C1 = "high_acid__high_inflammatory_remodeling",
  C2 = "high_or_intermediate_acid__low_inflammatory_remodeling",
  C3 = "intermediate_or_low_acid__high_inflammatory_remodeling",
  C4 = "low_acid__low_inflammatory_remodeling"
)

score_list <- list()
coverage_list <- list()
ig_label_list <- list()
ig_manifest_list <- list()
ig_center_list <- list()
sample_count_records <- list()

for (cancer in selected_cancers) {
  message("\n--- Processing ", cancer, " ---")
  current_data <- pancancer_data_list[[cancer]]
  expr_matrix <- current_data$expression_data
  acid_scores <- current_data$acid_scores
  clinical_data <- current_data$clinical_data

  if (is.null(expr_matrix) || is.null(acid_scores) || is.null(clinical_data)) {
    stop("癌种 ", cancer, " 缺失 expression_data、acid_scores 或 clinical_data。")
  }
  if (!"sample" %in% colnames(clinical_data)) {
    stop("癌种 ", cancer, " 的 clinical_data 缺失 sample 列。")
  }
  if (!"Group" %in% colnames(clinical_data)) {
    stop("癌种 ", cancer, " 的 clinical_data 缺失 Group 列，无法定义 AG。")
  }

  common_samples <- Reduce(
    intersect,
    list(colnames(expr_matrix), colnames(acid_scores), clinical_data$sample)
  )
  if (length(common_samples) < 10) {
    stop("癌种 ", cancer, " 共同样本数 < 10。")
  }

  clinical_use <- clinical_data %>%
    filter(sample %in% common_samples) %>%
    distinct(sample, .keep_all = TRUE)
  ordered_samples <- clinical_use$sample

  expr_use <- expr_matrix[, ordered_samples, drop = FALSE]
  acid_values <- as.numeric(acid_scores[1, ordered_samples])
  names(acid_values) <- ordered_samples

  if (any(is.na(acid_values))) {
    stop("癌种 ", cancer, " 的 Acid_Score 存在 NA。")
  }

  expr_genes <- rownames(expr_use)
  coverage_df <- bind_rows(lapply(feature_ids, function(feature_id) {
    genes_total <- gene_sets[[feature_id]]
    genes_present <- intersect(genes_total, expr_genes)
    tibble(
      Cancer = cancer,
      feature_id = feature_id,
      short_label = row_label_map[[feature_id]],
      genes_total = length(genes_total),
      genes_present = length(genes_present),
      coverage = length(genes_present) / length(genes_total)
    )
  }))
  coverage_list[[cancer]] <- coverage_df

  gene_sets_present <- lapply(feature_ids, function(feature_id) {
    intersect(gene_sets[[feature_id]], expr_genes)
  })
  names(gene_sets_present) <- feature_ids

  zero_coverage_features <- names(gene_sets_present)[lengths(gene_sets_present) == 0]
  if (length(zero_coverage_features) > 0) {
    stop(
      "癌种 ", cancer, " 存在 gene set 完全无覆盖: ",
      paste(zero_coverage_features, collapse = ", ")
    )
  }

  expr_log2 <- log2(expr_use + 1)
  ssgsea_param <- ssgseaParam(
    exprData = expr_log2,
    geneSets = gene_sets_present,
    normalize = TRUE
  )
  score_matrix <- gsva(ssgsea_param, verbose = FALSE)

  score_df <- as.data.frame(t(score_matrix)) %>%
    rownames_to_column("Sample") %>%
    mutate(
      Cancer = cancer,
      Acid_Score = acid_values[Sample],
      Acid_Group = clinical_use$Group[match(Sample, clinical_use$sample)],
      .after = Sample
    )

  unknown_groups <- setdiff(unique(score_df$Acid_Group), names(ag_map))
  if (length(unknown_groups) > 0) {
    stop("癌种 ", cancer, " 出现未知 Acid_Group: ", paste(unknown_groups, collapse = ", "))
  }

  score_list[[cancer]] <- score_df
  sample_count_records[[cancer]] <- tibble(Cancer = cancer, n_samples = nrow(score_df))

  immune_axis <- rowMeans(score_df[, immune_high_features, drop = FALSE], na.rm = TRUE)
  cluster_input <- as.matrix(score_df[, feature_ids, drop = FALSE])
  rownames(cluster_input) <- score_df$Sample

  dist_matrix <- dist(cluster_input, method = "euclidean")
  hclust_res <- hclust(dist_matrix, method = "ward.D")
  raw_labels <- cutree(hclust_res, k = 2)

  tmp_label <- tibble(
    Sample = names(raw_labels),
    Cancer = cancer,
    IG2_raw_cluster = as.integer(raw_labels),
    Acid_Score = score_df$Acid_Score[match(names(raw_labels), score_df$Sample)],
    immune_inflammatory_axis = immune_axis[match(names(raw_labels), score_df$Sample)]
  )

  cluster_stats <- tmp_label %>%
    group_by(IG2_raw_cluster) %>%
    summarise(
      n_samples = n(),
      mean_immune_inflammatory_axis = mean(immune_inflammatory_axis, na.rm = TRUE),
      median_immune_inflammatory_axis = median(immune_inflammatory_axis, na.rm = TRUE),
      mean_acid_score = mean(Acid_Score, na.rm = TRUE),
      median_acid_score = median(Acid_Score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_immune_inflammatory_axis)) %>%
    mutate(
      IG = factor(c("IG1", "IG2"), levels = ig_levels),
      IG_Label = unname(ig_label_map[as.character(IG)]),
      immune_rank_within_cancer = c(1L, 2L),
      Cancer = cancer,
      .before = 1
    )

  ig_label_list[[cancer]] <- tmp_label %>%
    select(Sample, Cancer, IG2_raw_cluster) %>%
    left_join(
      cluster_stats %>%
        select(Cancer, IG2_raw_cluster, IG, IG_Label, immune_rank_within_cancer),
      by = c("Cancer", "IG2_raw_cluster")
    )

  ig_manifest_list[[cancer]] <- cluster_stats %>%
    select(Cancer, IG, IG_Label, IG2_raw_cluster, immune_rank_within_cancer, everything())

  ig_center_list[[cancer]] <- score_df %>%
    select(Sample, all_of(feature_ids)) %>%
    inner_join(ig_label_list[[cancer]] %>% select(Sample, IG), by = "Sample") %>%
    group_by(IG) %>%
    summarise(across(all_of(feature_ids), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(Cancer = cancer, .before = 1)

  message("  samples=", nrow(score_df), "; IG cluster sizes=", paste(cluster_stats$n_samples, collapse = "/"))
}

score_df_all <- bind_rows(score_list) %>%
  arrange(Cancer, Sample)
coverage_summary <- bind_rows(coverage_list) %>%
  left_join(feature_metadata, by = "feature_id") %>%
  select(
    Cancer,
    display_order,
    feature_id,
    short_label = short_label.x,
    pathway,
    functional_group,
    genes_total,
    genes_present,
    coverage
  ) %>%
  arrange(Cancer, display_order)

ig_labels <- bind_rows(ig_label_list)
ig_manifest <- bind_rows(ig_manifest_list) %>%
  arrange(Cancer, immune_rank_within_cancer)
ig_centers <- bind_rows(ig_center_list) %>%
  arrange(Cancer, IG)
sample_count_summary <- bind_rows(sample_count_records) %>%
  arrange(Cancer)

label_df <- score_df_all %>%
  left_join(ig_labels, by = c("Sample", "Cancer")) %>%
  mutate(
    AG = unname(ag_map[Acid_Group]),
    AG_Label = unname(ag_label_map[AG]),
    AGIG = paste(AG, as.character(IG), sep = "_"),
    C_subtype = unname(c_map[AGIG]),
    C_label = unname(c_label_map[C_subtype]),
    immune_inflammatory_axis = rowMeans(across(all_of(immune_high_features)), na.rm = TRUE)
  )

if (any(is.na(label_df$AG))) {
  stop("AG 映射失败，存在未分配样本。")
}
if (any(is.na(label_df$IG))) {
  stop("IG 映射失败，存在未分配样本。")
}
if (any(is.na(label_df$C_subtype))) {
  stop("C_subtype 映射失败，存在未分配 AGIG。")
}

label_df <- label_df %>%
  mutate(
    AG = factor(AG, levels = ag_levels),
    IG = factor(IG, levels = ig_levels),
    AGIG = factor(AGIG, levels = agig_levels),
    C_subtype = factor(C_subtype, levels = c_levels)
  ) %>%
  arrange(Cancer, C_subtype, AG, IG, Sample)

atomic_label_df <- label_df %>%
  select(
    Sample,
    Cancer,
    Acid_Score,
    Acid_Group,
    AG,
    AG_Label,
    IG2_raw_cluster,
    IG,
    IG_Label,
    immune_rank_within_cancer,
    immune_inflammatory_axis,
    AGIG,
    C_subtype,
    C_label,
    all_of(feature_ids)
  )

counts_overall <- atomic_label_df %>%
  count(C_subtype, AGIG, AG, IG, C_label, name = "n_samples") %>%
  arrange(C_subtype, AG, IG)

counts_by_cancer_long <- atomic_label_df %>%
  count(Cancer, C_subtype, name = "n_samples") %>%
  complete(
    Cancer = sort(unique(atomic_label_df$Cancer)),
    C_subtype = factor(c_levels, levels = c_levels),
    fill = list(n_samples = 0L)
  ) %>%
  left_join(
    tibble(C_subtype = factor(c_levels, levels = c_levels), C_label = unname(c_label_map[c_levels])),
    by = "C_subtype"
  ) %>%
  arrange(Cancer, C_subtype)

counts_by_cancer_wide <- counts_by_cancer_long %>%
  select(Cancer, C_subtype, n_samples) %>%
  pivot_wider(names_from = C_subtype, values_from = n_samples, values_fill = 0) %>%
  mutate(
    total_samples = rowSums(across(all_of(c_levels))),
    all_groups_present = if_all(all_of(c_levels), ~ .x > 0),
    all_groups_n_ge_10 = if_all(all_of(c_levels), ~ .x >= 10),
    all_groups_n_ge_20 = if_all(all_of(c_levels), ~ .x >= 20),
    min_group_n = pmin(C1, C2, C3, C4)
  ) %>%
  arrange(Cancer)

agig_counts_by_cancer <- atomic_label_df %>%
  count(Cancer, AGIG, AG, IG, name = "n_samples") %>%
  complete(
    Cancer = sort(unique(atomic_label_df$Cancer)),
    AGIG = factor(agig_levels, levels = agig_levels),
    fill = list(n_samples = 0L)
  ) %>%
  mutate(
    AG = factor(sub("_IG[12]$", "", as.character(AGIG)), levels = ag_levels),
    IG = factor(sub("^AG[123]_", "", as.character(AGIG)), levels = ig_levels)
  ) %>%
  arrange(Cancer, AG, IG)

write_tsv(feature_metadata, file.path(result_dir, "lm12_feature_metadata.tsv"))
write_tsv(coverage_summary, file.path(result_dir, "lm12_gene_coverage.tsv"))
write_tsv(score_df_all, file.path(result_dir, "pancancer_lm12_scores_by_sample.tsv.gz"))
write_tsv(ig_manifest, file.path(result_dir, "ig2_cluster_manifest_by_cancer.tsv"))
write_tsv(ig_centers, file.path(result_dir, "ig2_cluster_feature_medians_by_cancer.tsv"))
write_csv(atomic_label_df, file.path(result_dir, "ag3ig2_c1c4_sample_labels.csv"))
write_tsv(counts_overall, file.path(result_dir, "ag3ig2_c1c4_counts_overall.tsv"))
write_tsv(counts_by_cancer_long, file.path(result_dir, "ag3ig2_c1c4_counts_by_cancer_long.tsv"))
write_tsv(counts_by_cancer_wide, file.path(result_dir, "ag3ig2_c1c4_counts_by_cancer_wide.tsv"))
write_tsv(agig_counts_by_cancer, file.path(result_dir, "ag3ig2_atomic_counts_by_cancer.tsv"))
write_tsv(sample_count_summary, file.path(result_dir, "sample_count_summary_by_cancer.tsv"))

save(
  score_list,
  score_df_all,
  feature_metadata,
  gene_sets,
  coverage_summary,
  ig_labels,
  ig_manifest,
  ig_centers,
  atomic_label_df,
  counts_overall,
  counts_by_cancer_long,
  counts_by_cancer_wide,
  agig_counts_by_cancer,
  file = file.path(result_dir, "ag3ig2_c1c4_subtype_build_objects.Rdata")
)

summary_lines <- c(
  paste("analysis", "AG3xIG2_C1C4_formal_subtype_build", sep = "\t"),
  paste("created_at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sep = "\t"),
  paste("test_mode", is_test, sep = "\t"),
  paste("n_cancers", length(selected_cancers), sep = "\t"),
  paste("n_samples", nrow(atomic_label_df), sep = "\t"),
  paste("n_features", length(feature_ids), sep = "\t"),
  paste("feature_ids", paste(feature_ids, collapse = ";"), sep = "\t"),
  paste("n_cancers_all_groups_present", sum(counts_by_cancer_wide$all_groups_present), sep = "\t"),
  paste("n_cancers_all_groups_n_ge_10", sum(counts_by_cancer_wide$all_groups_n_ge_10), sep = "\t"),
  paste("n_cancers_all_groups_n_ge_20", sum(counts_by_cancer_wide$all_groups_n_ge_20), sep = "\t"),
  paste("counts_overall", paste(paste(counts_overall$C_subtype, counts_overall$AGIG, counts_overall$n_samples, sep = ":"), collapse = ";"), sep = "\t")
)
writeLines(summary_lines, file.path(result_dir, "run_summary.txt"))

message("\nSaved outputs to: ", result_dir)
message("=== 01_build_ag3ig2_c1c4_subtypes.R complete ===")
}

# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/shared_definitions/02_validate_state_axes_from_labels.R
stage_shared_definitions_02_validate_state_axes_from_labels <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Task: c1c4_state_assignment; step: s02
# Inputs (DCF fields): assignment_labels,feature_metadata
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_state_assignment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
label_file <- unname(paths["assignment_labels"])
feature_metadata_file <- unname(paths["feature_metadata"])
out_dir <- unname(paths["output_dir"])
required_files <- c(label_file, feature_metadata_file)

# 防止输出覆盖输入、配置或既有结果；检查模式不会创建目录。
protected_paths <- c(dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

c_levels <- c("C1", "C2", "C3", "C4")
comparisons <- tibble(
  comparison = c("C1_vs_C2", "C3_vs_C4", "C2_vs_C3", "C1_vs_C4", "C1_vs_C3", "C2_vs_C4"),
  group_1 = c("C1", "C3", "C2", "C1", "C1", "C2"),
  group_2 = c("C2", "C4", "C3", "C4", "C3", "C4"),
  comparison_class = c(
    "high_acid_inflammatory_axis",
    "low_acid_inflammatory_axis",
    "decoupled_intermediate_states",
    "extreme_endpoint",
    "reference",
    "reference"
  )
)

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

safe_quantile <- function(x, p) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

kruskal_one <- function(df, score_col) {
  df <- df %>% filter(!is.na(.data[[score_col]]), !is.na(C_subtype))
  n <- nrow(df)
  n_groups <- n_distinct(df$C_subtype)
  min_group_n <- if (n_groups > 0) min(as.integer(table(df$C_subtype))) else NA_integer_
  if (n_groups < 2) {
    return(tibble(
      n = n,
      n_groups = n_groups,
      min_group_n = min_group_n,
      kruskal_statistic = NA_real_,
      kruskal_df = NA_real_,
      kruskal_p = NA_real_,
      epsilon2 = NA_real_
    ))
  }
  kt <- kruskal.test(reformulate("C_subtype", response = score_col), data = df)
  h <- as.numeric(kt$statistic)
  eps2 <- (h - n_groups + 1) / (n - n_groups)
  tibble(
    n = n,
    n_groups = n_groups,
    min_group_n = min_group_n,
    kruskal_statistic = h,
    kruskal_df = as.numeric(kt$parameter),
    kruskal_p = kt$p.value,
    epsilon2 = max(0, min(1, eps2))
  )
}

pairwise_one <- function(df, score_col) {
  bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
    comp <- comparisons[i, ]
    g1 <- comp$group_1
    g2 <- comp$group_2
    x <- df %>% filter(C_subtype == g1) %>% pull(all_of(score_col))
    y <- df %>% filter(C_subtype == g2) %>% pull(all_of(score_col))
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    p <- if (length(x) > 0 && length(y) > 0) {
      tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    tibble(
      comparison = comp$comparison,
      comparison_class = comp$comparison_class,
      group_1 = g1,
      group_2 = g2,
      n_group_1 = length(x),
      n_group_2 = length(y),
      median_group_1 = median(x, na.rm = TRUE),
      median_group_2 = median(y, na.rm = TRUE),
      median_diff_group_1_minus_group_2 = median_group_1 - median_group_2,
      wilcox_p = p
    )
  }))
}

summarise_score_by_cancer <- function(df, score_col) {
  df %>%
    group_by(Cancer, C_subtype, C_label) %>%
    summarise(
      n = sum(!is.na(.data[[score_col]])),
      mean_value = mean(.data[[score_col]], na.rm = TRUE),
      median_value = median(.data[[score_col]], na.rm = TRUE),
      q25_value = safe_quantile(.data[[score_col]], 0.25),
      q75_value = safe_quantile(.data[[score_col]], 0.75),
      .groups = "drop"
    ) %>%
    mutate(score = score_col, .before = 1) %>%
    arrange(score, Cancer, C_subtype)
}

write_metric_set <- function(df, score_col, prefix) {
  summary_by_cancer <- summarise_score_by_cancer(df, score_col)
  write_tsv(summary_by_cancer, file.path(out_dir, paste0(prefix, "_summary_by_cancer.tsv")))

  equalweight_summary <- summary_by_cancer %>%
    group_by(score, C_subtype, C_label) %>%
    summarise(
      n_cancers = n_distinct(Cancer[n > 0]),
      total_n = sum(n),
      mean_of_cancer_medians = mean(median_value, na.rm = TRUE),
      median_of_cancer_medians = median(median_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(score, C_subtype)
  write_tsv(equalweight_summary, file.path(out_dir, paste0(prefix, "_summary_equalweight.tsv")))

  kruskal_by_cancer <- df %>%
    group_by(Cancer) %>%
    group_modify(~ kruskal_one(.x, score_col)) %>%
    ungroup() %>%
    mutate(
      score = score_col,
      kruskal_p_adj_BH = p.adjust(kruskal_p, method = "BH"),
      .before = 1
    ) %>%
    arrange(kruskal_p)
  write_tsv(kruskal_by_cancer, file.path(out_dir, paste0(prefix, "_kruskal_by_cancer.tsv")))

  pairwise_by_cancer <- df %>%
    group_by(Cancer) %>%
    group_modify(~ pairwise_one(.x, score_col)) %>%
    ungroup() %>%
    group_by(Cancer) %>%
    mutate(wilcox_p_adj_within_cancer_BH = p.adjust(wilcox_p, method = "BH")) %>%
    ungroup() %>%
    mutate(score = score_col, .before = 1) %>%
    arrange(Cancer, comparison)
  write_tsv(pairwise_by_cancer, file.path(out_dir, paste0(prefix, "_pairwise_by_cancer.tsv")))

  pairwise_summary <- pairwise_by_cancer %>%
    group_by(score, comparison, comparison_class, group_1, group_2) %>%
    summarise(
      n_cancers_tested = sum(!is.na(wilcox_p)),
      n_p_lt_0.05 = sum(wilcox_p < 0.05, na.rm = TRUE),
      n_p_adj_lt_0.05 = sum(wilcox_p_adj_within_cancer_BH < 0.05, na.rm = TRUE),
      n_group_1_gt_group_2 = sum(median_diff_group_1_minus_group_2 > 0, na.rm = TRUE),
      n_group_2_gt_group_1 = sum(median_diff_group_1_minus_group_2 < 0, na.rm = TRUE),
      median_of_cancer_median_diff = median(median_diff_group_1_minus_group_2, na.rm = TRUE),
      paired_wilcox_p = {
        diff <- median_diff_group_1_minus_group_2[!is.na(median_diff_group_1_minus_group_2)]
        if (length(diff) > 0) {
          tryCatch(wilcox.test(diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(paired_wilcox_p_adj_BH = p.adjust(paired_wilcox_p, method = "BH")) %>%
    arrange(score, comparison)
  write_tsv(pairwise_summary, file.path(out_dir, paste0(prefix, "_pairwise_summary_by_contrast.tsv")))

  wide <- summary_by_cancer %>%
    select(Cancer, C_subtype, median_value) %>%
    pivot_wider(names_from = C_subtype, values_from = median_value) %>%
    filter(if_all(all_of(c_levels), ~ !is.na(.x))) %>%
    arrange(Cancer)

  paired_pairwise <- bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
    comp <- comparisons[i, ]
    g1 <- comp$group_1
    g2 <- comp$group_2
    diff <- wide[[g1]] - wide[[g2]]
    p <- if (length(diff[!is.na(diff)]) > 0) {
      tryCatch(wilcox.test(diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    tibble(
      score = score_col,
      comparison = comp$comparison,
      comparison_class = comp$comparison_class,
      group_1 = g1,
      group_2 = g2,
      n_cancers = length(diff),
      n_group_1_gt_group_2 = sum(diff > 0, na.rm = TRUE),
      n_group_2_gt_group_1 = sum(diff < 0, na.rm = TRUE),
      median_diff_group_1_minus_group_2 = median(diff, na.rm = TRUE),
      paired_wilcox_p = p
    )
  })) %>%
    mutate(paired_wilcox_p_adj_BH = p.adjust(paired_wilcox_p, method = "BH")) %>%
    arrange(score, comparison)
  write_tsv(paired_pairwise, file.path(out_dir, paste0(prefix, "_pan_cancer_paired_pairwise.tsv")))

  invisible(list(
    summary_by_cancer = summary_by_cancer,
    equalweight_summary = equalweight_summary,
    kruskal_by_cancer = kruskal_by_cancer,
    pairwise_by_cancer = pairwise_by_cancer,
    pairwise_summary = pairwise_summary,
    paired_pairwise = paired_pairwise
  ))
}

message("=== 6.0 C1-C4 axis validation ===")
message("Input labels: ", label_file)
message("Output dir: ", out_dir)

feature_metadata <- read_tsv(feature_metadata_file, show_col_types = FALSE) %>%
  arrange(display_order)
feature_ids <- feature_metadata$feature_id
inflammatory_features <- setdiff(feature_ids, "LM12_Antimicrobial_Peptides")

labels <- read_csv(label_file, show_col_types = FALSE)
required_cols <- c("Sample", "Cancer", "Acid_Score", "AG", "IG", "AGIG", "C_subtype", "C_label", feature_ids)
missing_cols <- setdiff(required_cols, names(labels))
if (length(missing_cols) > 0) {
  stop("标签表缺失必要列: ", paste(missing_cols, collapse = ", "))
}

labels <- labels %>%
  mutate(
    C_subtype = factor(as.character(C_subtype), levels = c_levels),
    AG = factor(as.character(AG), levels = c("AG1", "AG2", "AG3")),
    IG = factor(as.character(IG), levels = c("IG1", "IG2"))
  ) %>%
  group_by(Cancer) %>%
  mutate(
    Acid_Score_z = safe_z(Acid_Score),
    LM12_inflammatory_axis = rowMeans(across(all_of(inflammatory_features)), na.rm = TRUE),
    LM12_all_feature_axis = rowMeans(across(all_of(feature_ids)), na.rm = TRUE),
    LM12_inflammatory_axis_z = safe_z(LM12_inflammatory_axis),
    LM12_all_feature_axis_z = safe_z(LM12_all_feature_axis)
  ) %>%
  ungroup()

feature_z_df <- labels %>%
  select(Sample, Cancer, C_subtype, all_of(feature_ids)) %>%
  group_by(Cancer) %>%
  mutate(across(all_of(feature_ids), safe_z, .names = "{.col}_z")) %>%
  ungroup() %>%
  select(Sample, Cancer, C_subtype, ends_with("_z"))

axis_sample_scores <- labels %>%
  select(
    Sample, Cancer, AG, IG, AGIG, C_subtype, C_label,
    Acid_Score, Acid_Score_z,
    LM12_inflammatory_axis, LM12_inflammatory_axis_z,
    LM12_all_feature_axis, LM12_all_feature_axis_z,
    all_of(feature_ids)
  ) %>%
  left_join(feature_z_df, by = c("Sample", "Cancer", "C_subtype"))
write_tsv(axis_sample_scores, file.path(out_dir, "c1c4_axis_sample_scores.tsv.gz"))

message("[1/4] 样本量与状态中心")
counts_by_cancer <- labels %>%
  count(Cancer, C_subtype, name = "n") %>%
  complete(Cancer, C_subtype = factor(c_levels, levels = c_levels), fill = list(n = 0)) %>%
  group_by(Cancer) %>%
  mutate(
    total_samples = sum(n),
    all_groups_present = all(n > 0),
    all_groups_n_ge_10 = all(n >= 10),
    all_groups_n_ge_20 = all(n >= 20),
    min_group_n = min(n)
  ) %>%
  ungroup() %>%
  arrange(Cancer, C_subtype)
write_tsv(counts_by_cancer, file.path(out_dir, "c1c4_counts_by_cancer_long.tsv"))

centroids_by_cancer <- labels %>%
  group_by(Cancer, C_subtype, C_label) %>%
  summarise(
    n = n(),
    median_acid_z = median(Acid_Score_z, na.rm = TRUE),
    median_inflammatory_axis_z = median(LM12_inflammatory_axis_z, na.rm = TRUE),
    median_all_lm12_axis_z = median(LM12_all_feature_axis_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Cancer, C_subtype)
write_tsv(centroids_by_cancer, file.path(out_dir, "c1c4_state_centroids_by_cancer.tsv"))

centroids_equalweight <- centroids_by_cancer %>%
  group_by(C_subtype, C_label) %>%
  summarise(
    n_cancers = n_distinct(Cancer[n > 0]),
    total_n = sum(n),
    mean_median_acid_z = mean(median_acid_z, na.rm = TRUE),
    median_median_acid_z = median(median_acid_z, na.rm = TRUE),
    mean_median_inflammatory_axis_z = mean(median_inflammatory_axis_z, na.rm = TRUE),
    median_median_inflammatory_axis_z = median(median_inflammatory_axis_z, na.rm = TRUE),
    mean_median_all_lm12_axis_z = mean(median_all_lm12_axis_z, na.rm = TRUE),
    median_median_all_lm12_axis_z = median(median_all_lm12_axis_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(C_subtype)
write_tsv(centroids_equalweight, file.path(out_dir, "c1c4_state_centroids_equalweight.tsv"))

message("[2/4] Acid Score 统计")
acid_results <- write_metric_set(labels, "Acid_Score_z", "acid_score_z")

message("[3/4] LM12 inflammatory-remodeling axis 统计")
inflammatory_results <- write_metric_set(labels, "LM12_inflammatory_axis_z", "lm12_inflammatory_axis_z")

message("[4/4] 12 个 LM-Imm feature 统计")
feature_long <- labels %>%
  select(Sample, Cancer, C_subtype, C_label, all_of(feature_ids)) %>%
  pivot_longer(cols = all_of(feature_ids), names_to = "feature_id", values_to = "value") %>%
  group_by(Cancer, feature_id) %>%
  mutate(value_z = safe_z(value)) %>%
  ungroup() %>%
  left_join(feature_metadata, by = "feature_id")

lm12_feature_summary_by_cancer <- feature_long %>%
  group_by(feature_id, short_label, display_order, Cancer, C_subtype, C_label) %>%
  summarise(
    n = sum(!is.na(value_z)),
    mean_value_z = mean(value_z, na.rm = TRUE),
    median_value_z = median(value_z, na.rm = TRUE),
    q25_value_z = safe_quantile(value_z, 0.25),
    q75_value_z = safe_quantile(value_z, 0.75),
    .groups = "drop"
  ) %>%
  arrange(display_order, Cancer, C_subtype)
write_tsv(lm12_feature_summary_by_cancer, file.path(out_dir, "lm12_feature_summary_by_cancer.tsv"))

lm12_feature_summary_equalweight <- lm12_feature_summary_by_cancer %>%
  group_by(feature_id, short_label, display_order, C_subtype, C_label) %>%
  summarise(
    n_cancers = n_distinct(Cancer[n > 0]),
    total_n = sum(n),
    mean_of_cancer_medians_z = mean(median_value_z, na.rm = TRUE),
    median_of_cancer_medians_z = median(median_value_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(display_order, C_subtype)
write_tsv(lm12_feature_summary_equalweight, file.path(out_dir, "lm12_feature_summary_equalweight.tsv"))

lm12_feature_pairwise_by_cancer <- feature_long %>%
  group_by(feature_id, short_label, display_order, Cancer) %>%
  group_modify(~ pairwise_one(.x, "value_z")) %>%
  ungroup() %>%
  group_by(feature_id, Cancer) %>%
  mutate(wilcox_p_adj_within_feature_cancer_BH = p.adjust(wilcox_p, method = "BH")) %>%
  ungroup() %>%
  arrange(display_order, Cancer, comparison)
write_tsv(lm12_feature_pairwise_by_cancer, file.path(out_dir, "lm12_feature_pairwise_by_cancer.tsv"))

lm12_feature_pairwise_summary <- lm12_feature_pairwise_by_cancer %>%
  group_by(feature_id, short_label, display_order, comparison, comparison_class, group_1, group_2) %>%
  summarise(
    n_cancers_tested = sum(!is.na(wilcox_p)),
    n_p_lt_0.05 = sum(wilcox_p < 0.05, na.rm = TRUE),
    n_p_adj_lt_0.05 = sum(wilcox_p_adj_within_feature_cancer_BH < 0.05, na.rm = TRUE),
    n_group_1_gt_group_2 = sum(median_diff_group_1_minus_group_2 > 0, na.rm = TRUE),
    n_group_2_gt_group_1 = sum(median_diff_group_1_minus_group_2 < 0, na.rm = TRUE),
    median_of_cancer_median_diff = median(median_diff_group_1_minus_group_2, na.rm = TRUE),
    paired_wilcox_p = {
      diff <- median_diff_group_1_minus_group_2[!is.na(median_diff_group_1_minus_group_2)]
      if (length(diff) > 0) {
        tryCatch(wilcox.test(diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  group_by(feature_id) %>%
  mutate(paired_wilcox_p_adj_within_feature_BH = p.adjust(paired_wilcox_p, method = "BH")) %>%
  ungroup() %>%
  arrange(display_order, comparison)
write_tsv(lm12_feature_pairwise_summary, file.path(out_dir, "lm12_feature_pairwise_summary_by_contrast.tsv"))

axis_focus_summary <- bind_rows(
  acid_results$pairwise_summary %>%
    mutate(axis = "Acid_Score_z", feature_id = NA_character_, short_label = "Acid Score"),
  inflammatory_results$pairwise_summary %>%
    mutate(axis = "LM12_inflammatory_axis_z", feature_id = NA_character_, short_label = "LM12 inflammatory axis")
) %>%
  select(axis, feature_id, short_label, everything()) %>%
  arrange(axis, comparison)
write_tsv(axis_focus_summary, file.path(out_dir, "c1c4_axis_pairwise_summary_by_contrast.tsv"))

run_summary_lines <- c(
  "analysis\tC1C4_axis_validation",
  paste0("created_at\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("n_cancers\t", n_distinct(labels$Cancer)),
  paste0("n_samples\t", nrow(labels)),
  paste0("n_features\t", length(feature_ids)),
  paste0("feature_ids\t", paste(feature_ids, collapse = ";")),
  paste0("comparisons\t", paste(comparisons$comparison, collapse = ";")),
  paste0("acid_C1_vs_C4_n_group_1_gt_group_2\t", acid_results$pairwise_summary$n_group_1_gt_group_2[acid_results$pairwise_summary$comparison == "C1_vs_C4"]),
  paste0("acid_C1_vs_C4_n_group_2_gt_group_1\t", acid_results$pairwise_summary$n_group_2_gt_group_1[acid_results$pairwise_summary$comparison == "C1_vs_C4"]),
  paste0("lm12_axis_C1_vs_C2_n_group_1_gt_group_2\t", inflammatory_results$pairwise_summary$n_group_1_gt_group_2[inflammatory_results$pairwise_summary$comparison == "C1_vs_C2"]),
  paste0("lm12_axis_C3_vs_C4_n_group_1_gt_group_2\t", inflammatory_results$pairwise_summary$n_group_1_gt_group_2[inflammatory_results$pairwise_summary$comparison == "C3_vs_C4"])
)
writeLines(run_summary_lines, file.path(out_dir, "run_summary.txt"))

message("=== 02_validate_c1c4_subtype_axes.R complete ===")
message("Saved outputs to: ", out_dir)
}

# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/shared_definitions/03_export_canonical_labels.R
stage_shared_definitions_03_export_canonical_labels <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Construct IG and C1-C4 labels and export the canonical patient assignments.
# Task: c1c4_state_assignment; step: s03
# Inputs (DCF fields): assignment_labels,axis_scores,feature_metadata
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_state_assignment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
label_01_file <- unname(paths["assignment_labels"])
axis_score_file <- unname(paths["axis_scores"])
feature_metadata_file <- unname(paths["feature_metadata"])
out_dir <- unname(paths["output_dir"])
required_files <- c(label_01_file, axis_score_file, feature_metadata_file)

# 防止输出覆盖输入、配置或既有结果；检查模式不会创建目录。
protected_paths <- c(dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

canonical_version <- "AG3IG2_C1C4_v1_20260615"
created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

required_files <- c(label_01_file, axis_score_file, feature_metadata_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("缺失输入文件:\n", paste(missing_files, collapse = "\n"))
}

c_levels <- c("C1", "C2", "C3", "C4")
ag_levels <- c("AG1", "AG2", "AG3")
ig_levels <- c("IG1", "IG2")
c_map <- c(
  AG1_IG1 = "C1",
  AG1_IG2 = "C2",
  AG2_IG2 = "C2",
  AG2_IG1 = "C3",
  AG3_IG1 = "C3",
  AG3_IG2 = "C4"
)

message("=== 6.0 freeze canonical AG3 x IG2 -> C1-C4 labels ===")
message("Input 01 labels: ", label_01_file)
message("Input 02 axis scores: ", axis_score_file)
message("Output dir: ", out_dir)

feature_metadata <- read_tsv(feature_metadata_file, show_col_types = FALSE) %>%
  arrange(display_order)
feature_ids <- feature_metadata$feature_id
feature_z_cols <- paste0(feature_ids, "_z")

labels_01 <- read_csv(label_01_file, show_col_types = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    Cancer = as.character(Cancer),
    AG = as.character(AG),
    IG = as.character(IG),
    AGIG = as.character(AGIG),
    C_subtype = as.character(C_subtype)
  )

axis_scores <- read_tsv(axis_score_file, show_col_types = FALSE) %>%
  mutate(
    Sample = as.character(Sample),
    Cancer = as.character(Cancer),
    AG = as.character(AG),
    IG = as.character(IG),
    AGIG = as.character(AGIG),
    C_subtype = as.character(C_subtype)
  )

required_01_cols <- c(
  "Sample", "Cancer", "Acid_Group", "AG", "AG_Label", "IG2_raw_cluster",
  "IG", "IG_Label", "immune_rank_within_cancer", "AGIG", "C_subtype",
  "C_label", feature_ids
)
required_axis_cols <- c(
  "Sample", "Cancer", "AG", "IG", "AGIG", "C_subtype", "C_label",
  "Acid_Score", "Acid_Score_z", "LM12_inflammatory_axis",
  "LM12_inflammatory_axis_z", "LM12_all_feature_axis",
  "LM12_all_feature_axis_z", feature_ids, feature_z_cols
)

missing_01_cols <- setdiff(required_01_cols, names(labels_01))
missing_axis_cols <- setdiff(required_axis_cols, names(axis_scores))
if (length(missing_01_cols) > 0) {
  stop("01 标签表缺失必要列: ", paste(missing_01_cols, collapse = ", "))
}
if (length(missing_axis_cols) > 0) {
  stop("02 轴得分表缺失必要列: ", paste(missing_axis_cols, collapse = ", "))
}

message("[1/5] 输入一致性校验")
dup_01 <- labels_01 %>% count(Sample, Cancer) %>% filter(n > 1)
dup_axis <- axis_scores %>% count(Sample, Cancer) %>% filter(n > 1)
if (nrow(dup_01) > 0 || nrow(dup_axis) > 0) {
  stop("存在重复 Sample-Cancer 键，无法冻结 canonical label。")
}

compare_cols <- c("AG", "IG", "AGIG", "C_subtype")
joined_check <- labels_01 %>%
  select(Sample, Cancer, all_of(compare_cols)) %>%
  inner_join(
    axis_scores %>% select(Sample, Cancer, all_of(compare_cols)),
    by = c("Sample", "Cancer"),
    suffix = c("_01", "_axis")
  )

if (nrow(joined_check) != nrow(labels_01) || nrow(joined_check) != nrow(axis_scores)) {
  stop("01 标签表和 02 轴得分表的 Sample-Cancer 键不完全一致。")
}

for (col in compare_cols) {
  mismatch_n <- sum(joined_check[[paste0(col, "_01")]] != joined_check[[paste0(col, "_axis")]], na.rm = TRUE)
  if (mismatch_n > 0) {
    stop("01 与 02 的 ", col, " 不一致，mismatch_n = ", mismatch_n)
  }
}

mapped_c <- unname(c_map[labels_01$AGIG])
if (any(is.na(mapped_c)) || any(mapped_c != labels_01$C_subtype)) {
  bad <- labels_01 %>% filter(is.na(mapped_c) | mapped_c != C_subtype) %>% distinct(AGIG, C_subtype)
  stop("AGIG -> C_subtype 确定性映射校验失败:\n", paste(capture.output(print(bad)), collapse = "\n"))
}

message("[2/5] 构建 canonical full/minimal label 表")
metadata_cols <- c(
  "Sample", "Cancer", "Acid_Group", "AG_Label", "IG2_raw_cluster",
  "IG_Label", "immune_rank_within_cancer"
)
labels_meta <- labels_01 %>%
  select(all_of(metadata_cols))

canonical_labels <- axis_scores %>%
  left_join(labels_meta, by = c("Sample", "Cancer")) %>%
  mutate(
    canonical_version = canonical_version,
    C_subtype = factor(C_subtype, levels = c_levels),
    AG = factor(AG, levels = ag_levels),
    IG = factor(IG, levels = ig_levels),
    AGIG = factor(AGIG, levels = names(c_map)),
    freeze_date = created_at
  ) %>%
  arrange(Cancer, C_subtype, AG, IG, Sample) %>%
  mutate(canonical_row_id = row_number()) %>%
  select(
    canonical_version,
    canonical_row_id,
    freeze_date,
    Sample,
    Cancer,
    C_subtype,
    C_label,
    AG,
    AG_Label,
    Acid_Group,
    IG,
    IG_Label,
    AGIG,
    IG2_raw_cluster,
    immune_rank_within_cancer,
    Acid_Score,
    Acid_Score_z,
    LM12_inflammatory_axis,
    LM12_inflammatory_axis_z,
    LM12_all_feature_axis,
    LM12_all_feature_axis_z,
    all_of(feature_ids),
    all_of(feature_z_cols)
  )

canonical_minimal <- canonical_labels %>%
  select(
    canonical_version,
    canonical_row_id,
    freeze_date,
    Sample,
    Cancer,
    C_subtype,
    C_label,
    AG,
    AG_Label,
    Acid_Group,
    IG,
    IG_Label,
    AGIG,
    Acid_Score,
    Acid_Score_z,
    LM12_inflammatory_axis,
    LM12_inflammatory_axis_z
  )

message("[3/5] 生成样本量和字典")
counts_long <- canonical_labels %>%
  count(Cancer, C_subtype, C_label, name = "n_samples") %>%
  complete(
    Cancer,
    C_subtype = factor(c_levels, levels = c_levels),
    fill = list(n_samples = 0)
  ) %>%
  group_by(Cancer) %>%
  mutate(
    total_samples = sum(n_samples),
    all_groups_present = all(n_samples > 0),
    all_groups_n_ge_10 = all(n_samples >= 10),
    all_groups_n_ge_20 = all(n_samples >= 20),
    min_group_n = min(n_samples)
  ) %>%
  ungroup() %>%
  arrange(Cancer, C_subtype)

label_lookup <- canonical_labels %>%
  distinct(C_subtype, C_label)
counts_long <- counts_long %>%
  select(-C_label) %>%
  left_join(label_lookup, by = "C_subtype") %>%
  select(Cancer, C_subtype, C_label, everything())

counts_wide <- counts_long %>%
  select(Cancer, C_subtype, n_samples) %>%
  pivot_wider(names_from = C_subtype, values_from = n_samples) %>%
  mutate(
    total_samples = rowSums(across(all_of(c_levels))),
    min_group_n = do.call(pmin, c(across(all_of(c_levels)), list(na.rm = TRUE))),
    all_groups_present = if_all(all_of(c_levels), ~ .x > 0),
    all_groups_n_ge_10 = if_all(all_of(c_levels), ~ .x >= 10),
    all_groups_n_ge_20 = if_all(all_of(c_levels), ~ .x >= 20)
  ) %>%
  arrange(Cancer)

counts_overall <- canonical_labels %>%
  count(C_subtype, C_label, name = "n_samples") %>%
  complete(C_subtype = factor(c_levels, levels = c_levels), fill = list(n_samples = 0)) %>%
  select(C_subtype, n_samples) %>%
  left_join(label_lookup, by = "C_subtype") %>%
  select(C_subtype, C_label, n_samples)

subtype_dictionary <- tribble(
  ~field, ~level, ~label, ~definition, ~construction_rule, ~interpretation_note,
  "AG", "AG1", "acid_high", "Within-cancer high-acid group", "clinical_data$Group == High", "The acid axis is fixed as a within-cancer three-level stratification.",
  "AG", "AG2", "acid_intermediate", "Within-cancer intermediate-acid group", "clinical_data$Group == Medium", "Intermediate acid state is retained instead of forcing a binary contrast.",
  "AG", "AG3", "acid_low", "Within-cancer low-acid group", "clinical_data$Group == Low", "Low-acid endpoint.",
  "IG", "IG1", "lm12_inflammatory_remodeling_high", "Higher LM12 inflammatory-remodeling program", "Within-cancer hierarchical clustering k=2 on 12 LM-Imm scores; clusters are oriented by the LM12 inflammatory axis mean.", "IG is not interpreted as immune infiltration abundance or anti-tumor immune strength.",
  "IG", "IG2", "lm12_inflammatory_remodeling_low", "Lower LM12 inflammatory-remodeling program", "Within-cancer hierarchical clustering k=2 on 12 LM-Imm scores; clusters are oriented by the LM12 inflammatory axis mean.", "IG is not interpreted as immune infiltration abundance or anti-tumor immune strength.",
  "C_subtype", "C1", "high_acid__high_inflammatory_remodeling", "High acid plus high inflammatory-remodeling program", "C1 = AG1_IG1", "Primary endpoint contrast with C4.",
  "C_subtype", "C2", "high_or_intermediate_acid__low_inflammatory_remodeling", "High/intermediate acid plus low inflammatory-remodeling program", "C2 = AG1_IG2 + AG2_IG2", "Useful for disentangling acid and inflammatory-remodeling axes against C1 or C3.",
  "C_subtype", "C3", "intermediate_or_low_acid__high_inflammatory_remodeling", "Intermediate/low acid plus high inflammatory-remodeling program", "C3 = AG2_IG1 + AG3_IG1", "Useful for disentangling acid and inflammatory-remodeling axes against C2 or C4.",
  "C_subtype", "C4", "low_acid__low_inflammatory_remodeling", "Low acid plus low inflammatory-remodeling program", "C4 = AG3_IG2", "Primary endpoint contrast with C1."
)

message("[4/5] 写出冻结文件")
full_csv <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_labels.csv")
minimal_csv <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_labels_minimal.csv")
rdata_file <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_labels.Rdata")
dictionary_file <- file.path(out_dir, "final_ag3ig2_c1c4_subtype_dictionary.tsv")
counts_long_file <- file.path(out_dir, "final_ag3ig2_c1c4_counts_by_cancer_long.tsv")
counts_wide_file <- file.path(out_dir, "final_ag3ig2_c1c4_counts_by_cancer_wide.tsv")
counts_overall_file <- file.path(out_dir, "final_ag3ig2_c1c4_counts_overall.tsv")
manifest_file <- file.path(out_dir, "final_ag3ig2_c1c4_freeze_manifest.tsv")
checksums_file <- file.path(out_dir, "final_ag3ig2_c1c4_file_checksums.tsv")

write_csv(canonical_labels, full_csv)
write_csv(canonical_minimal, minimal_csv)
write_tsv(subtype_dictionary, dictionary_file)
write_tsv(counts_long, counts_long_file)
write_tsv(counts_wide, counts_wide_file)
write_tsv(counts_overall, counts_overall_file)

freeze_metadata <- list(
  canonical_version = canonical_version,
  created_at = created_at,
  n_samples = nrow(canonical_labels),
  n_cancers = n_distinct(canonical_labels$Cancer),
  c_levels = c_levels,
  c_map = c_map,
  input_label_01_file = label_01_file,
  input_axis_score_file = axis_score_file,
  feature_metadata_file = feature_metadata_file,
  no_write_to_data_dir = TRUE
)

subtype_labels <- canonical_labels
subtype_labels_minimal <- canonical_minimal
save(
  subtype_labels,
  subtype_labels_minimal,
  subtype_dictionary,
  counts_long,
  counts_wide,
  counts_overall,
  freeze_metadata,
  file = rdata_file
)

manifest <- tibble(
  key = c(
    "canonical_version",
    "created_at",
    "n_samples",
    "n_cancers",
    "n_features",
    "n_cancers_all_groups_present",
    "n_cancers_all_groups_n_ge_10",
    "n_cancers_all_groups_n_ge_20",
    "input_label_01_file",
    "input_label_01_md5",
    "input_axis_score_file",
    "input_axis_score_md5",
    "feature_metadata_file",
    "feature_metadata_md5",
    "write_to_data_dir",
    "legacy_data_label_note"
  ),
  value = c(
    canonical_version,
    created_at,
    as.character(nrow(canonical_labels)),
    as.character(n_distinct(canonical_labels$Cancer)),
    as.character(length(feature_ids)),
    as.character(sum(counts_wide$all_groups_present)),
    as.character(sum(counts_wide$all_groups_n_ge_10)),
    as.character(sum(counts_wide$all_groups_n_ge_20)),
    label_01_file,
    unname(tools::md5sum(label_01_file)),
    axis_score_file,
    unname(tools::md5sum(axis_score_file)),
    feature_metadata_file,
    unname(tools::md5sum(feature_metadata_file)),
    "FALSE",
    "data/final_pancancer_subtype_labels.csv is legacy LM-Imm C1/C2/C3 and is not used by the current C1-C4 mainline."
  )
)
write_tsv(manifest, manifest_file)

message("[5/5] 输出校验和运行摘要")
output_files <- c(
  full_csv,
  minimal_csv,
  rdata_file,
  dictionary_file,
  counts_long_file,
  counts_wide_file,
  counts_overall_file,
  manifest_file
)
checksums <- tibble(
  file = output_files,
  md5 = unname(tools::md5sum(output_files)),
  size_bytes = file.info(output_files)$size
)
write_tsv(checksums, checksums_file)

summary_lines <- c(
  "analysis\tFreeze_AG3IG2_C1C4_canonical_labels",
  paste0("created_at\t", created_at),
  paste0("canonical_version\t", canonical_version),
  paste0("n_samples\t", nrow(canonical_labels)),
  paste0("n_cancers\t", n_distinct(canonical_labels$Cancer)),
  paste0("n_features\t", length(feature_ids)),
  paste0("counts_overall\t", paste(paste0(counts_overall$C_subtype, ":", counts_overall$n_samples), collapse = ";")),
  paste0("n_cancers_all_groups_present\t", sum(counts_wide$all_groups_present)),
  paste0("n_cancers_all_groups_n_ge_10\t", sum(counts_wide$all_groups_n_ge_10)),
  paste0("n_cancers_all_groups_n_ge_20\t", sum(counts_wide$all_groups_n_ge_20)),
  paste0("full_csv\t", full_csv),
  paste0("minimal_csv\t", minimal_csv),
  paste0("rdata_file\t", rdata_file),
  "write_to_data_dir\tFALSE"
)
writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))

message("=== 04_freeze_c1c4_labels.R complete ===")
message("Saved outputs to: ", out_dir)
}

# Purpose: Score external bulk cohorts and test the fixed LM12 pathway enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/external_bulk/01_build_external_scores.R
stage_external_bulk_01_build_external_scores <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Score external bulk cohorts and test the fixed LM12 pathway enrichment.
# Task: external_bulk_enrichment; step: s01
# Inputs (DCF fields): expression_dir,metadata_dir,final_genes,lm12_manifest,lm12_genes
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task external_bulk_enrichment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
geo_source_dir <- unname(paths["expression_dir"])
metadata_dir <- unname(paths["metadata_dir"])
acid_gene_file <- unname(paths["final_genes"])
lm_imm_12_manifest_file <- unname(paths["lm12_manifest"])
lm_imm_12_gene_long_file <- unname(paths["lm12_genes"])
out_dir <- unname(paths["output_dir"])
manifest_file <- file.path(metadata_dir, "external_cohort_manifest.tsv")
gse40419_annot_file <- file.path(metadata_dir, "GSE40419_sample_annotation.tsv")
gse127165_annot_file <- file.path(metadata_dir, "GSE127165_sample_annotation.tsv")
gse96058_annot_file <- file.path(metadata_dir, "GSE96058_expression_metadata_manifest.tsv")
expr_file_map <- c(
  GSE40419 = file.path(geo_source_dir, "GSE40419_LC-87_RPKM_expression.txt.gz"),
  GSE127165 = file.path(geo_source_dir, "GSE127165_114_samples_RNA_seq_FPKM.txt.gz"),
  GSE96058 = file.path(geo_source_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
)
gse40419_file <- unname(expr_file_map["GSE40419"])
gse127165_file <- unname(expr_file_map["GSE127165"])
gse96058_file <- unname(expr_file_map["GSE96058"])
required_files <- c(acid_gene_file, manifest_file, gse40419_annot_file,
  gse127165_annot_file, gse96058_annot_file, lm_imm_12_manifest_file,
  lm_imm_12_gene_long_file, expr_file_map)

# Inspect required fixed filenames without reading or decompressing their contents.
input_dirs <- c(geo_source_dir, metadata_dir)
if (any(!dir.exists(input_dirs)) || any(file.access(input_dirs, mode = 5L) != 0L)) {
  stop("Input directories must be readable and searchable.", call. = FALSE)
}
if (any(!file_test("-f", required_files)) || any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All fixed-name inputs must be readable nonempty regular files.", call. = FALSE)
}
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)


# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(input_dirs, dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named_object` <- bind_stage_helper(shared_load_named_object_026, environment())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GSVA)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readr)
})
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

gene_alias_map <- c(
  "PKM2" = "PKM"
)

normalize_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

apply_gene_alias_map <- function(gene_symbols) {
  alias_hit <- gene_symbols %in% names(gene_alias_map)
  gene_symbols[alias_hit] <- unname(gene_alias_map[gene_symbols[alias_hit]])
  list(
    gene_symbols = gene_symbols,
    alias_applied_n = sum(alias_hit, na.rm = TRUE)
  )
}

aggregate_duplicate_rows_mean <- function(expr_mat) {
  duplicate_n <- sum(duplicated(rownames(expr_mat)))
  if (duplicate_n == 0) {
    return(list(expr_mat = expr_mat, duplicate_n = 0L))
  }
  summed <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(
    matrix(1, nrow = nrow(expr_mat), ncol = ncol(expr_mat)),
    group = rownames(expr_mat),
    reorder = FALSE
  )
  list(expr_mat = summed / counts, duplicate_n = duplicate_n)
}

build_expression_qc <- function(cohort, gene_n_raw, gene_n_nonempty, gene_n_unique, duplicate_gene_n, sample_n, alias_applied_n) {
  data.frame(
    Cohort = cohort,
    Gene_N_Raw = gene_n_raw,
    Gene_N_NonEmpty = gene_n_nonempty,
    Gene_N_Unique = gene_n_unique,
    Duplicate_Gene_N = duplicate_gene_n,
    Alias_Applied_N = alias_applied_n,
    Sample_N = sample_n,
    stringsAsFactors = FALSE
  )
}

prepare_gse40419_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  gene_symbols <- normalize_gene_symbol(raw_df$gene)
  alias_res <- apply_gene_alias_map(gene_symbols)
  expr_df <- raw_df[, -(1:6), drop = FALSE]
  expr_mat <- as.matrix(data.frame(lapply(expr_df, as.numeric), check.names = FALSE))
  rownames(expr_mat) <- alias_res$gene_symbols
  keep <- nzchar(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  agg <- aggregate_duplicate_rows_mean(expr_mat)
  list(
    expr_mat = agg$expr_mat,
    qc_df = build_expression_qc(
      cohort = "GSE40419",
      gene_n_raw = nrow(raw_df),
      gene_n_nonempty = sum(keep),
      gene_n_unique = nrow(agg$expr_mat),
      duplicate_gene_n = agg$duplicate_n,
      sample_n = ncol(agg$expr_mat),
      alias_applied_n = alias_res$alias_applied_n
    )
  )
}

prepare_gse127165_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  ensembl_ids <- normalize_gene_symbol(sub("\\..*$", "", raw_df[[1]]))
  symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  gene_symbols <- normalize_gene_symbol(unname(symbol_map[ensembl_ids]))
  alias_res <- apply_gene_alias_map(gene_symbols)
  expr_df <- raw_df[, -1, drop = FALSE]
  expr_mat <- as.matrix(data.frame(lapply(expr_df, as.numeric), check.names = FALSE))
  rownames(expr_mat) <- alias_res$gene_symbols
  keep <- nzchar(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  agg <- aggregate_duplicate_rows_mean(expr_mat)
  list(
    expr_mat = agg$expr_mat,
    qc_df = build_expression_qc(
      cohort = "GSE127165",
      gene_n_raw = nrow(raw_df),
      gene_n_nonempty = sum(keep),
      gene_n_unique = nrow(agg$expr_mat),
      duplicate_gene_n = agg$duplicate_n,
      sample_n = ncol(agg$expr_mat),
      alias_applied_n = alias_res$alias_applied_n
    )
  )
}

prepare_gse96058_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = ",", data.table = FALSE)
  gene_symbols <- normalize_gene_symbol(raw_df[[1]])
  alias_res <- apply_gene_alias_map(gene_symbols)
  expr_df <- raw_df[, -1, drop = FALSE]
  expr_mat <- as.matrix(data.frame(lapply(expr_df, as.numeric), check.names = FALSE))
  rownames(expr_mat) <- alias_res$gene_symbols
  keep <- nzchar(rownames(expr_mat))
  expr_mat <- expr_mat[keep, , drop = FALSE]
  agg <- aggregate_duplicate_rows_mean(expr_mat)
  list(
    expr_mat = agg$expr_mat,
    qc_df = build_expression_qc(
      cohort = "GSE96058",
      gene_n_raw = nrow(raw_df),
      gene_n_nonempty = sum(keep),
      gene_n_unique = nrow(agg$expr_mat),
      duplicate_gene_n = agg$duplicate_n,
      sample_n = ncol(agg$expr_mat),
      alias_applied_n = alias_res$alias_applied_n
    )
  )
}

calc_coverage_tables <- function(cohort, expr_mat, acid_genes, lm_imm_12_manifest, lm_imm_12_gene_long) {
  expr_genes <- rownames(expr_mat)

  acid_detail_df <- data.frame(
    Cohort = cohort,
    Gene = acid_genes,
    Signature = "Acid_Score",
    Present = acid_genes %in% expr_genes,
    stringsAsFactors = FALSE
  )

  acid_summary_df <- data.frame(
    Cohort = cohort,
    Signature = "Acid_Score",
    Signature_Name = "Acid_Score",
    Gene_Set_Size = length(acid_genes),
    Covered_Gene_N = sum(acid_detail_df$Present),
    Coverage_Ratio = sum(acid_detail_df$Present) / length(acid_genes),
    stringsAsFactors = FALSE
  )

  lm_imm_12_summary_df <- lm_imm_12_gene_long %>%
    group_by(feature_id, pathway) %>%
    summarise(
      Gene_Set_Size = n_distinct(gene),
      Covered_Gene_N = sum(unique(gene) %in% expr_genes),
      .groups = "drop"
    ) %>%
    left_join(
      lm_imm_12_manifest %>% dplyr::select(feature_id, short_label, pathway),
      by = c("feature_id", "pathway")
    ) %>%
    mutate(
      Cohort = cohort,
      Signature = feature_id,
      Signature_Name = pathway,
      Coverage_Ratio = Covered_Gene_N / Gene_Set_Size
    ) %>%
    dplyr::select(Cohort, Signature, Signature_Name, Gene_Set_Size, Covered_Gene_N, Coverage_Ratio)

  list(
    acid_detail_df = acid_detail_df,
    signature_summary_df = bind_rows(acid_summary_df, lm_imm_12_summary_df)
  )
}

final_selected_genes <- load_named_object(acid_gene_file, "final_selected_genes")
if (!exists("final_selected_genes")) {
  stop("final_selected_genes.rdata 中未找到 final_selected_genes。")
}
acid_genes <- unique(as.character(final_selected_genes$Gene))

cohort_manifest <- fread(manifest_file, sep = "\t", data.table = FALSE)
gse40419_annot <- fread(gse40419_annot_file, sep = "\t", data.table = FALSE)
gse127165_annot <- fread(gse127165_annot_file, sep = "\t", data.table = FALSE)
gse96058_annot <- fread(gse96058_annot_file, sep = "\t", data.table = FALSE)
lm_imm_12_manifest <- read_tsv(lm_imm_12_manifest_file, show_col_types = FALSE)
lm_imm_12_gene_long <- read_tsv(lm_imm_12_gene_long_file, show_col_types = FALSE)

lm_imm_12_gene_sets <- split(
  lm_imm_12_gene_long$gene,
  lm_imm_12_gene_long$feature_id
)
score_col_order <- c("Acid_Score", lm_imm_12_manifest$feature_id)
signatures_for_scoring <- c(list(Acid_Score = acid_genes), lm_imm_12_gene_sets)

annotation_map <- list(
  GSE40419 = gse40419_annot,
  GSE127165 = gse127165_annot,
  GSE96058 = gse96058_annot
)

expr_list <- list(
  GSE40419 = prepare_gse40419_expr(gse40419_file),
  GSE127165 = prepare_gse127165_expr(gse127165_file),
  GSE96058 = prepare_gse96058_expr(gse96058_file)
)

all_score_rows <- list()
expression_qc_list <- list()
coverage_detail_list <- list()
coverage_summary_list <- list()

for (cohort in names(expr_list)) {
  expr_mat <- expr_list[[cohort]]$expr_mat
  expression_qc_list[[cohort]] <- expr_list[[cohort]]$qc_df

  coverage_res <- calc_coverage_tables(
    cohort = cohort,
    expr_mat = expr_mat,
    acid_genes = acid_genes,
    lm_imm_12_manifest = lm_imm_12_manifest,
    lm_imm_12_gene_long = lm_imm_12_gene_long
  )
  coverage_detail_list[[cohort]] <- coverage_res$acid_detail_df
  coverage_summary_list[[cohort]] <- coverage_res$signature_summary_df

  expr_log2 <- log2(pmax(expr_mat, 0) + 1)
  ssgsea_param <- ssgseaParam(
    exprData = expr_log2,
    geneSets = signatures_for_scoring,
    normalize = TRUE
  )
  score_mat <- gsva(ssgsea_param, verbose = FALSE)
  score_df <- as.data.frame(t(score_mat), check.names = FALSE, stringsAsFactors = FALSE)
  score_df$Sample_ID <- rownames(score_df)

  annot_df <- annotation_map[[cohort]]
  score_df <- annot_df %>%
    left_join(score_df, by = "Sample_ID") %>%
    dplyr::select(
      any_of(c(
        "Cohort", "Cancer", "Sample_ID", "Sample_Type", "Pair_ID",
        "Is_Replicate", "Primary_Sample_ID", "Has_Clinical_Metadata"
      )),
      all_of(score_col_order)
    )

  all_score_rows[[cohort]] <- score_df

  write.table(
    score_df,
    file = file.path(out_dir, paste0(cohort, "_score_matrix.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

external_score_matrix <- bind_rows(all_score_rows) %>%
  arrange(match(Cohort, c("GSE40419", "GSE127165", "GSE96058")), Sample_ID)

write.table(
  external_score_matrix,
  file = file.path(out_dir, "external_score_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bind_rows(expression_qc_list),
  file = file.path(out_dir, "external_expression_preparation_qc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bind_rows(coverage_detail_list),
  file = file.path(out_dir, "acid_gene_coverage_detail.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  bind_rows(coverage_summary_list),
  file = file.path(out_dir, "signature_coverage_qc.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("外部 GEO score matrix 已保存到: ", out_dir)
}

# Purpose: Score external bulk cohorts and test the fixed LM12 pathway enrichment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/external_bulk/02_run_lm12_gene_level_gsea.R
stage_external_bulk_02_run_lm12_gene_level_gsea <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Score external bulk cohorts and test the fixed LM12 pathway enrichment.
# Task: external_bulk_enrichment; step: s02
# Inputs (DCF fields): expression_dir,score_dir,lm12_manifest,lm12_genes
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task external_bulk_enrichment --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE)
geo_source_dir <- unname(paths["expression_dir"])
score_dir <- unname(paths["score_dir"])
manifest_file <- unname(paths["lm12_manifest"])
gene_long_file <- unname(paths["lm12_genes"])
out_dir <- unname(paths["output_dir"])
expr_file_map <- c(
  GSE40419 = file.path(geo_source_dir, "GSE40419_LC-87_RPKM_expression.txt.gz"),
  GSE127165 = file.path(geo_source_dir, "GSE127165_114_samples_RNA_seq_FPKM.txt.gz"),
  GSE96058 = file.path(geo_source_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
)
score_file_map <- setNames(
  file.path(score_dir, paste0(names(expr_file_map), "_score_matrix.tsv")),
  names(expr_file_map)
)
required_files <- c(manifest_file, gene_long_file, expr_file_map, score_file_map)

# Inspect required fixed filenames without reading or decompressing their contents.
input_dirs <- c(geo_source_dir, score_dir)
if (any(!dir.exists(input_dirs)) || any(file.access(input_dirs, mode = 5L) != 0L)) {
  stop("Input directories must be readable and searchable.", call. = FALSE)
}
if (any(!file_test("-f", required_files)) || any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All fixed-name inputs must be readable nonempty regular files.", call. = FALSE)
}
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)


# Protect inputs and existing results; metadata checks create no directories.
protected_paths <- c(input_dirs, dirname(required_files), config_path)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(limma)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
set.seed(20260804)
if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}

gene_alias_map <- c("PKM2" = "PKM")

normalize_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  hit <- x %in% names(gene_alias_map)
  x[hit] <- unname(gene_alias_map[x[hit]])
  x
}

aggregate_duplicate_rows_mean <- function(expr_mat) {
  if (!anyDuplicated(rownames(expr_mat))) return(expr_mat)
  summed <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(
    matrix(1, nrow = nrow(expr_mat), ncol = ncol(expr_mat)),
    group = rownames(expr_mat), reorder = FALSE
  )
  summed / counts
}

prepare_gse40419_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  gene_symbols <- normalize_gene_symbol(raw_df$gene)
  expr_mat <- as.matrix(data.frame(lapply(raw_df[, -(1:6), drop = FALSE], as.numeric), check.names = FALSE))
  rownames(expr_mat) <- gene_symbols
  expr_mat <- expr_mat[nzchar(rownames(expr_mat)), , drop = FALSE]
  aggregate_duplicate_rows_mean(expr_mat)
}

prepare_gse127165_expr <- function(path) {
  raw_df <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", data.table = FALSE)
  ensembl_ids <- sub("\\..*$", "", as.character(raw_df[[1]]))
  symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
  )
  gene_symbols <- normalize_gene_symbol(unname(symbol_map[ensembl_ids]))
  expr_mat <- as.matrix(data.frame(lapply(raw_df[, -1, drop = FALSE], as.numeric), check.names = FALSE))
  rownames(expr_mat) <- gene_symbols
  expr_mat <- expr_mat[nzchar(rownames(expr_mat)), , drop = FALSE]
  aggregate_duplicate_rows_mean(expr_mat)
}

prepare_gse96058_expr <- function(path) {
  # 该矩阵较大，保持data.table数值列并尽早释放原对象，避免data.frame/lapply产生多份全矩阵副本。
  raw_dt <- fread(cmd = paste("gzip -dc", shQuote(path)), sep = ",", data.table = TRUE)
  gene_symbols <- normalize_gene_symbol(raw_dt[[1]])
  raw_dt[[1]] <- NULL
  expr_mat <- as.matrix(raw_dt)
  rm(raw_dt)
  gc(verbose = FALSE)
  rownames(expr_mat) <- gene_symbols
  expr_mat <- expr_mat[nzchar(rownames(expr_mat)), , drop = FALSE]
  aggregate_duplicate_rows_mean(expr_mat)
}

expr_loader_map <- list(
  GSE40419 = prepare_gse40419_expr,
  GSE127165 = prepare_gse127165_expr,
  GSE96058 = prepare_gse96058_expr
)

manifest_df <- read_tsv(manifest_file, show_col_types = FALSE) %>% arrange(display_order)
gene_long_df <- read_tsv(gene_long_file, show_col_types = FALSE)
if (
  nrow(manifest_df) != 12 ||
    anyDuplicated(manifest_df$feature_id) > 0 ||
    !all(manifest_df$display_order == seq_len(12))
) {
  stop("LM12 manifest必须为12条唯一feature且display_order严格为1..12。")
}

lm12_gene_sets <- split(gene_long_df$gene, gene_long_df$feature_id)
lm12_gene_sets <- lapply(lm12_gene_sets[manifest_df$feature_id], function(x) {
  unique(normalize_gene_symbol(x[nzchar(x)]))
})

assign_high_low <- function(score_df) {
  score_df <- score_df %>% arrange(desc(Acid_Score), Sample_ID)
  high_n <- ceiling(nrow(score_df) / 2)
  score_df$Acid_Group <- c(rep("High_Acid", high_n), rep("Low_Acid", nrow(score_df) - high_n))
  score_df
}

run_one_analysis <- function(cohort, expr_log2, score_df, analysis_set) {
  common_samples <- intersect(score_df$Sample_ID, colnames(expr_log2))
  score_work <- score_df %>%
    filter(Sample_ID %in% common_samples, Sample_Type == "Tumor", is.finite(Acid_Score)) %>%
    distinct(Sample_ID, .keep_all = TRUE)
  score_work <- assign_high_low(score_work)

  expr_work <- expr_log2[, score_work$Sample_ID, drop = FALSE]
  group_factor <- factor(score_work$Acid_Group, levels = c("Low_Acid", "High_Acid"))
  design <- model.matrix(~ group_factor)
  fit <- eBayes(lmFit(expr_work, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$Gene <- rownames(tt)

  rank_df <- tt %>%
    transmute(Gene, logFC, AveExpr, moderated_t = t, P_Value = P.Value, FDR = adj.P.Val, B) %>%
    filter(nzchar(Gene), is.finite(moderated_t)) %>%
    arrange(desc(moderated_t), Gene)
  if (anyDuplicated(rank_df$Gene) > 0) stop(cohort, ": ranked genes存在重复。")
  ranks <- rank_df$moderated_t
  names(ranks) <- rank_df$Gene
  ranks <- sort(ranks, decreasing = TRUE)

  coverage_df <- manifest_df %>%
    transmute(
      display_order, feature_id, short_label, pathway, direction_label, functional_group,
      gene_set_size = vapply(lm12_gene_sets[feature_id], length, integer(1)),
      covered_gene_n = vapply(lm12_gene_sets[feature_id], function(gs) sum(gs %in% names(ranks)), integer(1)),
      coverage_ratio = covered_gene_n / gene_set_size
    )

  fgsea_res <- fgseaMultilevel(
    pathways = lm12_gene_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  ) %>%
    as.data.frame()
  fgsea_res$leadingEdge <- vapply(fgsea_res$leadingEdge, paste, collapse = ",", character(1))

  result_df <- manifest_df %>%
    dplyr::select(display_order, feature_id, short_label, pathway, direction_label, functional_group) %>%
    left_join(fgsea_res, by = c("feature_id" = "pathway")) %>%
    left_join(coverage_df %>% dplyr::select(feature_id, gene_set_size, covered_gene_n, coverage_ratio), by = "feature_id") %>%
    mutate(
      Cohort = cohort,
      Cancer = unique(score_work$Cancer),
      analysis_set = analysis_set,
      High_Acid_N = sum(score_work$Acid_Group == "High_Acid"),
      Low_Acid_N = sum(score_work$Acid_Group == "Low_Acid"),
      expected_nes_sign = ifelse(direction_label == "High_Acid_Positive", 1, -1),
      direction_consistent = sign(NES) == expected_nes_sign
    ) %>%
    dplyr::select(
      Cohort, Cancer, analysis_set, display_order, feature_id, short_label, pathway,
      direction_label, functional_group, High_Acid_N, Low_Acid_N,
      gene_set_size, covered_gene_n, coverage_ratio,
      size, ES, NES, pval, padj, log2err, leadingEdge,
      expected_nes_sign, direction_consistent
    ) %>%
    arrange(display_order)

  if (nrow(result_df) != 12 || any(!is.finite(result_df$NES))) {
    stop(cohort, " / ", analysis_set, ": 未产生12条有限NES结果。")
  }

  prefix <- paste0(cohort, "_lm12_gene_level_gsea", ifelse(analysis_set == "primary", "", paste0("_", analysis_set)))
  write_tsv(score_work %>% dplyr::select(Cohort, Cancer, Sample_ID, Acid_Score, Acid_Group, any_of(c("Is_Replicate", "Primary_Sample_ID"))), file.path(out_dir, paste0(prefix, "_sample_groups.tsv")))
  write_tsv(coverage_df, file.path(out_dir, paste0(prefix, "_coverage.tsv")))
  write_tsv(result_df, file.path(out_dir, paste0(prefix, "_results.tsv")))
  rank_con <- gzfile(file.path(out_dir, paste0(prefix, "_ranked_t.tsv.gz")), open = "wt")
  write.table(rank_df, rank_con, sep = "\t", row.names = FALSE, quote = FALSE)
  close(rank_con)
  saveRDS(
    list(
      cohort = cohort,
      cancer = unique(score_work$Cancer),
      analysis_set = analysis_set,
      sample_groups = score_work,
      pathway_order = manifest_df$feature_id,
      pathway_gene_sets = lm12_gene_sets,
      fgsea_results = result_df,
      ranked_stats = ranks
    ),
    file.path(out_dir, paste0(prefix, "_objects.rds"))
  )

  list(results = result_df, coverage = coverage_df, sample_groups = score_work, rank_n = length(ranks))
}

primary_results <- list()
qa_rows <- list()
sensitivity_results <- list()

for (cohort in names(expr_file_map)) {
  message("读取并处理: ", cohort)
  expr_mat <- expr_loader_map[[cohort]](expr_file_map[[cohort]])
  expr_log2 <- log2(pmax(expr_mat, 0) + 1)
  rm(expr_mat)
  gc(verbose = FALSE)

  score_df <- read_tsv(score_file_map[[cohort]], show_col_types = FALSE)
  score_primary <- score_df
  if (cohort == "GSE96058") {
    score_primary <- score_df %>% filter(!as.logical(Is_Replicate))
  }

  primary <- run_one_analysis(cohort, expr_log2, score_primary, "primary")
  primary_results[[cohort]] <- primary$results
  qa_rows[[cohort]] <- tibble(
    Cohort = cohort,
    Cancer = unique(primary$results$Cancer),
    analysis_set = "primary",
    Sample_N = nrow(primary$sample_groups),
    High_Acid_N = sum(primary$sample_groups$Acid_Group == "High_Acid"),
    Low_Acid_N = sum(primary$sample_groups$Acid_Group == "Low_Acid"),
    Ranked_Gene_N = primary$rank_n,
    Pathway_N = nrow(primary$results),
    Direction_Consistent_N = sum(primary$results$direction_consistent),
    FDR_lt_0_05_N = sum(primary$results$padj < 0.05)
  )

  if (cohort == "GSE96058") {
    sensitivity <- run_one_analysis(cohort, expr_log2, score_df, "all_samples_sensitivity")
    sensitivity_results[[cohort]] <- sensitivity$results
    qa_rows[[paste0(cohort, "_all")]] <- tibble(
      Cohort = cohort,
      Cancer = unique(sensitivity$results$Cancer),
      analysis_set = "all_samples_sensitivity",
      Sample_N = nrow(sensitivity$sample_groups),
      High_Acid_N = sum(sensitivity$sample_groups$Acid_Group == "High_Acid"),
      Low_Acid_N = sum(sensitivity$sample_groups$Acid_Group == "Low_Acid"),
      Ranked_Gene_N = sensitivity$rank_n,
      Pathway_N = nrow(sensitivity$results),
      Direction_Consistent_N = sum(sensitivity$results$direction_consistent),
      FDR_lt_0_05_N = sum(sensitivity$results$padj < 0.05)
    )
    rm(sensitivity)
    gc(verbose = FALSE)
  }

  rm(expr_log2, score_df, score_primary)
  gc(verbose = FALSE)
}

combined_primary <- bind_rows(primary_results) %>%
  arrange(match(Cohort, names(expr_file_map)), display_order)
combined_qa <- bind_rows(qa_rows)

if (
  nrow(combined_primary) != 36 ||
    anyDuplicated(combined_primary[c("Cohort", "feature_id")]) > 0 ||
    !all(combined_qa$Pathway_N == 12)
) {
  stop("三队列LM12 GSEA最终QA失败。")
}

write_tsv(combined_primary, file.path(out_dir, "external_lm_imm_12_gene_level_gsea_results.tsv"))
write_tsv(combined_qa, file.path(out_dir, "external_lm_imm_12_gene_level_gsea_qa.tsv"))
if (length(sensitivity_results) > 0) {
  write_tsv(bind_rows(sensitivity_results), file.path(out_dir, "external_lm_imm_12_gene_level_gsea_sensitivity_results.tsv"))
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "external_lm_imm_12_gene_level_gsea_sessionInfo.txt"))

message("外部LM12 gene-level GSEA完成。")
message("主分析结果: ", nrow(combined_primary), " rows")
message(paste0(combined_qa$Cohort, "/", combined_qa$analysis_set, ": n=", combined_qa$Sample_N, ", direction=", combined_qa$Direction_Consistent_N, "/12, FDR<0.05=", combined_qa$FDR_lt_0_05_N, collapse = "\n"))
message("输出目录: ", out_dir)
}
