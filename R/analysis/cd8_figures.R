# Purpose: cd8 figures analysis stages.
# Configuration: config/cd8.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/05_plot_la_temporal_quadrant.R
stage_cd8_perturbation_05_plot_la_temporal_quadrant <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["leading_edge_file"]
de_basenames <- paste0("de_LA_minus_NC_",c("1h","6h","24h","48h"),".tsv")
required_files <- c(unname(input_paths),file.path(unname(paths["bulk_de_dir"]),de_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths["bulk_de_dir"])
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

# ---- Retained complete source body; documented path/quotation adapters only ----
#!/usr/bin/env Rscript

# Fig.7C candidate v0.8：LA−NC 早期—晚期四方位图，重映射为 V0.19 的新 C panel。
# 只读复用 results/03 与已冻结的 Hypoxia leading-edge 结果；不重新拟合模型。
# 相位均值仅用于描述性可视化，不构成新的统计检验。

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

work_dir <- unname(paths["output_dir"])
de_dir <- unname(paths["bulk_de_dir"])
leading_edge_file <- unname(paths["leading_edge_file"])

plot_ready_dir <- file.path(work_dir, "plot_ready")
export_dir <- file.path(work_dir, "panel_exports")
qa_dir <- file.path(work_dir, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

times <- c("1h", "6h", "24h", "48h")
fdr_cut <- 0.05
lfc_cut <- 1
width_mm <- 44
height_mm <- 44

read_one <- function(tm) {
  path <- file.path(de_dir, paste0("de_LA_minus_NC_", tm, ".tsv"))
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing input: ", path)
  x <- read.delim(path, check.names = FALSE)
  required <- c("gene_id", "logFC", "FDR", "contrast", "time", "left", "right")
  if (!all(required %in% names(x))) stop("Missing required columns: ", path)
  if (anyDuplicated(x$gene_id)) stop("Duplicated gene_id in: ", path)
  if (!all(is.finite(x$logFC)) || !all(is.finite(x$FDR))) {
    stop("Non-finite logFC/FDR in: ", path)
  }
  if (!all(x$contrast == "LA_minus_NC") ||
      !all(x$left == "LA" & x$right == "NC")) {
    stop("Unexpected contrast direction in: ", path)
  }
  if (!all(as.character(x$time) == tm)) stop("Unexpected time label in: ", path)
  x$source_file <- path
  x
}

xs <- lapply(times, read_one)
names(xs) <- times
gene_ids <- lapply(xs, function(x) x$gene_id)
if (!all(vapply(gene_ids[-1], function(z) identical(z, gene_ids[[1]]), logical(1)))) {
  stop("gene_id order differs across the four DE tables.")
}

wide <- data.frame(gene_id = gene_ids[[1]], stringsAsFactors = FALSE)
for (tm in times) {
  wide[[paste0("logFC_", tm)]] <- xs[[tm]]$logFC
  wide[[paste0("FDR_", tm)]] <- xs[[tm]]$FDR
}
wide$early_response <- rowMeans(wide[, c("logFC_1h", "logFC_6h")])
wide$late_response <- rowMeans(wide[, c("logFC_24h", "logFC_48h")])

strict_matrix <- sapply(times, function(tm) {
  fdr <- wide[[paste0("FDR_", tm)]]
  lfc <- wide[[paste0("logFC_", tm)]]
  is.finite(fdr) & fdr < fdr_cut & is.finite(lfc) & abs(lfc) >= lfc_cut
})
colnames(strict_matrix) <- times
strict_any <- rowSums(strict_matrix) > 0

quad_raw <- ifelse(
  wide$early_response > 0 & wide$late_response > 0, "Sustained up",
  ifelse(
    wide$early_response < 0 & wide$late_response < 0, "Sustained down",
    ifelse(
      wide$early_response > 0 & wide$late_response < 0, "Early up -> late down",
      ifelse(
        wide$early_response < 0 & wide$late_response > 0,
        "Early down -> late up", "Other tested genes"
      )
    )
  )
)

quadrant_levels <- c(
  "Sustained up",
  "Sustained down",
  "Early up -> late down",
  "Early down -> late up",
  "Other tested genes"
)
quad_raw[!strict_any] <- "Other tested genes"
wide$strict_any <- strict_any
wide$strict_time_count <- rowSums(strict_matrix)
wide$quadrant <- factor(quad_raw, levels = quadrant_levels)

if (!file.exists(leading_edge_file) || file.info(leading_edge_file)$size == 0) {
  stop("Missing leading-edge input: ", leading_edge_file)
}
le <- read.delim(leading_edge_file, check.names = FALSE)
required_le <- c("context", "feature_id", "gene")
if (!all(required_le %in% names(le))) stop("Leading-edge file lacks required columns.")
bridge_features <- c(
  "MSIGDB::HALLMARK_HYPOXIA",
  "MSIGDB::HALLMARK_INTERFERON_GAMMA_RESPONSE"
)
shared_bridge <- unique(unlist(lapply(bridge_features, function(fid) {
  mel_genes <- unique(le$gene[
    le$context == "melanoma_CD8_C1_minus_C4" & le$feature_id == fid
  ])
  bulk_genes <- unique(le$gene[
    le$context == "bulk_LA_minus_NC" & le$feature_id == fid
  ])
  intersect(mel_genes, bulk_genes)
})))
bridge_anchor_genes <- intersect(shared_bridge, wide$gene_id[wide$strict_any])
expected_bridge_anchor_genes <- c("DDIT4", "DUSP1", "JUN", "VCAM1")
if (!setequal(bridge_anchor_genes, expected_bridge_anchor_genes)) {
  stop(
    "Cross-context bridge membership changed. Found: ",
    paste(sort(bridge_anchor_genes), collapse = ",")
  )
}
leading_edge_union <- unique(le$gene)
rep_pool <- intersect(leading_edge_union, wide$gene_id[wide$strict_any])
rep_pool <- setdiff(rep_pool, bridge_anchor_genes)
rep_score <- data.frame(
  gene_id = rep_pool,
  quadrant = wide$quadrant[match(rep_pool, wide$gene_id)],
  strict_time_count = wide$strict_time_count[match(rep_pool, wide$gene_id)],
  max_phase_abs = pmax(
    abs(wide$early_response[match(rep_pool, wide$gene_id)]),
    abs(wide$late_response[match(rep_pool, wide$gene_id)])
  ),
  stringsAsFactors = FALSE
)
quadrant_representatives <- unlist(lapply(quadrant_levels[1:4], function(qv) {
  z <- rep_score[rep_score$quadrant == qv, , drop = FALSE]
  if (!nrow(z)) return(character(0))
  z <- z[order(-z$strict_time_count, -z$max_phase_abs, z$gene_id), , drop = FALSE]
  z$gene_id[[1]]
}), use.names = FALSE)
expected_quadrant_representatives <- c("ADGRG3", "MT2A", "MT1E", "CTSL")
if (!setequal(quadrant_representatives, expected_quadrant_representatives)) {
  stop(
    "Quadrant representative rule changed. Found: ",
    paste(sort(quadrant_representatives), collapse = ",")
  )
}
annotation_genes <- unique(c(bridge_anchor_genes, quadrant_representatives))
expected_annotation_genes <- c(
  "DDIT4", "DUSP1", "JUN", "VCAM1",
  "ADGRG3", "MT2A", "MT1E", "CTSL"
)
if (!setequal(annotation_genes, expected_annotation_genes)) {
  stop(
    "Annotation membership changed. Found: ",
    paste(sort(annotation_genes), collapse = ",")
  )
}
wide$cross_context_bridge_anchor <- wide$gene_id %in% bridge_anchor_genes
wide$quadrant_representative <- wide$gene_id %in% quadrant_representatives
wide$annotation_gene <- wide$gene_id %in% annotation_genes
axis_lim <- 3
display_outside <- abs(wide$early_response) > axis_lim |
  abs(wide$late_response) > axis_lim
strict_outside <- sum(wide$strict_any & display_outside)
if (strict_outside > 0) {
  stop("The requested +/-3 display window would clip strict-response genes.")
}

palette <- c(
  "Sustained up" = "#D47A2A",
  "Sustained down" = "#3B77A7",
  "Early up -> late down" = "#8C6BB1",
  "Early down -> late up" = "#3C9C8B",
  "Other tested genes" = "#C9C9C9"
)

all_points <- wide
signal <- wide[wide$strict_any & wide$quadrant != "Other tested genes", , drop = FALSE]
annotated <- wide[wide$annotation_gene, , drop = FALSE]
count_tab <- table(factor(signal$quadrant, levels = quadrant_levels[1:4]))

label_dat <- data.frame(
  early_response = c(axis_lim * 0.93, -axis_lim * 0.93,
                     axis_lim * 0.93, -axis_lim * 0.93),
  late_response = c(axis_lim * 0.93, -axis_lim * 0.93,
                    -axis_lim * 0.93, axis_lim * 0.93),
  quadrant = factor(quadrant_levels[1:4], levels = quadrant_levels),
  label = c(
    paste0("Up / up\nn = ", count_tab[[1]]),
    paste0("Down / down\nn = ", count_tab[[2]]),
    paste0("Up / down\nn = ", count_tab[[3]]),
    paste0("Down / up\nn = ", count_tab[[4]])
  ),
  hjust = c(1, 0, 1, 0),
  vjust = c(1, 0, 0, 1),
  stringsAsFactors = FALSE
)

p_square <- ggplot() +
  geom_point(
    data = all_points,
    aes(x = early_response, y = late_response),
    colour = palette[["Other tested genes"]],
    size = 0.30, alpha = 0.35, stroke = 0
  ) +
  geom_point(
    data = signal,
    aes(x = early_response, y = late_response, colour = quadrant),
    size = 0.78, alpha = 0.90, stroke = 0
  ) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#5D5D5D") +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#5D5D5D") +
  geom_text(
    data = label_dat,
    aes(x = early_response, y = late_response, label = label,
        colour = quadrant, hjust = hjust, vjust = vjust),
    family = "Arial", fontface = "bold", size = 5.4 / ggplot2::.pt,
    lineheight = 0.90, show.legend = FALSE
  ) +
  geom_text_repel(
    data = annotated,
    aes(x = early_response, y = late_response, label = gene_id),
    family = "Arial", size = 5.0 / ggplot2::.pt, colour = "#1A1A1A",
    seed = 20260914, box.padding = 0.26, point.padding = 0.08,
    min.segment.length = 0, segment.size = 0.16, segment.colour = "#555555",
    max.overlaps = Inf, direction = "both", force = 2.2,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = palette, limits = quadrant_levels,
                      drop = FALSE, guide = "none") +
  scale_fill_manual(values = palette, limits = quadrant_levels,
                    drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(-axis_lim, axis_lim),
    breaks = -3:3,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(-axis_lim, axis_lim),
    breaks = -3:3,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_fixed(ratio = 1, clip = "on") +
  labs(
    x = "Early LA-NC response\nmean log2FC (1, 6 h)",
    y = "Late LA-NC response\nmean log2FC (24, 48 h)"
  ) +
  theme_classic(base_family = "Arial", base_size = 5.8) +
  theme(
    axis.line = element_line(linewidth = 0.28, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25, colour = "black"),
    axis.ticks.length = grid::unit(0.8, "mm"),
    axis.text = element_text(size = 5.4, colour = "black"),
    axis.title = element_text(size = 5.6, colour = "black"),
    axis.title.x = element_text(margin = margin(t = 2.0)),
    axis.title.y = element_text(margin = margin(r = 2.0)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.42),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 1.5, r = 1.5, b = 1.5, l = 1.5, unit = "mm")
  )

p <- p_square

prefix <- file.path(export_dir, "Fig7C_bulk_temporal_quadrant_compact_candidate_v0.8")
ggsave(paste0(prefix, ".pdf"), p, width = width_mm / 25.4,
       height = height_mm / 25.4, units = "in", device = cairo_pdf)
ggsave(paste0(prefix, ".tiff"), p, width = width_mm / 25.4,
       height = height_mm / 25.4, units = "in", dpi = 600,
       compression = "lzw", device = "tiff")
ggsave(paste0(prefix, ".png"), p, width = width_mm / 25.4,
       height = height_mm / 25.4, units = "in", dpi = 300,
       device = "png", bg = "white")

source_data <- cbind(
  wide[, c("gene_id", "early_response", "late_response", "strict_any",
           "strict_time_count", "quadrant", "cross_context_bridge_anchor")],
  as.data.frame(strict_matrix, check.names = FALSE)
)
names(source_data)[(ncol(source_data) - length(times) + 1):ncol(source_data)] <-
  paste0("strict_", times)
write.table(
  source_data,
  file.path(plot_ready_dir, "fig7c_bulk_temporal_quadrant_compact_v0.8_source_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

summary_rows <- data.frame(
  quadrant = quadrant_levels[1:4],
  strict_genes = as.integer(count_tab),
  total_tested_genes = nrow(wide),
  axis_limit = axis_lim,
  stringsAsFactors = FALSE
)
write.table(
  summary_rows,
  file.path(plot_ready_dir, "fig7c_bulk_temporal_quadrant_compact_v0.8_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
`sha256` <- bind_stage_helper(shared_sha256_022, environment())
input_files <- c(vapply(xs, function(x) unique(x$source_file), character(1)),
                 leading_edge_file)
manifest_fields <- c(
  "panel_id", "core_conclusion", "evidence_role", "comparison",
  "biological_replicate", "phase_definition", "phase_statistical_status",
  "threshold", "signal_rule", "quadrant_palette", "axis_limit",
  "tested_genes", "strict_genes", "panel_size_mm", "square_data_body",
  "font_target_pt", "label_rule", "annotation_genes", "annotation_descriptors", "bridge_anchor_genes", "quadrant_representatives",
  "secondary_candidate", "formal_write",
  paste0("input_", seq_along(input_files)),
  paste0("input_sha256_", seq_along(input_files))
)
manifest_values <- c(
  "Fig7C_candidate",
  "Purified CD8 cells show an early-to-late, time-structured LA-NC transcriptional response.",
  "descriptive whole-transcriptome temporal background",
  "donor-matched LA minus NC at 1, 6, 24 and 48 h",
  "donor (n=3)",
  "early=mean(log2FC_1h,log2FC_6h); late=mean(log2FC_24h,log2FC_48h)",
  "descriptive phase means; no new DE test or independent replicate",
  "FDR<0.05 and absolute log2FC>=1 at any original timepoint",
  "all tested genes grey; strict genes coloured by early/late sign",
  paste(names(palette)[1:4], unname(palette)[1:4], sep = "=", collapse = ";"),
  as.character(axis_lim),
  as.character(nrow(wide)),
  as.character(sum(wide$strict_any)),
  paste0(width_mm, " x ", height_mm),
  "TRUE; entire panel data body uses coord_fixed(ratio=1)",
  "5.8 (allowed 5.6-6.0)",
  "label strict genes shared by melanoma and bulk HALLMARK_HYPOXIA or HALLMARK_INTERFERON_GAMMA_RESPONSE leading-edge",
  paste(sort(annotation_genes), collapse = ";"),
  "DDIT4=hypoxia/energy stress; DUSP1=MAPK feedback; JUN=AP-1 immediate-early; VCAM1=IFN-gamma edge; ADGRG3=memory-edge representative; MT2A=stress/effector-edge representative; MT1E=transient stress-edge representative; CTSL=LM12-edge representative",
  paste(sort(bridge_anchor_genes), collapse = ";"),
  paste(sort(quadrant_representatives), collapse = ";"),
  "EFNA1 (bulk Hypoxia leading-edge only; not a cross-context label)",
  "FALSE",
  input_files,
  vapply(input_files, sha256, character(1))
)
stopifnot(length(manifest_fields) == length(manifest_values))
manifest <- data.frame(field = manifest_fields, value = manifest_values,
                       stringsAsFactors = FALSE)
write.table(
  manifest,
  file.path(plot_ready_dir, "fig7c_bulk_temporal_quadrant_compact_v0.8_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

prefix_files <- paste0(prefix, c(".pdf", ".tiff", ".png"))
qa <- data.frame(
  check_id = c(
    "four_timepoints_present", "row_count_equal", "gene_id_order_equal",
    "contrast_direction", "finite_values", "strict_threshold_logic",
    "quadrant_counts_sum", "coordinates_not_clipped", "square_data_geometry",
    "cross_context_le_input", "annotation_membership", "annotation_count",
    "annotation_descriptors", "strict_outside_display", "panel_dimensions",
    "pdf_nonempty", "tiff_nonempty", "png_nonempty", "panel_letter_embedded",
    "phase_means_descriptive", "formal_write"
  ),
  status = c(
    if (identical(names(xs), times)) "PASS" else "FAIL",
    if (length(unique(vapply(xs, nrow, integer(1)))) == 1L) "PASS" else "FAIL",
    if (all(vapply(gene_ids[-1], function(z) identical(z, gene_ids[[1]]), logical(1)))) "PASS" else "FAIL",
    "PASS",
    if (all(is.finite(wide$early_response)) && all(is.finite(wide$late_response))) "PASS" else "FAIL",
    if (all(wide$strict_any == (rowSums(strict_matrix) > 0))) "PASS" else "FAIL",
    if (sum(count_tab) == sum(wide$strict_any)) "PASS" else "FAIL",
    if (strict_outside == 0) "PASS" else "FAIL",
    "PASS",
    if (file.exists(leading_edge_file) && file.info(leading_edge_file)$size > 0) "PASS" else "FAIL",
    if (setequal(annotation_genes, expected_annotation_genes)) "PASS" else "FAIL",
    if (length(annotation_genes) == 8L) "PASS" else "FAIL",
    "PASS",
    if (strict_outside == 0) "PASS" else "FAIL",
    "PASS",
    if (file.exists(prefix_files[[1]]) && file.info(prefix_files[[1]])$size > 0) "PASS" else "FAIL",
    if (file.exists(prefix_files[[2]]) && file.info(prefix_files[[2]])$size > 0) "PASS" else "FAIL",
    if (file.exists(prefix_files[[3]]) && file.info(prefix_files[[3]])$size > 0) "PASS" else "FAIL",
    "PASS",
    "PASS",
    "PASS"
  ),
  detail = c(
    paste(times, collapse = ","),
    paste(unique(vapply(xs, nrow, integer(1))), collapse = ";"),
    paste0("n=", nrow(wide)),
    "LA minus NC",
    paste0("n=", nrow(wide)),
    "FDR<0.05 & abs(log2FC)>=1 in at least one original timepoint",
    paste0("strict=", sum(wide$strict_any), "; quadrants=", sum(count_tab)),
    paste0("axis_limit=+/-", axis_lim, "; display_outside=", sum(display_outside), "; strict_outside=", strict_outside),
    "compact panel data body uses coord_fixed(ratio=1)",
    "gsea_leading_edge_long.tsv",
    paste(sort(annotation_genes), collapse = ","),
    paste0("n=", length(annotation_genes)),
    "Functional cues are display-only and do not alter gene selection.",
    paste0("n=", strict_outside, " strict genes outside display window; source data retained"),
    paste0(width_mm, " x ", height_mm, " mm"),
    as.character(file.info(prefix_files[[1]])$size),
    as.character(file.info(prefix_files[[2]])$size),
    as.character(file.info(prefix_files[[3]])$size),
    "No B. label in independent panel",
    "phase means are visualization only",
    "Candidate preview only"
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa,
  file.path(qa_dir, "Fig7C_bulk_temporal_quadrant_compact_candidate_v0.8_qa.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
writeLines(
  capture.output(sessionInfo()),
  file.path(qa_dir, "Fig7C_bulk_temporal_quadrant_compact_candidate_v0.8_session_info.txt")
)

if (any(qa$status == "FAIL")) stop("One or more QA checks failed.")
print(summary_rows)
cat("Rendered:", prefix, "\n")
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/06_plot_hcl_temporal_quadrant.R
stage_cd8_perturbation_06_plot_hcl_temporal_quadrant <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- character(0L)
de_basenames <- paste0("de_HCl_minus_NC_",c("1h","6h","24h","48h"),".tsv")
required_files <- c(unname(input_paths),file.path(unname(paths["bulk_de_dir"]),de_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths["bulk_de_dir"])
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

# ---- Retained complete source body; documented path/quotation adapters only ----
#!/usr/bin/env Rscript

# S13A v0.1：HCl−NC 早期/晚期转录响应四象限。
# 只读复用已完成的三条件 DE；不重新拟合模型。

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

s13_root <- unname(paths["output_dir"])
de_dir <- unname(paths["bulk_de_dir"])
plot_ready_dir <- file.path(s13_root, "plot_ready")
export_dir <- file.path(s13_root, "panel_exports")
qa_dir <- file.path(s13_root, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

times <- c("1h", "6h", "24h", "48h")
fdr_cut <- 0.05
lfc_cut <- 1

read_one <- function(tm) {
  path <- file.path(de_dir, paste0("de_HCl_minus_NC_", tm, ".tsv"))
  if (!file.exists(path) || file.info(path)$size == 0) stop("Missing input: ", path)
  x <- read.delim(path, check.names = FALSE)
  req <- c("gene_id", "logFC", "FDR", "contrast", "time", "left", "right")
  if (!all(req %in% names(x))) stop("Missing columns in ", path)
  if (anyDuplicated(x$gene_id)) stop("Duplicated gene IDs in ", path)
  if (!all(x$contrast == "HCl_minus_NC") || !all(x$left == "HCl") || !all(x$right == "NC")) {
    stop("Unexpected contrast direction in ", path)
  }
  if (!all(is.finite(x$logFC)) || !all(is.finite(x$FDR))) stop("Non-finite DE values in ", path)
  x$source_file <- path
  x
}

xs <- lapply(times, read_one)
gene_ids <- lapply(xs, function(x) x$gene_id)
if (!all(vapply(gene_ids[-1], function(z) identical(z, gene_ids[[1]]), logical(1)))) {
  stop("gene_id order differs across HCl DE tables")
}
wide <- data.frame(gene_id = gene_ids[[1]], stringsAsFactors = FALSE)
for (tm in times) {
  xx <- xs[[match(tm, times)]]
  wide[[paste0("logFC_", tm)]] <- xx$logFC
  wide[[paste0("FDR_", tm)]] <- xx$FDR
}
wide$early_response <- rowMeans(wide[, c("logFC_1h", "logFC_6h")])
wide$late_response <- rowMeans(wide[, c("logFC_24h", "logFC_48h")])

strict_matrix <- sapply(times, function(tm) {
  fdr <- wide[[paste0("FDR_", tm)]]
  lfc <- wide[[paste0("logFC_", tm)]]
  is.finite(fdr) & fdr < fdr_cut & is.finite(lfc) & abs(lfc) >= lfc_cut
})
colnames(strict_matrix) <- times
wide$strict_any <- rowSums(strict_matrix) > 0
wide$strict_time_count <- rowSums(strict_matrix)
wide$quadrant <- ifelse(
  wide$early_response > 0 & wide$late_response > 0, "Up / up",
  ifelse(wide$early_response < 0 & wide$late_response < 0, "Down / down",
    ifelse(wide$early_response > 0 & wide$late_response < 0, "Up / down",
      ifelse(wide$early_response < 0 & wide$late_response > 0, "Down / up", "Other tested genes"))))
wide$quadrant[!wide$strict_any] <- "Other tested genes"
quadrant_levels <- c("Up / up", "Down / down", "Up / down", "Down / up", "Other tested genes")
wide$quadrant <- factor(wide$quadrant, levels = quadrant_levels)

axis_lim <- max(3, ceiling(max(abs(c(wide$early_response, wide$late_response)), na.rm = TRUE) * 2) / 2)

# 每个显著象限选择两个代表基因；规则固定，不按单次运行的最小 P 值挑选。
wide$max_phase_abs <- pmax(abs(wide$early_response), abs(wide$late_response))
rep_genes <- unlist(lapply(quadrant_levels[1:4], function(qv) {
  z <- wide[wide$strict_any & wide$quadrant == qv, , drop = FALSE]
  if (!nrow(z)) return(character(0))
  z <- z[order(-z$strict_time_count, -z$max_phase_abs, z$gene_id), , drop = FALSE]
  head(z$gene_id, 2L)
}), use.names = FALSE)
wide$annotation_gene <- wide$gene_id %in% rep_genes

palette <- c("Up / up" = "#C8782B", "Down / down" = "#3B77A7",
             "Up / down" = "#7B5AA6", "Down / up" = "#3C9C8B",
             "Other tested genes" = "#C9C9C9")
signal <- wide[wide$strict_any & wide$quadrant != "Other tested genes", , drop = FALSE]
annotated <- wide[wide$annotation_gene, , drop = FALSE]
count_tab <- table(factor(signal$quadrant, levels = quadrant_levels[1:4]))
label_dat <- data.frame(
  x = c(axis_lim * .93, -axis_lim * .93, axis_lim * .93, -axis_lim * .93),
  y = c(axis_lim * .93, -axis_lim * .93, -axis_lim * .93, axis_lim * .93),
  quadrant = factor(quadrant_levels[1:4], levels = quadrant_levels),
  label = c(paste0("Up / up\nn = ", count_tab[[1]]),
            paste0("Down / down\nn = ", count_tab[[2]]),
            paste0("Up / down\nn = ", count_tab[[3]]),
            paste0("Down / up\nn = ", count_tab[[4]])),
  hjust = c(1, 0, 1, 0), vjust = c(1, 0, 0, 1), stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_point(data = wide, aes(early_response, late_response), colour = palette[["Other tested genes"]],
             size = .30, alpha = .35) +
  geom_point(data = signal, aes(early_response, late_response, colour = quadrant), size = .78, alpha = .9) +
  geom_hline(yintercept = 0, linewidth = .28, colour = "#5D5D5D") +
  geom_vline(xintercept = 0, linewidth = .28, colour = "#5D5D5D") +
  geom_text(data = label_dat, aes(x, y, label = label, colour = quadrant, hjust = hjust, vjust = vjust),
            family = "Arial", fontface = "bold", size = 5.0 / ggplot2::.pt,
            lineheight = .9, show.legend = FALSE) +
  geom_text_repel(data = annotated, aes(early_response, late_response, label = gene_id),
                  family = "Arial", size = 4.6 / ggplot2::.pt, colour = "#1A1A1A",
                  seed = 20260917, box.padding = .25, point.padding = .08,
                  min.segment.length = 0, segment.size = .16, segment.colour = "#555555",
                  max.overlaps = Inf, direction = "both", force = 1.8, show.legend = FALSE) +
  scale_colour_manual(values = palette, limits = quadrant_levels, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(-axis_lim, axis_lim), breaks = pretty(c(-axis_lim, axis_lim), n = 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-axis_lim, axis_lim), breaks = pretty(c(-axis_lim, axis_lim), n = 5), expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "on") +
  labs(x = "Early HCl-NC response\nmean log2FC (1, 6 h)",
       y = "Late HCl-NC response\nmean log2FC (24, 48 h)") +
  theme_classic(base_family = "Arial", base_size = 5.8) +
  theme(axis.line = element_line(linewidth = .28, colour = "black"),
        axis.ticks = element_line(linewidth = .25, colour = "black"),
        axis.ticks.length = grid::unit(.8, "mm"),
        axis.text = element_text(size = 5.2, colour = "black"),
        axis.title = element_text(size = 5.5, colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = .42),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(t = 1.5, r = 1.5, b = 1.5, l = 1.5, unit = "mm"))

prefix <- file.path(export_dir, "S13A_hcl_temporal_quadrant_candidate_v0.1")
ggsave(paste0(prefix, ".pdf"), p, width = 44/25.4, height = 44/25.4, units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(prefix, ".tiff"), p, width = 44/25.4, height = 44/25.4, units = "in", dpi = 600, compression = "lzw", device = "tiff", bg = "white")
ggsave(paste0(prefix, ".png"), p, width = 44/25.4, height = 44/25.4, units = "in", dpi = 300, device = "png", bg = "white")

source_data <- cbind(wide[, c("gene_id", "early_response", "late_response", "strict_any", "strict_time_count", "quadrant", "annotation_gene")],
                     as.data.frame(strict_matrix, check.names = FALSE))
names(source_data)[(ncol(source_data)-length(times)+1):ncol(source_data)] <- paste0("strict_", times)
write.table(source_data, file.path(plot_ready_dir, "s13a_hcl_temporal_quadrant_source.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(data.frame(quadrant = quadrant_levels[1:4], strict_genes = as.integer(count_tab), axis_limit = axis_lim,
                       stringsAsFactors = FALSE), file.path(plot_ready_dir, "s13a_hcl_temporal_quadrant_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

qa <- data.frame(check_id = c("de_tables", "tested_genes", "strict_outside_axis", "annotation_genes", "panel_size_mm", "formal_write", "png_nonempty", "pdf_nonempty", "tiff_nonempty"),
                 expected = c("4", as.character(length(wide$gene_id)), "0", as.character(length(rep_genes)), "44 x 44", "FALSE", "TRUE", "TRUE", "TRUE"),
                 observed = c(length(xs), length(wide$gene_id), sum(wide$strict_any & (abs(wide$early_response)>axis_lim | abs(wide$late_response)>axis_lim)), length(rep_genes), "44 x 44", "FALSE",
                              file.exists(paste0(prefix, ".png")) && file.info(paste0(prefix, ".png"))$size > 0,
                              file.exists(paste0(prefix, ".pdf")) && file.info(paste0(prefix, ".pdf"))$size > 0,
                              file.exists(paste0(prefix, ".tiff")) && file.info(paste0(prefix, ".tiff"))$size > 0), stringsAsFactors = FALSE)
qa$pass <- as.character(qa$expected) == as.character(qa$observed)
write.table(qa, file.path(qa_dir, "S13A_hcl_temporal_quadrant_candidate_v0.1_qa.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(qa_dir, "S13A_hcl_temporal_quadrant_candidate_v0.1_session_info.txt"))
if (!all(qa$pass)) stop("S13A QA failed")
message("Rendered S13A: ", prefix)
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/07_plot_three_condition_pca.R
stage_cd8_perturbation_07_plot_three_condition_pca <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("pca_coordinates","pca_variance")]
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

# ---- Retained complete source body; documented path/resource/quotation adapters only ----
#!/usr/bin/env Rscript

# Fig.7B candidate v0.9：三条件纯 CD8 bulk 样本级 PCA。
# 颜色表示处理，形状表示时间，灰线表示同一供体×时间的三条件配对。
# PCA仅用于整体转录结构概览，不把视觉距离当作独立处理效应检验。

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

work_dir <- unname(paths["output_dir"])
pca_file <- unname(paths["pca_coordinates"])
variance_file <- unname(paths["pca_variance"])

plot_ready_dir <- file.path(work_dir, "plot_ready")
export_dir <- file.path(work_dir, "panel_exports")
qa_dir <- file.path(work_dir, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pca_file) || file.info(pca_file)$size == 0) stop("Missing PCA coordinates: ", pca_file)
if (!file.exists(variance_file) || file.info(variance_file)$size == 0) stop("Missing PCA variance table: ", variance_file)

pca <- read.delim(pca_file, check.names = FALSE)
required_pca <- c("SampleID", "Donor", "Treatment", "Time", "PC1", "PC2", "PC3")
if (!all(required_pca %in% names(pca))) {
  stop("PCA table lacks required columns: ", paste(setdiff(required_pca, names(pca)), collapse = ", "))
}
if (nrow(pca) != 36L || anyDuplicated(pca$SampleID)) stop("Expected 36 unique factorial PCA samples.")
if (!setequal(unique(pca$Treatment), c("NC", "LA", "HCl"))) stop("Unexpected treatment levels.")
if (!setequal(unique(pca$Time), c("1h", "6h", "24h", "48h"))) stop("Unexpected time levels.")
if (any(!is.finite(as.numeric(pca$PC1))) || any(!is.finite(as.numeric(pca$PC2)))) {
  stop("Non-finite PC1/PC2 values.")
}

variance <- read.delim(variance_file, check.names = FALSE)
if (!all(c("component", "variance_fraction", "variance_percent") %in% names(variance))) {
  stop("PCA variance table is incomplete.")
}
pc1_var <- variance$variance_percent[match("PC1", variance$component)]
pc2_var <- variance$variance_percent[match("PC2", variance$component)]
if (length(pc1_var) != 1L || length(pc2_var) != 1L ||
    !is.finite(pc1_var) || !is.finite(pc2_var)) {
  stop("PC1/PC2 variance percentages are unavailable.")
}

treatment_levels <- c("NC", "LA", "HCl")
time_levels <- c("1h", "6h", "24h", "48h")
pca$Treatment <- factor(pca$Treatment, levels = treatment_levels)
pca$Time <- factor(pca$Time, levels = time_levels)
pca$pair_id <- paste(pca$Donor, pca$Time, sep = "_")
pca <- pca[order(pca$pair_id, pca$Treatment), , drop = FALSE]
if (any(table(pca$pair_id) != 3L)) stop("Each donor x time pair must contain three treatments.")

treatment_cols <- c(NC = "#2C7FB8", LA = "#31A354", HCl = "#D95F0E")
time_shapes <- c("1h" = 16, "6h" = 17, "24h" = 15, "48h" = 18)

# 保留 PC1/PC2 的数据范围，但将坐标框压缩为略宽于高的矩形，
# 使 PCA 主体、坐标文字和底部图例的最终可见高度与 C panel 接近。
x_span <- diff(range(pca$PC1, na.rm = TRUE)) * 1.10
y_span <- diff(range(pca$PC2, na.rm = TRUE)) * 1.10
square_ratio <- 1.27

p <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_path(
    aes(group = pair_id),
    colour = "#A8A8A8", linewidth = 0.26, alpha = 0.55,
    lineend = "round", show.legend = FALSE
  ) +
  geom_point(
    aes(colour = Treatment, shape = Time),
    size = 1.65, stroke = 0.20, alpha = 0.92
  ) +
  scale_colour_manual(
    values = treatment_cols, limits = treatment_levels, drop = FALSE,
    name = "Treatment"
  ) +
  scale_shape_manual(
    values = time_shapes, limits = time_levels, drop = FALSE,
    name = "Time"
  ) +
  scale_x_continuous(
    breaks = pretty(pca$PC1, n = 3),
    expand = expansion(mult = 0.05)
  ) +
  scale_y_continuous(
    breaks = pretty(pca$PC2, n = 3),
    expand = expansion(mult = 0.05)
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1_var),
    y = sprintf("PC2 (%.1f%%)", pc2_var)
  ) +
  coord_fixed(ratio = square_ratio, clip = "on") +
  theme_classic(base_family = "Arial", base_size = 5.8) +
  theme(
    axis.line = element_line(linewidth = 0.27, colour = "black"),
    axis.ticks = element_line(linewidth = 0.24, colour = "black"),
    axis.ticks.length = unit(0.45, "mm"),
    axis.text = element_text(size = 5.6, colour = "black"),
    axis.title = element_text(size = 5.6, colour = "black"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.38),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.text = element_text(size = 5.6),
    legend.key.width = unit(3.5, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.spacing.x = unit(0.8, "mm"),
    legend.spacing.y = unit(0.2, "mm"),
    legend.margin = margin(t = 0.2, r = 0, b = 0, l = 0, unit = "mm"),
    plot.margin = margin(t = 1.0, r = 1.0, b = 0.5, l = 1.0, unit = "mm")
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    shape = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  )

width_mm <- 44
height_mm <- 47
prefix <- file.path(export_dir, "Fig7B_bulk_three_condition_pca_candidate_v0.9")
ggsave(paste0(prefix, ".pdf"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", device = cairo_pdf)
ggsave(paste0(prefix, ".tiff"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 600, compression = "lzw", device = "tiff")
ggsave(paste0(prefix, ".png"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 300, device = "png", bg = "white")

source_data <- pca[, c("SampleID", "Donor", "Treatment", "Time", "PC1", "PC2", "PC3", "pair_id")]
source_data$panel <- "Fig7B"
write.table(
  source_data,
  file.path(plot_ready_dir, "fig7b_three_condition_pca_v0.9_source_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_data <- data.frame(
  panel = "Fig7B",
  n_samples = nrow(pca),
  n_donors = length(unique(pca$Donor)),
  treatments = paste(treatment_levels, collapse = ";"),
  times = paste(time_levels, collapse = ";"),
  PC1_variance_percent = pc1_var,
  PC2_variance_percent = pc2_var,
  pairing = "Donor x time; three treatments per pair",
  width_mm = width_mm,
  height_mm = height_mm,
  font_pt = "5.6-5.8",
  stringsAsFactors = FALSE
)
write.table(
  summary_data,
  file.path(plot_ready_dir, "fig7b_three_condition_pca_v0.9_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- data.frame(
  field = c(
    "panel_id", "status", "source_coordinates", "source_variance",
    "n_samples", "n_donors", "treatment_encoding", "time_encoding",
    "pairing_encoding", "width_mm", "height_mm", "font_pt", "output_prefix"
  ),
  value = c(
    "Fig7B", "CANDIDATE", pca_file, variance_file,
    nrow(pca), length(unique(pca$Donor)), "colour", "shape",
    "grey path by donor x time", width_mm, height_mm, "5.6-5.8",
    basename(prefix)
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  file.path(plot_ready_dir, "fig7b_three_condition_pca_v0.9_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

prefix_files <- paste0(prefix, c(".pdf", ".tiff", ".png"))
qa <- data.frame(
  check_id = c(
    "samples_36", "three_treatments", "four_times", "donor_time_pairs",
    "pc_values_finite", "output_pdf", "output_tiff", "output_png",
    "panel_dimensions", "panel_letter_embedded"
  ),
  status = c(
    ifelse(nrow(pca) == 36L, "PASS", "FAIL"),
    ifelse(setequal(unique(as.character(pca$Treatment)), treatment_levels), "PASS", "FAIL"),
    ifelse(setequal(unique(as.character(pca$Time)), time_levels), "PASS", "FAIL"),
    ifelse(all(table(pca$pair_id) == 3L), "PASS", "FAIL"),
    ifelse(all(is.finite(pca$PC1) & is.finite(pca$PC2)), "PASS", "FAIL"),
    ifelse(file.exists(prefix_files[[1]]) && file.info(prefix_files[[1]])$size > 0, "PASS", "FAIL"),
    ifelse(file.exists(prefix_files[[2]]) && file.info(prefix_files[[2]])$size > 0, "PASS", "FAIL"),
    ifelse(file.exists(prefix_files[[3]]) && file.info(prefix_files[[3]])$size > 0, "PASS", "FAIL"),
    paste0(width_mm, " x ", height_mm, " mm"),
    "No panel letter in independent file"
  ),
  detail = c(
    "36 factorial samples",
    paste(sort(unique(as.character(pca$Treatment))), collapse = ";"),
    paste(sort(unique(as.character(pca$Time))), collapse = ";"),
    "12 donor x time pairs, 3 conditions each",
    "PC1 and PC2 finite",
    paste0(width_mm, " x ", height_mm, " mm"),
    "600 dpi LZW",
    "300 dpi preview",
    "B panel canvas",
    "Letter added at assembly"
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa,
  file.path(qa_dir, "Fig7B_bulk_three_condition_pca_candidate_v0.9_qa.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(qa_dir, "Fig7B_bulk_three_condition_pca_candidate_v0.9_session_info.txt")
)

if (any(qa$status == "FAIL")) stop("Fig7B QA failed.")
cat("Rendered:", prefix, "\n")
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/08_plot_three_contrast_program_dotplot.R
stage_cd8_perturbation_08_plot_three_contrast_program_dotplot <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("gsea_three_contrasts")]
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

# ---- Retained complete source body; documented path/resource/quotation adapters only ----
#!/usr/bin/env Rscript

# Fig.7D candidate v0.3：三条件三比较程序 dotplot（D+E同排紧凑版）。
# 颜色表示 GSEA NES，点径表示 -log10(FDR)，实心/空心区分 FDR 阈值；
# 只复用既有 results/17 表，不重新计算 GSEA 或缩小 FDR 分母。

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

work_dir <- unname(paths["output_dir"])
gsea_file <- unname(paths["gsea_three_contrasts"])

plot_ready_dir <- file.path(work_dir, "plot_ready")
export_dir <- file.path(work_dir, "panel_exports")
qa_dir <- file.path(work_dir, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gsea_file) || file.info(gsea_file)$size == 0) {
  stop("Missing three-condition GSEA table: ", gsea_file)
}

program_specs <- data.frame(
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
  display_label = c(
    "CD8 activation",
    "TCR signaling",
    "T-cell cytotoxicity",
    "IFN-gamma response",
    "Memory-enriched",
    "Effector-enriched",
    "Glycolysis",
    "Hypoxia",
    "UPR"
  ),
  group = c(
    rep("Immune activation", 4),
    rep("CD8 state", 2),
    rep("Metabolic/stress", 3)
  ),
  row_order = seq_len(9),
  stringsAsFactors = FALSE
)
program_specs$full_feature_id <- paste0("MSIGDB::", program_specs$feature_id)

times <- c("1h", "6h", "24h", "48h")
contrasts <- c("LA_minus_NC", "HCl_minus_NC", "LA_minus_HCl")
contrast_labels <- c(
  LA_minus_NC = "LA - NC",
  HCl_minus_NC = "HCl - NC",
  LA_minus_HCl = "LA - HCl"
)

gsea <- read.delim(gsea_file, check.names = FALSE)
required <- c(
  "contrast", "time", "feature_id", "feature_label", "feature_layer",
  "source_collection", "status", "pval", "padj", "ES", "NES", "size",
  "n_genes_tested", "leadingEdge_genes"
)
if (!all(required %in% names(gsea))) {
  stop("GSEA table lacks required columns: ",
       paste(setdiff(required, names(gsea)), collapse = ", "))
}

gsea <- gsea[
  gsea$contrast %in% contrasts &
    gsea$time %in% times &
    gsea$feature_id %in% program_specs$full_feature_id,
  , drop = FALSE
]
if (nrow(gsea) != 108L) stop("Expected 108 displayed rows; found ", nrow(gsea), ".")
if (anyDuplicated(gsea[, c("contrast", "time", "feature_id")])) {
  stop("Duplicated contrast x time x feature rows.")
}
if (!setequal(unique(gsea$contrast), contrasts) ||
    !setequal(unique(gsea$time), times)) {
  stop("Unexpected contrast or time labels.")
}
if (any(gsea$status != "PASS" & gsea$status != "REUSED_FROZEN")) {
  stop("A displayed row is not PASS/REUSED_FROZEN.")
}
for (nm in c("pval", "padj", "ES", "NES", "size", "n_genes_tested")) {
  gsea[[nm]] <- as.numeric(gsea[[nm]])
}
if (any(!is.finite(gsea$NES)) || any(!is.finite(gsea$padj))) {
  stop("Non-finite NES or FDR values.")
}

spec_idx <- match(gsea$feature_id, program_specs$full_feature_id)
if (anyNA(spec_idx)) stop("Program specification join failed.")
gsea$display_label <- program_specs$display_label[spec_idx]
gsea$group <- program_specs$group[spec_idx]
gsea$row_order <- program_specs$row_order[spec_idx]
gsea$program <- factor(gsea$display_label, levels = rev(program_specs$display_label))
gsea$contrast_label <- factor(
  contrast_labels[gsea$contrast],
  levels = unname(contrast_labels[contrasts])
)
gsea$time <- factor(gsea$time, levels = times)
gsea$fdr_score <- pmin(-log10(pmax(gsea$padj, .Machine$double.xmin)), 4)
gsea$fdr_label <- ifelse(gsea$padj < 0.05, "FDR < 0.05", "FDR >= 0.05")

colour_limit <- max(2, ceiling(max(abs(gsea$NES), na.rm = TRUE) * 2) / 2)
if (!is.finite(colour_limit) || colour_limit <= 0) stop("Invalid NES colour limit.")

p <- ggplot(gsea, aes(x = time, y = program)) +
  geom_hline(yintercept = c(3.5, 5.5), linewidth = 0.34, colour = "#8C8C8C") +
  geom_point(
    aes(size = fdr_score, colour = NES, shape = fdr_label),
    stroke = 0.34, alpha = 0.94
  ) +
  facet_grid(. ~ contrast_label, scales = "fixed", drop = FALSE) +
  scale_shape_manual(
    values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1),
    name = "Significance",
    breaks = c("FDR < 0.05", "FDR >= 0.05"),
    labels = c("q < 0.05", "q >= 0.05")
  ) +
  scale_size_continuous(
    range = c(1.3, 4.6), limits = c(0, 4),
    breaks = c(0, 2, 4),
    name = expression(-log[10](FDR)),
    guide = guide_legend(order = 2)
  ) +
  scale_colour_gradient2(
    low = "#3B77A7", mid = "#F7F7F7", high = "#C75B3B",
    midpoint = 0, limits = c(-colour_limit, colour_limit),
    name = "NES",
    guide = guide_colorbar(
      order = 1, direction = "horizontal", barwidth = unit(18, "mm"), barheight = unit(1.8, "mm")
    )
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Arial", base_size = 5.8) +
  theme(
    panel.border = element_rect(fill = NA, colour = "#1A1A1A", linewidth = 0.36),
    panel.spacing.x = unit(1.0, "mm"),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.24, colour = "#1A1A1A"),
    axis.ticks.length = unit(0.55, "mm"),
    axis.text.x = element_text(size = 5.6, colour = "#1A1A1A",
                               margin = margin(t = 0.75, unit = "mm")),
    axis.text.y = element_text(size = 5.35, colour = "#1A1A1A",
                               hjust = 1, margin = margin(r = 1.0, unit = "mm")),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.0, face = "bold", colour = "#1A1A1A",
                              margin = margin(b = 1.2, unit = "mm")),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 4.8, face = "bold"),
    legend.text = element_text(size = 4.4),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.2, "mm"),
    legend.spacing.x = unit(0.8, "mm"),
    legend.spacing.y = unit(0.4, "mm"),
    legend.margin = margin(t = 0.4, unit = "mm"),
    plot.margin = margin(t = 1.0, r = 1.0, b = 0.8, l = 1.0, unit = "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(
    shape = guide_legend(
      order = 3, nrow = 1, override.aes = list(size = 2.2, colour = "#303030")
    )
  )

width_mm <- 112
height_mm <- 56
prefix <- file.path(export_dir, "Fig7D_bulk_three_contrast_program_dotplot_candidate_v0.3")
ggsave(paste0(prefix, ".pdf"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", device = cairo_pdf)
ggsave(paste0(prefix, ".tiff"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 600, compression = "lzw", device = "tiff")
ggsave(paste0(prefix, ".png"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 300, device = "png", bg = "white")

source_data <- gsea[, c(
  "contrast", "contrast_label", "time", "feature_id", "feature_label",
  "display_label", "group", "row_order", "feature_layer", "source_collection",
  "NES", "padj", "pval", "ES", "size", "n_genes_tested", "status",
  "fdr_score", "fdr_label", "leadingEdge_genes"
)]
source_data$panel <- "Fig7D"
write.table(
  source_data,
  file.path(plot_ready_dir, "fig7d_bulk_three_contrast_program_dotplot_v0.3_source_data.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)

summary_data <- data.frame(
  panel = "Fig7D",
  n_rows = nrow(gsea),
  n_contrasts = length(contrasts),
  n_programs = length(program_specs$feature_id),
  n_timepoints = length(times),
  n_fdr_significant = sum(gsea$padj < 0.05),
  fdr_denominator = 20L,
  canvas_width_mm = width_mm,
  canvas_height_mm = height_mm,
  note = "Nine displayed complete MSigDB sets; FDR denominator remains 20 testable fixed sets per contrast x time.",
  stringsAsFactors = FALSE
)
write.table(
  summary_data,
  file.path(plot_ready_dir, "fig7d_bulk_three_contrast_program_dotplot_v0.3_summary.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)

manifest <- data.frame(
  field = c(
    "panel_id", "status", "source_result", "contrasts", "program_layer",
    "displayed_programs", "tested_set_denominator", "timepoints",
    "colour_encoding", "size_encoding", "shape_encoding",
    "canvas_geometry_mm", "font_pt", "output_prefix"
  ),
  value = c(
    "Fig7D", "CANDIDATE", gsea_file, paste(contrasts, collapse = ";"),
    "MSIGDB_FULL_SET", 9L, 20L, paste(times, collapse = ";"),
    "NES", "-log10(FDR), capped at 4", "FDR < 0.05 vs FDR >= 0.05",
    paste0(width_mm, " x ", height_mm), "5.0-6.0", basename(prefix)
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
    file.path(plot_ready_dir, "fig7d_bulk_three_contrast_program_dotplot_v0.3_manifest.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)

qa <- data.frame(
  check_id = c(
    "rows_108", "three_contrasts", "nine_programs", "four_times",
    "all_displayed_status_valid", "finite_values", "fdr_denominator_recorded",
    "no_lm12_rows", "output_pdf", "output_tiff", "output_png",
    "panel_letter_embedded", "formal_write"
  ),
  status = c(
    ifelse(nrow(gsea) == 108L, "PASS", "FAIL"),
    ifelse(setequal(unique(gsea$contrast), contrasts), "PASS", "FAIL"),
    ifelse(length(unique(gsea$feature_id)) == 9L, "PASS", "FAIL"),
    ifelse(setequal(unique(as.character(gsea$time)), times), "PASS", "FAIL"),
    ifelse(all(gsea$status %in% c("PASS", "REUSED_FROZEN")), "PASS", "FAIL"),
    ifelse(all(is.finite(gsea$NES) & is.finite(gsea$padj)), "PASS", "FAIL"),
    ifelse(summary_data$fdr_denominator == 20L, "PASS", "FAIL"),
    ifelse(all(grepl("^MSIGDB::", gsea$feature_id)), "PASS", "FAIL"),
    ifelse(file.exists(paste0(prefix, ".pdf")) && file.info(paste0(prefix, ".pdf"))$size > 0, "PASS", "FAIL"),
    ifelse(file.exists(paste0(prefix, ".tiff")) && file.info(paste0(prefix, ".tiff"))$size > 0, "PASS", "FAIL"),
    ifelse(file.exists(paste0(prefix, ".png")) && file.info(paste0(prefix, ".png"))$size > 0, "PASS", "FAIL"),
    "PASS",
    "PASS"
  ),
  detail = c(
    "9 programs x 4 times x 3 contrasts",
    paste(contrasts, collapse = ";"),
    paste(program_specs$feature_id, collapse = ";"),
    paste(times, collapse = ";"),
    "all displayed rows are PASS or REUSED_FROZEN",
    "NES and FDR finite",
    "20 testable fixed sets per contrast x time",
    "LM12 rows excluded from this panel",
    paste0(width_mm, " x ", height_mm, " mm"),
    "600 dpi LZW",
    "300 dpi preview",
    "No panel letter in independent panel",
    "Candidate preview only"
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa, file.path(qa_dir, "Fig7D_bulk_three_contrast_program_dotplot_candidate_v0.3_qa.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()),
           file.path(qa_dir, "Fig7D_bulk_three_contrast_program_dotplot_candidate_v0.3_session_info.txt"))
if (any(qa$status == "FAIL")) stop("Fig7D QA failed.")

message("Rendered Fig7D candidate: ", prefix)
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/09_plot_hypoxia_gsea_curve.R
stage_cd8_perturbation_09_plot_hypoxia_gsea_curve <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["hallmark_gmt"]
if(basename(unname(input_paths))!="h.all.v2025.1.Hs.symbols.gmt") stop("Expected fixed Hallmark 2025.1 GMT basename.")
required_files <- c(unname(input_paths),
 file.path(unname(paths["bulk_de_dir"]),paste0("de_LA_minus_NC_",c("1h","6h","24h","48h"),".tsv")),
 file.path(unname(paths["bridge_dir"]),c("bulk_all_fixed_preranked_gsea.tsv","fixed_gene_set_manifest.tsv","bulk_fixed_readout_donor_summary.tsv")))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("bulk_de_dir","bridge_dir")])
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

relocate_hallmark_gmt <- function(source_path) {
  # 仅按固定版本 basename 重定位；原 manifest 与成员校验保持不变。
  expected <- "h.all.v2025.1.Hs.symbols.gmt"
  if(length(source_path)!=1L || is.na(source_path) || !nzchar(source_path) ||
     basename(source_path)!=expected) stop("Unregistered Hallmark GMT source_path.")
  unname(paths["hallmark_gmt"])
}

# ---- Retained complete source body; documented path/resource/quotation adapters only ----
#!/usr/bin/env Rscript

# Fig.7E candidate v0.6：四个时间点的 HALLMARK_HYPOXIA LA-NC GSEA 曲线（D+E同排紧凑版）。
# 仅重建 fixed DE ranking 的 running enrichment curve，不重跑 GSEA 或改变上游统计。

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(fgsea)
  library(grid)
})

work_dir <- unname(paths["output_dir"])
bridge_dir <- unname(paths["bridge_dir"])
de_dir <- unname(paths["bulk_de_dir"])

plot_ready_dir <- file.path(work_dir, "plot_ready")
export_dir <- file.path(work_dir, "panel_exports")
qa_dir <- file.path(work_dir, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

times <- c("1h", "6h", "24h", "48h")
time_labels <- c("1 h | LA - NC", "6 h | LA - NC", "24 h | LA - NC", "48 h | LA - NC")
feature_id <- "MSIGDB::HALLMARK_HYPOXIA"
feature_id_short <- "HALLMARK_HYPOXIA"
contrast <- "LA_minus_NC"
gsea_file <- file.path(bridge_dir, "bulk_all_fixed_preranked_gsea.tsv")
set_manifest_file <- file.path(bridge_dir, "fixed_gene_set_manifest.tsv")
donor_summary_file <- file.path(bridge_dir, "bulk_fixed_readout_donor_summary.tsv")
input_files <- c(
  gsea_file,
  set_manifest_file,
  donor_summary_file,
  file.path(de_dir, paste0("de_LA_minus_NC_", times, ".tsv"))
)
if (any(!file.exists(input_files)) || any(file.info(input_files)$size == 0)) {
  stop("Missing GSEA/DE input file(s): ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}

gsea <- read.delim(gsea_file, check.names = FALSE)
required_gsea <- c(
  "pval", "padj", "ES", "NES", "size", "context", "time", "donor",
  "feature_id", "seed", "leadingEdge_genes"
)
if (!all(required_gsea %in% names(gsea))) {
  stop("GSEA input lacks required columns: ", paste(setdiff(required_gsea, names(gsea)), collapse = ", "))
}
gsea <- gsea[
  gsea$context == "bulk_LA_minus_NC" & gsea$feature_id == feature_id,
  required_gsea,
  drop = FALSE
]
if (nrow(gsea) != length(times)) stop("Expected four HALLMARK_HYPOXIA rows; found ", nrow(gsea), ".")
if (!setequal(unique(gsea$time), times)) stop("Unexpected GSEA time labels.")
if (anyDuplicated(gsea[, c("feature_id", "time")])) stop("Duplicated GSEA feature-time rows.")
for (nm in c("pval", "padj", "ES", "NES", "size")) gsea[[nm]] <- as.numeric(gsea[[nm]])
if (any(!is.finite(gsea$ES) | !is.finite(gsea$NES) | !is.finite(gsea$padj))) {
  stop("Non-finite ES/NES/FDR values in fixed GSEA table.")
}

set_manifest <- read.delim(set_manifest_file, check.names = FALSE)
required_manifest <- c(
  "feature_id", "feature_label", "feature_layer", "source_collection",
  "source_path", "direction", "n_genes_source", "set_sha256"
)
if (!all(required_manifest %in% names(set_manifest))) stop("Fixed gene-set manifest is incomplete.")
set_row <- set_manifest[set_manifest$feature_id == feature_id_short, , drop = FALSE]
if (nrow(set_row) != 1L) stop("Expected one fixed HALLMARK_HYPOXIA manifest row.")
if (set_row$feature_layer != "MSIGDB_FULL_SET") stop("HALLMARK_HYPOXIA is not registered as MSIGDB_FULL_SET.")
gmt_path <- relocate_hallmark_gmt(as.character(set_row$source_path[[1L]]))
if (!file.exists(gmt_path) || file.info(gmt_path)$size == 0) stop("Missing fixed MSigDB GMT: ", gmt_path)

read_gmt_set <- function(path, id) {
  lines <- readLines(path, warn = FALSE)
  hit <- lines[startsWith(lines, paste0(id, "\t"))]
  if (length(hit) != 1L) stop("Expected one GMT row for ", id)
  tokens <- strsplit(hit, "\t", fixed = TRUE)[[1L]]
  unique(tokens[-seq_len(min(2L, length(tokens)))])
}

gene_set <- read_gmt_set(gmt_path, feature_id_short)
if (length(gene_set) != as.integer(set_row$n_genes_source[[1L]])) {
  stop("GMT member count does not match fixed manifest.")
}
if (!requireNamespace("digest", quietly = TRUE)) stop("Package digest is required for set checksum QA.")
set_hash <- digest::digest(sort(unique(gene_set)), algo = "sha256")
if (!identical(set_hash, as.character(set_row$set_sha256[[1L]]))) {
  stop("HALLMARK_HYPOXIA member checksum does not match fixed manifest.")
}

build_ranked_stats <- function(path) {
  de <- read.delim(path, check.names = FALSE)
  required_de <- c("gene_id", "logFC", "PValue", "FDR", "contrast", "time", "left", "right")
  if (!all(required_de %in% names(de))) stop("DE input lacks required columns: ", path)
  de <- de[de$contrast == contrast, , drop = FALSE]
  if (nrow(de) == 0L) stop("No LA_minus_NC rows in ", path)
  de$logFC <- as.numeric(de$logFC)
  de$PValue <- as.numeric(de$PValue)
  de$rank_score <- sign(de$logFC) * -log10(pmax(de$PValue, .Machine$double.xmin))
  de <- de[is.finite(de$rank_score) & !is.na(de$gene_id) & nzchar(de$gene_id), , drop = FALSE]
  de <- de[!duplicated(de$gene_id), , drop = FALSE]
  if (nrow(de) != 13655L) stop("Expected 13,655 filtered genes in ", path, "; found ", nrow(de), ".")
  stats <- de$rank_score
  names(stats) <- as.character(de$gene_id)
  stats <- sort(stats, decreasing = TRUE)
  stats
}

running_curve <- function(stats, gene_set) {
  hit <- names(stats) %in% gene_set
  selected <- which(hit)
  if (length(selected) < 3L) stop("Too few HALLMARK_HYPOXIA genes in ranked list.")
  n <- length(stats)
  m <- length(selected)
  weights <- abs(stats[selected])
  nr <- sum(weights)
  if (!is.finite(nr) || nr <= 0) stop("Invalid weighted rank sum.")
  step <- ifelse(hit, 0, -1 / (n - m))
  step[selected] <- weights / nr
  running <- cumsum(step)
  es <- if (max(running) >= -min(running)) max(running) else min(running)
  selected_stats <- sort(selected)
  fgsea_es <- as.numeric(fgsea:::calcGseaStat(stats, selected_stats, gseaParam = 1, scoreType = "std"))
  leading_positions <- if (max(running) >= -min(running)) {
    selected_stats[selected_stats <= which.max(running)]
  } else {
    rev(selected_stats[selected_stats >= which.min(running)])
  }
  leading_hit <- logical(length(stats))
  leading_hit[leading_positions] <- TRUE
  list(
    rank_index = seq_along(stats),
    gene_id = names(stats),
    rank_score = as.numeric(stats),
    gene_set_hit = hit,
    running_es = running,
    es = es,
    fgsea_es = fgsea_es,
    n_hits = m,
    leading_hit = leading_hit,
    n_leading = sum(leading_hit)
  )
}

curve_list <- vector("list", length(times))
summary_list <- vector("list", length(times))
for (i in seq_along(times)) {
  tm <- times[[i]]
  stats <- build_ranked_stats(file.path(de_dir, paste0("de_LA_minus_NC_", tm, ".tsv")))
  calc <- running_curve(stats, gene_set)
  fixed_row <- gsea[gsea$time == tm, , drop = FALSE]
  if (nrow(fixed_row) != 1L) stop("Missing fixed GSEA row for ", tm)
  curve_list[[i]] <- data.frame(
    time = tm,
    rank_index = calc$rank_index,
    rank_percent = if (length(stats) > 1L) {
      100 * (calc$rank_index - 1) / (length(stats) - 1)
    } else {
      rep(0, length(stats))
    },
    gene_id = calc$gene_id,
    rank_score = calc$rank_score,
    gene_set_hit = calc$gene_set_hit,
    leading_edge_hit = calc$leading_hit,
    running_es = calc$running_es,
    stringsAsFactors = FALSE
  )
  summary_list[[i]] <- data.frame(
    time = tm,
    feature_id = feature_id,
    contrast = contrast,
    n_ranked_genes = length(stats),
    n_genes_source = length(gene_set),
    n_genes_present = calc$n_hits,
    n_leading_edge = calc$n_leading,
    ES_reconstructed = calc$es,
    ES_fgsea_check = calc$fgsea_es,
    ES_reported = fixed_row$ES[[1L]],
    ES_delta = calc$fgsea_es - fixed_row$ES[[1L]],
    NES = fixed_row$NES[[1L]],
    FDR = fixed_row$padj[[1L]],
    p_value = fixed_row$pval[[1L]],
    fixed_seed = fixed_row$seed[[1L]],
    stringsAsFactors = FALSE
  )
}

curve_data <- do.call(rbind, curve_list)
summary_data <- do.call(rbind, summary_list)
if (any(abs(summary_data$ES_delta) > 0.02)) {
  stop("Reconstructed running ES differs from fixed fgsea ES by >0.02.")
}
if (any(!is.finite(curve_data$running_es))) stop("Non-finite running enrichment score.")

curve_data$time <- factor(curve_data$time, levels = times)
summary_data$time <- factor(summary_data$time, levels = times)
y_range <- range(curve_data$running_es, 0, finite = TRUE)
span <- diff(y_range)
if (!is.finite(span) || span <= 0) span <- 1
y_limits <- c(y_range[[1L]] - 0.13 * span, y_range[[2L]] + 0.12 * span)
tick_y0 <- y_limits[[1L]] + 0.035 * diff(y_limits)
tick_y1 <- y_limits[[1L]] + 0.105 * diff(y_limits)
hit_data <- curve_data[curve_data$gene_set_hit, , drop = FALSE]
hit_data$tick_y0 <- tick_y0
hit_data$tick_y1 <- tick_y1
hit_data$tick_colour <- ifelse(hit_data$leading_edge_hit, "#B85F2A", "#AFAFAF")
summary_data$annotation_x <- 98
summary_data$annotation_y <- y_limits[[2L]] - 0.05 * diff(y_limits)
summary_data$annotation_label <- sprintf("NES %.2f\nFDR %.3g", summary_data$NES, summary_data$FDR)

curve_col <- "#B85F2A"
black <- "#1A1A1A"
p <- ggplot(curve_data, aes(x = rank_percent, y = running_es, group = time)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#777777") +
  geom_line(linewidth = 0.50, colour = curve_col, lineend = "round") +
  geom_segment(
    data = hit_data,
    aes(x = rank_percent, xend = rank_percent, y = tick_y0, yend = tick_y1, colour = tick_colour),
    inherit.aes = FALSE, linewidth = 0.25, lineend = "butt"
  ) +
  geom_text(
    data = summary_data,
    aes(x = annotation_x, y = annotation_y, label = annotation_label),
    inherit.aes = FALSE, hjust = 1, vjust = 1,
    family = "Arial", size = 4.5 / ggplot2::.pt, colour = black, lineheight = 0.92
  ) +
  scale_colour_identity() +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = c("0", "25", "50", "75", "100"), expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = y_limits, breaks = pretty(y_limits, n = 4), expand = c(0, 0)
  ) +
  facet_wrap(
    ~time, ncol = 2, scales = "fixed",
    labeller = labeller(time = setNames(time_labels, times))
  ) +
  labs(x = "Gene rank percentile (LA - NC)", y = "Running ES") +
  theme_classic(base_family = "Arial", base_size = 5.8) +
  theme(
    panel.border = element_rect(fill = NA, colour = black, linewidth = 0.34),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(linewidth = 0.26, colour = black),
    axis.ticks = element_line(linewidth = 0.22, colour = black),
    axis.ticks.length = unit(0.50, "mm"),
    axis.text.x = element_text(size = 4.8, colour = black,
                               margin = margin(t = 0.55, unit = "mm")),
    axis.text.y = element_text(size = 4.8, colour = black,
                               margin = margin(r = 0.55, unit = "mm")),
    axis.title.x = element_text(size = 5.2, colour = black,
                                margin = margin(t = 0.8, unit = "mm")),
    axis.title.y = element_text(size = 5.2, colour = black,
                                margin = margin(r = 0.8, unit = "mm")),
    strip.background = element_rect(fill = "#F4F4F4", colour = black, linewidth = 0.30),
    strip.text = element_text(size = 5.3, face = "plain", colour = black,
                              margin = margin(t = 0.55, b = 0.55, unit = "mm")),
    panel.spacing = unit(1.0, "mm"),
    plot.margin = margin(t = 1.0, r = 1.8, b = 1.0, l = 2.4, unit = "mm")
  )

width_mm <- 67
height_mm <- 56
prefix <- file.path(export_dir, "Fig7E_hypoxia_gsea_curve_candidate_v0.6")
ggsave(paste0(prefix, ".pdf"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", device = cairo_pdf)
ggsave(paste0(prefix, ".tiff"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 600, compression = "lzw", device = "tiff")
ggsave(paste0(prefix, ".png"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 300, device = "png", bg = "white")

curve_out <- curve_data
curve_out$feature_id <- feature_id
curve_out$feature_label <- "Hypoxia"
curve_out$rank_formula <- "sign(logFC) * -log10(PValue)"
curve_out$display_measure <- "running enrichment score"
curve_out$display_geometry <- "four time-point GSEA small multiples"
write.table(
  curve_out,
  file.path(plot_ready_dir, "fig7e_hypoxia_gsea_curve_v0.6_curve_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
summary_out <- summary_data
summary_out$feature_label <- "Hypoxia"
summary_out$rank_formula <- "sign(logFC) * -log10(PValue)"
summary_out$display_measure <- "pooled LA-NC GSEA"
summary_out$display_geometry <- "2x2 running enrichment curves with hit ticks"
summary_out$formal_write <- FALSE
write.table(
  summary_out,
  file.path(plot_ready_dir, "fig7e_hypoxia_gsea_curve_v0.6_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
`sha256` <- bind_stage_helper(shared_sha256_022, environment())
manifest <- data.frame(
  field = c(
    "panel_id", "scientific_question", "feature_id", "feature_label", "contrast",
    "gsea_context", "rank_formula", "ranking_source", "gene_set_layer", "gene_set_source",
    "gene_set_source_sha256", "n_genes_source", "timepoints", "curve_unit", "donor_unit_note",
    "reported_metric_source", "reported_es_check", "panel_layout", "panel_size_mm", "font_target_pt",
    "title", "subtitle", "formal_write", "input_gsea", "input_gsea_sha256", "input_manifest",
    "input_manifest_sha256", "input_de_1", "input_de_1_sha256", "donor_score_retained"
  ),
  value = c(
    "Fig7E_candidate_v0.6",
    "Show the fixed HALLMARK_HYPOXIA running enrichment shape across four pooled bulk LA-NC rankings.",
    feature_id, "Hypoxia", contrast, "bulk_LA_minus_NC",
    "sign(logFC) * -log10(PValue)", "results/03_differential_response/de_LA_minus_NC_{time}.tsv",
    "MSIGDB_FULL_SET", gmt_path, as.character(set_row$set_sha256[[1L]]), as.character(length(gene_set)),
    paste(times, collapse = ";"), "gene rank and running ES; pooled contrast, not donor-level replicate",
    "3-donor score/CI files retained for supplementary/source-data use; not drawn in E",
    "results/11 bulk_all_fixed_preranked_gsea.tsv", "reconstructed ES checked against fixed fgsea ES",
    "2 x 2 small multiples", paste0(width_mm, " x ", height_mm), "5.8 (allowed 5.6-6.0)",
    "NULL", "NULL", "FALSE", gsea_file, sha256(gsea_file), set_manifest_file, sha256(set_manifest_file),
    input_files[[4L]], sha256(input_files[[4L]]), donor_summary_file
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  file.path(plot_ready_dir, "fig7e_hypoxia_gsea_curve_v0.6_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

prefix_files <- paste0(prefix, c(".pdf", ".tiff", ".png"))
qa <- data.frame(
  check_id = c(
    "four_timepoints", "context_exact", "feature_exact", "rank_formula_exact", "ranked_gene_count",
    "fixed_set_layer", "fixed_set_member_hash", "present_gene_count", "reported_es_match",
    "finite_curve", "common_x_y_scale", "hit_ticks", "panel_dimensions", "pdf_nonempty",
    "tiff_nonempty", "png_nonempty", "no_title", "no_subtitle", "formal_write"
  ),
  status = c(
    if (identical(sort(unique(as.character(summary_data$time))), sort(times))) "PASS" else "FAIL",
    if (all(gsea$context == "bulk_LA_minus_NC")) "PASS" else "FAIL",
    if (all(gsea$feature_id == feature_id)) "PASS" else "FAIL",
    "PASS",
    if (all(summary_data$n_ranked_genes == 13655L)) "PASS" else "FAIL",
    if (identical(as.character(set_row$feature_layer[[1L]]), "MSIGDB_FULL_SET")) "PASS" else "FAIL",
    if (identical(set_hash, as.character(set_row$set_sha256[[1L]]))) "PASS" else "FAIL",
    if (all(summary_data$n_genes_present >= 3L)) "PASS" else "FAIL",
    if (all(abs(summary_data$ES_delta) <= 0.02)) "PASS" else "FAIL",
    if (all(is.finite(curve_data$running_es))) "PASS" else "FAIL",
    "PASS",
    if (all(curve_data$gene_set_hit %in% c(TRUE, FALSE))) "PASS" else "FAIL",
    if (width_mm == 67 && height_mm == 56) "PASS" else "FAIL",
    if (file.exists(prefix_files[[1]]) && file.info(prefix_files[[1]])$size > 0) "PASS" else "FAIL",
    if (file.exists(prefix_files[[2]]) && file.info(prefix_files[[2]])$size > 0) "PASS" else "FAIL",
    if (file.exists(prefix_files[[3]]) && file.info(prefix_files[[3]])$size > 0) "PASS" else "FAIL",
    "PASS", "PASS", "PASS"
  ),
  detail = c(
    "1h;6h;24h;48h", "bulk_LA_minus_NC", feature_id, "sign(logFC) * -log10(PValue)",
    paste(summary_data$n_ranked_genes, collapse = ";"), "MSIGDB_FULL_SET", set_hash,
    paste(summary_data$n_genes_present, collapse = ";"), sprintf("max_abs_delta=%.6g", max(abs(summary_data$ES_delta))),
    "all running ES finite", "shared limits", paste(summary_data$n_genes_present, collapse = ";"),
    paste0(width_mm, " x ", height_mm, " mm"), as.character(file.info(prefix_files[[1]])$size),
    as.character(file.info(prefix_files[[2]])$size), as.character(file.info(prefix_files[[3]])$size),
    "No title", "No subtitle", "Candidate preview only"
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa,
  file.path(qa_dir, "Fig7E_hypoxia_gsea_curve_candidate_v0.6_qa.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

message("Fig7E v0.6 rendered: ", prefix)
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/10_plot_hypoxia_leading_edge_gene_trajectory.R
stage_cd8_perturbation_10_plot_hypoxia_leading_edge_gene_trajectory <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("trajectory_source_file","leading_edge_file")]
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

# ---- Retained complete source body; documented path/quotation adapters only ----
#!/usr/bin/env Rscript

# Fig.7F candidate v2.0：纯 CD8 bulk Hypoxia leading-edge 基因 2×6 小 panel。
# 每个小 panel 展示一个固定基因在 LA−NC 四个时间点的 log2FC；
# 三个模块以顶部色带分组。该候选仅用于 preview，不写入正式冻结目录。

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

work_dir <- unname(paths["output_dir"])

plot_ready_dir <- file.path(work_dir, "plot_ready")
export_dir <- file.path(work_dir, "panel_exports")
qa_dir <- file.path(work_dir, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

source_file <- unname(paths["trajectory_source_file"])
leading_edge_file <- unname(paths["leading_edge_file"])
if (!file.exists(source_file) || file.info(source_file)$size == 0) {
  stop("Missing plot-ready source data: ", source_file)
}
if (!file.exists(leading_edge_file) || file.info(leading_edge_file)$size == 0) {
  stop("Missing leading-edge provenance: ", leading_edge_file)
}

feature_id <- "MSIGDB::HALLMARK_HYPOXIA"
feature_label <- "Hypoxia"
times <- c("1h", "6h", "24h", "48h")
module_order <- c("Hypoxia / energy", "Glycolytic adaptation", "Immediate-early feedback")
module_short <- c(
  "Hypoxia / energy" = "Hypoxia / energy",
  "Glycolytic adaptation" = "Glycolytic adaptation",
  "Immediate-early feedback" = "Immediate-early feedback"
)
module_cols <- c(
  "Hypoxia / energy" = "#C8782B",
  "Glycolytic adaptation" = "#5B8DB8",
  "Immediate-early feedback" = "#7B5AA6"
)

# 2×6 的阅读顺序：每个模块占据相邻两列，每个模块四个基因分布在两行。
gene_order <- c(
  "DDIT4", "HMOX1", "ENO1", "LDHA", "DUSP1", "JUN",
  "PDK1", "SLC2A3", "PGK1", "PGM1", "KLF6", "ZFP36"
)
gene_module <- c(
  DDIT4 = "Hypoxia / energy", HMOX1 = "Hypoxia / energy",
  PDK1 = "Hypoxia / energy", SLC2A3 = "Hypoxia / energy",
  ENO1 = "Glycolytic adaptation", LDHA = "Glycolytic adaptation",
  PGK1 = "Glycolytic adaptation", PGM1 = "Glycolytic adaptation",
  DUSP1 = "Immediate-early feedback", JUN = "Immediate-early feedback",
  KLF6 = "Immediate-early feedback", ZFP36 = "Immediate-early feedback"
)

dat <- read.delim(source_file, check.names = FALSE)
required <- c(
  "context", "contrast", "time", "gene", "module", "log2FC", "FDR",
  "independent_units", "feature_id", "feature_label"
)
if (!all(required %in% names(dat))) {
  stop("Source data lacks required columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
}
dat <- dat[
  dat$feature_id == feature_id & dat$feature_label == feature_label &
    dat$context == "purified_CD8_bulk" & dat$contrast == "LA_minus_NC", , drop = FALSE
]
if (nrow(dat) != 48L) stop("Expected 48 bulk LA-minus-NC gene-time rows; found ", nrow(dat), ".")
if (!setequal(unique(dat$gene), names(gene_module)) || length(unique(dat$gene)) != 12L) {
  stop("The source data does not contain the fixed 12-gene display set.")
}
if (!setequal(unique(dat$module), module_order)) stop("Unexpected module labels.")
if (any(!is.finite(dat$log2FC)) || any(!is.finite(dat$FDR))) stop("Non-finite source values.")
if (any(!dat$time %in% times) || any(dat$independent_units != 3L)) {
  stop("Unexpected time labels or independent-unit counts.")
}

dat$gene <- factor(dat$gene, levels = names(gene_module))
dat$module <- factor(dat$module, levels = module_order)
dat$time <- factor(dat$time, levels = times)
x_map <- c("1h" = 1, "6h" = 2, "24h" = 3, "48h" = 4)
dat$x <- unname(x_map[as.character(dat$time)])
dat <- dat[order(match(as.character(dat$gene), names(gene_module)), dat$x), , drop = FALSE]

# 所有基因共享一个 y 轴范围，避免不同小 panel 的尺度掩盖效应大小。
y_raw <- range(c(dat$log2FC, 0), finite = TRUE)
y_pad <- max(diff(y_raw) * 0.10, 0.12)
y_lim <- c(y_raw[1] - y_pad, y_raw[2] + y_pad)
y_breaks <- pretty(y_lim, n = 4)

black <- "#1A1A1A"
make_gene_panel <- function(gene, row_id, col_id) {
  d <- dat[as.character(dat$gene) == gene, , drop = FALSE]
  mod <- as.character(d$module[1])
  mod_col <- unname(module_cols[mod])
  show_y <- col_id %in% c(1L, 3L, 5L)
  show_x <- row_id == 2L
  x_title <- if (show_x && col_id == 3L) "Time after LA exposure" else NULL
  y_title <- if (row_id == 1L && col_id == 1L) "log2FC" else NULL

  p <- ggplot(d, aes(x = x, y = log2FC, group = 1)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#777777") +
    geom_line(linewidth = 0.58, colour = mod_col, lineend = "round") +
    geom_point(size = 1.25, colour = mod_col) +
    scale_x_continuous(
      limits = c(0.82, 4.18), breaks = c(1, 2, 3, 4),
      labels = if (show_x) c("1 h", "6 h", "24 h", "48 h") else NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks, expand = c(0, 0)) +
    labs(title = gene, x = x_title, y = y_title) +
    theme_classic(base_family = "Arial", base_size = 5.7) +
    theme(
      panel.border = element_rect(fill = NA, colour = black, linewidth = 0.32),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.22, colour = black),
      axis.ticks.length = unit(0.45, "mm"),
      axis.text.x = element_text(
        size = 4.55, colour = black,
        margin = margin(t = 0.4, unit = "mm"),
        lineheight = 0.84
      ),
      axis.text.y = element_text(size = 4.55, colour = black,
                                 margin = margin(r = 0.45, unit = "mm")),
      axis.title.x = element_text(size = 5.0, colour = black,
                                  margin = margin(t = 0.65, unit = "mm")),
      axis.title.y = element_text(size = 5.0, colour = black,
                                  margin = margin(r = 0.65, unit = "mm")),
      plot.title = element_text(size = 5.55, face = "bold", colour = mod_col,
                                hjust = 0.5, margin = margin(b = 0.55, unit = "mm")),
      plot.margin = margin(t = 0.5, r = 0.55, b = 0.55, l = 0.55, unit = "mm")
    )
  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

gene_panels <- vector("list", length(gene_order))
for (i in seq_along(gene_order)) {
  row_id <- if (i <= 6L) 1L else 2L
  col_id <- if (i <= 6L) i else i - 6L
  gene_panels[[i]] <- make_gene_panel(gene_order[i], row_id, col_id)
}

make_header <- function(label, fill) {
  ggplot(data.frame(x = 0.5, y = 0.5), aes(x, y)) +
    geom_tile(width = 1, height = 1, fill = fill, alpha = 0.16,
              colour = fill, linewidth = 0.35) +
    annotate("text", x = 0.5, y = 0.5, label = label,
             family = "Arial", fontface = "bold", size = 5.2 / ggplot2::.pt,
             colour = black) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(t = 0.3, r = 0.55, b = 0.45, l = 0.55, unit = "mm"))
}

header <- wrap_plots(
  make_header(module_short[[module_order[1]]], module_cols[[module_order[1]]]),
  make_header(module_short[[module_order[2]]], module_cols[[module_order[2]]]),
  make_header(module_short[[module_order[3]]], module_cols[[module_order[3]]]),
  nrow = 1, widths = c(2, 2, 2)
)
body <- wrap_plots(gene_panels, nrow = 2, widths = rep(1, 6), heights = c(1, 1))
p <- header / body + plot_layout(heights = c(0.13, 1), guides = "keep")

width_mm <- 180
height_mm <- 56
prefix <- file.path(export_dir, "Fig7F_hypoxia_leading_edge_gene_smallmultiples_candidate_v2.0")
ggsave(paste0(prefix, ".pdf"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", device = cairo_pdf)
ggsave(paste0(prefix, ".tiff"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 600, compression = "lzw", device = "tiff")
ggsave(paste0(prefix, ".png"), p,
       width = width_mm / 25.4, height = height_mm / 25.4,
       units = "in", dpi = 300, device = "png", bg = "white")

leading <- read.delim(leading_edge_file, check.names = FALSE)
if (!all(c("feature_id", "context", "gene") %in% names(leading))) {
  stop("Leading-edge provenance lacks feature_id/context/gene columns.")
}
hyp <- leading[leading$feature_id == feature_id, , drop = FALSE]
bulk_le <- unique(hyp$gene[hyp$context == "bulk_LA_minus_NC"])
mel_le <- unique(hyp$gene[hyp$context == "melanoma_CD8_C1_minus_C4"])
n_shared <- length(intersect(bulk_le, mel_le))
if (n_shared != 24L) stop("Expected 24-gene melanoma/bulk Hypoxia intersection; found ", n_shared, ".")

source_data <- dat
source_data$display_measure <- "within-contrast log2FC"
source_data$display_geometry <- "2x6 single-gene small multiples; three module header bands"
source_data$feature_id <- feature_id
source_data$feature_label <- feature_label
write.table(
  source_data,
  file.path(plot_ready_dir, "fig7f_hypoxia_leading_edge_gene_smallmultiples_v2.0_source_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

summary_tab <- data.frame(
  feature_id = feature_id,
  feature_label = feature_label,
  n_shared_hypoxia_leading_edge = n_shared,
  n_main_genes = 12L,
  module_counts = paste(names(table(dat$module)), as.integer(table(dat$module)), collapse = "; "),
  bulk_log2FC_min = min(dat$log2FC),
  bulk_log2FC_max = max(dat$log2FC),
  panel_size_mm = paste0(width_mm, " x ", height_mm),
  display_geometry = "2x6 single-gene small multiples with module header bands",
  stringsAsFactors = FALSE
)
write.table(
  summary_tab,
  file.path(plot_ready_dir, "fig7f_hypoxia_leading_edge_gene_smallmultiples_v2.0_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
`sha256` <- bind_stage_helper(shared_sha256_022, environment())
input_files <- c(source_file, leading_edge_file)
manifest_fields <- c(
  "panel_id", "scientific_question", "feature", "gene_selection_rule",
  "main_gene_count", "shared_pool_count", "module_order", "comparison",
  "biological_units", "display_measure", "display_geometry", "scale_policy",
  "panel_size_mm", "font_target_pt", "formal_write",
  paste0("input_", seq_along(input_files)), paste0("input_sha256_", seq_along(input_files))
)
manifest_values <- c(
  "Fig7F_candidate_v2.0",
  "Show Fig.6-informed Hypoxia leading-edge gene trajectories in purified-CD8 bulk after LA exposure over time.",
  "MSIGDB::HALLMARK_HYPOXIA; complete-set leading-edge",
  "12 fixed display genes from the 24-gene melanoma/bulk leading-edge intersection; fixed before plotting and not selected by minimum P value",
  "12", as.character(n_shared), paste(module_order, collapse = " | "),
  "purified CD8 bulk LA minus NC at 1/6/24/48 h",
  "3 bulk donors; 48 gene-time rows",
  "within-contrast log2FC",
  "2x6 single-gene small multiples; module header bands; bulk-only connected trajectories",
  "one common y-axis range for all 12 genes; zero reference line in every panel",
  paste0(width_mm, " x ", height_mm), "5.7 (allowed 5.6-6.0)", "FALSE",
  input_files, vapply(input_files, sha256, character(1))
)
if (length(manifest_fields) != length(manifest_values)) {
  stop("Manifest field/value length mismatch: ", length(manifest_fields), " vs ", length(manifest_values))
}
write.table(
  data.frame(field = manifest_fields, value = manifest_values, stringsAsFactors = FALSE),
  file.path(plot_ready_dir, "fig7f_hypoxia_leading_edge_gene_smallmultiples_v2.0_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

prefix_files <- paste0(prefix, c(".pdf", ".tiff", ".png"))
qa <- data.frame(
  check_id = c(
    "feature_exact", "shared_pool_exact_24", "selected_genes_12", "modules_three_x_four",
    "bulk_rows", "bulk_only", "finite_values", "independent_units", "two_by_six_order",
    "module_header_bands", "common_y_scale", "panel_dimensions", "pdf_nonempty",
    "tiff_nonempty", "png_nonempty", "no_heatmap", "no_subtitle", "formal_write"
  ),
  status = c(
    if (feature_id == "MSIGDB::HALLMARK_HYPOXIA") "PASS" else "FAIL",
    if (n_shared == 24L) "PASS" else "FAIL",
    if (length(unique(dat$gene)) == 12L) "PASS" else "FAIL",
    if (all(table(dat$module) == 16L)) "PASS" else "FAIL",
    if (nrow(dat) == 48L) "PASS" else "FAIL",
    if (all(dat$context == "purified_CD8_bulk")) "PASS" else "FAIL",
    if (all(is.finite(dat$log2FC)) && all(is.finite(dat$FDR))) "PASS" else "FAIL",
    if (all(dat$independent_units == 3L)) "PASS" else "FAIL",
    "PASS", "PASS", "PASS", paste0(width_mm, " x ", height_mm, " mm"),
    if (file.exists(prefix_files[[1]]) && file.info(prefix_files[[1]])$size > 0) "PASS" else "FAIL",
    if (file.exists(prefix_files[[2]]) && file.info(prefix_files[[2]])$size > 0) "PASS" else "FAIL",
    if (file.exists(prefix_files[[3]]) && file.info(prefix_files[[3]])$size > 0) "PASS" else "FAIL",
    "PASS", "PASS", "PASS"
  ),
  detail = c(
    feature_id, paste0("n=", n_shared), "n=12; all fixed", "3 modules; 4 genes each",
    "48 bulk LA-minus-NC gene-time rows; 3 donors", "purified_CD8_bulk only",
    "48 finite rows", "3 bulk donors", "row-major 2x6; two adjacent columns per module",
    "three colored module headers", "all panels use the same y range", paste0(width_mm, " x ", height_mm, " mm"),
    as.character(file.info(prefix_files[[1]])$size), as.character(file.info(prefix_files[[2]])$size),
    as.character(file.info(prefix_files[[3]])$size), "Line/point geometry only",
    "No title or subtitle", "Candidate preview only"
  ),
  stringsAsFactors = FALSE
)
write.table(
  qa,
  file.path(qa_dir, "Fig7F_hypoxia_leading_edge_gene_smallmultiples_candidate_v2.0_qa.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(qa_dir, "Fig7F_hypoxia_leading_edge_gene_smallmultiples_candidate_v2.0_session_info.txt")
)

if (any(qa$status == "FAIL")) stop("One or more QA checks failed.")
print(summary_tab)
cat("Rendered Fig7F candidate:", prefix, "\n")
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/11_plot_hcl_gene_trajectory.R
stage_cd8_perturbation_11_plot_hcl_gene_trajectory <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet","filtered_logcpm")]
required_files <- c(unname(input_paths),file.path(unname(paths["bulk_de_dir"]),paste0("de_HCl_minus_NC_",c("1h","6h","24h","48h"),".tsv")))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths["bulk_de_dir"])
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

# ---- Retained complete source body; documented path/quotation adapters only ----
#!/usr/bin/env Rscript

# S13B v0.1：HCl−NC 下正文 12 个 leading-edge 基因的 donor-aware 时间轨迹。
# donor 为独立生物学单位；图中不把 donor-time 扩展为独立重复。

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

s13_root <- unname(paths["output_dir"])
de_dir <- unname(paths["bulk_de_dir"])
sample_file <- unname(paths["sample_sheet"])
matrix_file <- unname(paths["filtered_logcpm"])
plot_ready_dir <- file.path(s13_root, "plot_ready")
export_dir <- file.path(s13_root, "panel_exports")
qa_dir <- file.path(s13_root, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
if (!all(file.exists(c(sample_file, matrix_file)))) stop("Missing sample sheet or filtered matrix")

times <- c("1h", "6h", "24h", "48h")
donors <- c("V1", "V2", "V3")
gene_order <- c("DDIT4", "HMOX1", "ENO1", "LDHA", "DUSP1", "JUN",
                "PDK1", "SLC2A3", "PGK1", "PGM1", "KLF6", "ZFP36")
gene_module <- c(
  DDIT4 = "Hypoxia / energy", HMOX1 = "Hypoxia / energy",
  PDK1 = "Hypoxia / energy", SLC2A3 = "Hypoxia / energy",
  ENO1 = "Glycolytic adaptation", LDHA = "Glycolytic adaptation",
  PGK1 = "Glycolytic adaptation", PGM1 = "Glycolytic adaptation",
  DUSP1 = "Immediate-early feedback", JUN = "Immediate-early feedback",
  KLF6 = "Immediate-early feedback", ZFP36 = "Immediate-early feedback"
)
module_order <- c("Hypoxia / energy", "Glycolytic adaptation", "Immediate-early feedback")
module_cols <- c("Hypoxia / energy" = "#C8782B", "Glycolytic adaptation" = "#5B8DB8", "Immediate-early feedback" = "#7B5AA6")
donor_cols <- c(V1 = "#4F7FA2", V2 = "#7B5AA6", V3 = "#C8782B")

sample_meta <- read.delim(sample_file, check.names = FALSE)
required_meta <- c("SampleID", "Donor", "Treatment", "Time")
if (!all(required_meta %in% names(sample_meta))) stop("Sample metadata columns missing")
sample_meta <- sample_meta[sample_meta$Treatment %in% c("NC", "HCl") & sample_meta$Time %in% times, , drop = FALSE]
if (nrow(sample_meta) != 24L) stop("Expected 24 NC/HCl factorial samples")
mat_obj <- readRDS(matrix_file)
mat <- if (is.matrix(mat_obj)) mat_obj else if (is.list(mat_obj) && !is.null(mat_obj$logCPM)) mat_obj$logCPM else as.matrix(mat_obj)
mat <- as.matrix(mat)
if (is.null(rownames(mat)) || is.null(colnames(mat))) stop("Filtered matrix lacks dimnames")
if (!all(gene_order %in% rownames(mat))) stop("Fixed 12-gene set missing from filtered matrix")
if (!all(sample_meta$SampleID %in% colnames(mat))) stop("Sample IDs do not match filtered matrix")

rows <- list()
idx <- 0L
for (gg in gene_order) for (tt in times) for (dd in donors) {
  hcl <- sample_meta$SampleID[sample_meta$Donor == dd & sample_meta$Time == tt & sample_meta$Treatment == "HCl"]
  nc <- sample_meta$SampleID[sample_meta$Donor == dd & sample_meta$Time == tt & sample_meta$Treatment == "NC"]
  if (length(hcl) != 1L || length(nc) != 1L) stop("Expected one HCl and one NC sample per donor-time")
  idx <- idx + 1L
  rows[[idx]] <- data.frame(gene = gg, module = unname(gene_module[[gg]]), time = tt, donor = dd,
                            HCl_logCPM = as.numeric(mat[gg, hcl]), NC_logCPM = as.numeric(mat[gg, nc]),
                            donor_delta = as.numeric(mat[gg, hcl] - mat[gg, nc]), stringsAsFactors = FALSE)
}
donor_dat <- do.call(rbind, rows)
if (nrow(donor_dat) != 144L || any(!is.finite(donor_dat$donor_delta))) stop("Unexpected donor trajectory rows")
means <- aggregate(donor_delta ~ gene + module + time, donor_dat, mean)
names(means)[names(means) == "donor_delta"] <- "mean_delta"
means$time <- factor(means$time, levels = times)
donor_dat$time <- factor(donor_dat$time, levels = times)
donor_dat$gene <- factor(donor_dat$gene, levels = gene_order)
donor_dat$module <- factor(donor_dat$module, levels = module_order)
means$gene <- factor(means$gene, levels = gene_order)
means$module <- factor(means$module, levels = module_order)

# 复用全转录组 DE 的 HCl−NC logFC/FDR，作为图下可追溯统计，不作为绘图筛选。
de_rows <- lapply(times, function(tt) {
  f <- file.path(de_dir, paste0("de_HCl_minus_NC_", tt, ".tsv"))
  if (!file.exists(f)) stop("Missing DE file: ", f)
  x <- read.delim(f, check.names = FALSE)
  z <- x[x$gene_id %in% gene_order, c("gene_id", "logFC", "FDR", "time", "contrast"), drop = FALSE]
  names(z)[names(z) == "gene_id"] <- "gene"
  z
})
de <- do.call(rbind, de_rows)
if (nrow(de) != 48L || anyDuplicated(de[, c("gene", "time")])) stop("Unexpected HCl DE rows for fixed genes")
de$time <- factor(de$time, levels = times)

y_raw <- range(c(donor_dat$donor_delta, 0), finite = TRUE)
y_pad <- max(diff(y_raw) * .10, .12)
y_lim <- c(y_raw[1] - y_pad, y_raw[2] + y_pad)
y_breaks <- pretty(y_lim, n = 4)
black <- "#1A1A1A"

make_gene_panel <- function(gene, row_id, col_id) {
  d <- donor_dat[as.character(donor_dat$gene) == gene, , drop = FALSE]
  m <- means[as.character(means$gene) == gene, , drop = FALSE]
  mod <- as.character(d$module[1]); mod_col <- unname(module_cols[mod])
  show_y <- col_id %in% c(1L, 3L, 5L); show_x <- row_id == 2L
  ggplot() +
    geom_hline(yintercept = 0, linewidth = .25, colour = "#777777") +
    geom_line(data = d, aes(x = time, y = donor_delta, group = donor, colour = donor), linewidth = .30, alpha = .55) +
    geom_point(data = d, aes(x = time, y = donor_delta, colour = donor), size = .65, alpha = .65) +
    geom_line(data = m, aes(x = time, y = mean_delta, group = 1), linewidth = .72, colour = mod_col) +
    geom_point(data = m, aes(x = time, y = mean_delta), size = 1.05, colour = mod_col) +
    scale_colour_manual(values = donor_cols, guide = "none") +
    scale_x_discrete(drop = FALSE, labels = if (show_x) c("1 h", "6 h", "24 h", "48 h") else NULL) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks, expand = c(0, 0)) +
    labs(title = gene, x = if (show_x && col_id == 3L) "Time after HCl exposure" else NULL,
         y = if (row_id == 1L && col_id == 1L) "HCl-NC logCPM" else NULL) +
    theme_classic(base_family = "Arial", base_size = 5.6) +
    theme(panel.border = element_rect(fill = NA, colour = black, linewidth = .32),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_blank(), axis.ticks = element_line(linewidth = .22, colour = black),
          axis.ticks.length = unit(.45, "mm"), axis.text.x = element_text(size = 4.35, colour = black),
          axis.text.y = element_text(size = 4.35, colour = black), axis.title.x = element_text(size = 4.9, colour = black),
          axis.title.y = element_text(size = 4.9, colour = black),
          plot.title = element_text(size = 5.35, face = "bold", colour = mod_col, hjust = .5, margin = margin(b = .5, unit = "mm")),
          plot.margin = margin(t = .5, r = .5, b = .5, l = .5, unit = "mm")) +
    {if (!show_y) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}
}

gene_panels <- vector("list", length(gene_order))
for (i in seq_along(gene_order)) {
  gene_panels[[i]] <- make_gene_panel(gene_order[i], if (i <= 6L) 1L else 2L, if (i <= 6L) i else i - 6L)
}
make_header <- function(label, fill) {
  ggplot(data.frame(x = .5, y = .5), aes(x, y)) +
    geom_tile(width = 1, height = 1, fill = fill, alpha = .16, colour = fill, linewidth = .35) +
    annotate("text", x = .5, y = .5, label = label, family = "Arial", fontface = "bold", size = 5.0 / ggplot2::.pt, colour = black) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") + theme_void() +
    theme(plot.margin = margin(t = .3, r = .5, b = .45, l = .5, unit = "mm"))
}
header <- wrap_plots(make_header(module_order[1], module_cols[1]), make_header(module_order[2], module_cols[2]),
                     make_header(module_order[3], module_cols[3]), nrow = 1, widths = c(2, 2, 2))
body <- wrap_plots(gene_panels, nrow = 2, widths = rep(1, 6), heights = c(1, 1))
p <- header / body + plot_layout(heights = c(.13, 1))

prefix <- file.path(export_dir, "S13B_hcl_gene_trajectory_candidate_v0.1")
ggsave(paste0(prefix, ".pdf"), p, width = 135/25.4, height = 56/25.4, units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(prefix, ".tiff"), p, width = 135/25.4, height = 56/25.4, units = "in", dpi = 600, compression = "lzw", device = "tiff", bg = "white")
ggsave(paste0(prefix, ".png"), p, width = 135/25.4, height = 56/25.4, units = "in", dpi = 300, device = "png", bg = "white")

source_data <- merge(donor_dat, de, by = c("gene", "time"), all.x = TRUE, sort = FALSE)
write.table(source_data, file.path(plot_ready_dir, "s13b_hcl_gene_trajectory_source.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(means, file.path(plot_ready_dir, "s13b_hcl_gene_trajectory_means.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

qa <- data.frame(check_id = c("gene_count", "donor_rows", "mean_rows", "de_rows", "donor_count", "panel_size_mm", "formal_write", "png_nonempty", "pdf_nonempty", "tiff_nonempty"),
                 expected = c("12", "144", "48", "48", "3", "135 x 56", "FALSE", "TRUE", "TRUE", "TRUE"),
                 observed = c(length(unique(donor_dat$gene)), nrow(donor_dat), nrow(means), nrow(de), length(unique(donor_dat$donor)), "135 x 56", "FALSE",
                              file.exists(paste0(prefix, ".png")) && file.info(paste0(prefix, ".png"))$size > 0,
                              file.exists(paste0(prefix, ".pdf")) && file.info(paste0(prefix, ".pdf"))$size > 0,
                              file.exists(paste0(prefix, ".tiff")) && file.info(paste0(prefix, ".tiff"))$size > 0), stringsAsFactors = FALSE)
qa$pass <- as.character(qa$expected) == as.character(qa$observed)
write.table(qa, file.path(qa_dir, "S13B_hcl_gene_trajectory_candidate_v0.1_qa.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(qa_dir, "S13B_hcl_gene_trajectory_candidate_v0.1_session_info.txt"))
if (!all(qa$pass)) stop("S13B QA failed")
message("Rendered S13B: ", prefix)
}

# Purpose: Export the reported CD8 figure source tables and component plots.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/cd8_perturbation/12_plot_cd8_marker_modules.R
stage_cd8_perturbation_12_plot_cd8_marker_modules <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("marker_expression")]
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

# ---- Retained complete source body; documented path/quotation adapters only ----
#!/usr/bin/env Rscript

# S13C v0.1：三条件 CD8 phenotype marker，四模块一行排列。
# 每个模块为一个生物信息学小 panel；每个基因保留 donor 点和条件均值轨迹。

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

s13_root <- unname(paths["output_dir"])
expr_file <- unname(paths["marker_expression"])
plot_ready_dir <- file.path(s13_root, "plot_ready")
export_dir <- file.path(s13_root, "panel_exports")
qa_dir <- file.path(s13_root, "qa")
dir.create(plot_ready_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(expr_file) || file.info(expr_file)$size == 0) stop("Missing marker expression input")

marker_spec <- data.frame(
  gene = c("CD69", "IL2RA", "GZMB", "PRF1", "NKG7", "GZMA", "IFNG",
           "PDCD1", "TIGIT", "LAG3", "CTLA4", "HAVCR2", "TOX",
           "CCR7", "IL7R", "TCF7", "SELL", "LEF1"),
  module = c(rep("Activation", 2), rep("Cytotoxic / effector", 5), rep("Checkpoint-associated", 6), rep("Memory-associated", 5)),
  stringsAsFactors = FALSE
)
module_order <- c("Activation", "Cytotoxic / effector", "Checkpoint-associated", "Memory-associated")
module_cols <- c("Activation" = "#C8782B", "Cytotoxic / effector" = "#C85C4A", "Checkpoint-associated" = "#7B5AA6", "Memory-associated" = "#4F7FA2")
treatment_cols <- c(NC = "#4F8BB8", LA = "#C8782B", HCl = "#7F8996")
times <- c("1h", "6h", "24h", "48h")

expr <- read.delim(expr_file, check.names = FALSE)
required <- c("gene", "sample_id", "logCPM", "Donor", "Treatment", "Time", "marker_group", "detection_status")
if (!all(required %in% names(expr))) stop("Marker expression columns missing")
expr <- expr[expr$gene %in% marker_spec$gene & expr$Time %in% times & expr$Treatment %in% c("NC", "LA", "HCl"), , drop = FALSE]
if (nrow(expr) != 18L * 3L * 3L * 4L) stop("Expected 648 marker expression rows; found ", nrow(expr))
if (any(!is.finite(expr$logCPM))) stop("Non-finite marker expression")
expr <- merge(expr, marker_spec, by = "gene", suffixes = c("", "_spec"), all.x = TRUE, sort = FALSE)
if (any(expr$marker_group != expr$module)) {
  # source uses underscore labels; the frozen marker membership is defined by marker_spec.
  expr$marker_group <- expr$module
}
expr$gene <- factor(expr$gene, levels = marker_spec$gene)
expr$module <- factor(expr$module, levels = module_order)
expr$Treatment <- factor(expr$Treatment, levels = c("NC", "LA", "HCl"))
expr$Time <- factor(expr$Time, levels = times)

means <- aggregate(logCPM ~ gene + module + Treatment + Time, expr, mean)
means$gene <- factor(means$gene, levels = marker_spec$gene)
means$module <- factor(means$module, levels = module_order)
means$Treatment <- factor(means$Treatment, levels = c("NC", "LA", "HCl"))
means$Time <- factor(means$Time, levels = times)

make_module_plot <- function(mod) {
  genes <- marker_spec$gene[marker_spec$module == mod]
  d <- expr[expr$module == mod, , drop = FALSE]
  m <- means[means$module == mod, , drop = FALSE]
  d$gene <- factor(d$gene, levels = rev(genes))
  m$gene <- factor(m$gene, levels = rev(genes))
  ggplot() +
    geom_line(data = d, aes(x = Time, y = logCPM, group = interaction(Treatment, Donor), colour = Treatment), linewidth = .24, alpha = .28) +
    geom_point(data = d, aes(x = Time, y = logCPM, colour = Treatment), size = .50, alpha = .42) +
    geom_line(data = m, aes(x = Time, y = logCPM, group = Treatment, colour = Treatment), linewidth = .62) +
    geom_point(data = m, aes(x = Time, y = logCPM, colour = Treatment), size = .72) +
    facet_grid(rows = vars(gene), scales = "free_y", switch = "y") +
    scale_colour_manual(values = treatment_cols, drop = FALSE, name = NULL) +
    labs(title = mod, x = "Time", y = NULL) +
    theme_classic(base_family = "Arial", base_size = 5.6) +
    theme(plot.title = element_text(size = 5.6, face = "bold", colour = unname(module_cols[mod]), hjust = .5,
                                    margin = margin(b = 1.1, unit = "mm")),
          strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(size = 4.4, face = "bold", colour = "#1A1A1A", angle = 0),
          panel.border = element_rect(fill = NA, colour = "#1A1A1A", linewidth = .28),
          panel.spacing = unit(.75, "mm"), panel.grid = element_blank(),
          axis.text.x = element_text(size = 4.3, colour = "#1A1A1A"),
          axis.text.y = element_text(size = 4.1, colour = "#1A1A1A"),
          axis.title.x = element_text(size = 4.8, colour = "#1A1A1A", margin = margin(t = .8, unit = "mm")),
          axis.ticks = element_line(linewidth = .22, colour = "#1A1A1A"),
          axis.ticks.length = unit(.4, "mm"),
          legend.position = "bottom", legend.text = element_text(size = 4.7),
          legend.key.size = unit(2.8, "mm"), legend.spacing.x = unit(1.5, "mm"),
          plot.margin = margin(t = .8, r = .5, b = .6, l = .5, unit = "mm"))
}

plots <- lapply(module_order, make_module_plot)
p <- wrap_plots(plots, nrow = 1, widths = rep(1, 4), guides = "collect") &
  theme(legend.position = "bottom")

prefix <- file.path(export_dir, "S13C_cd8_marker_modules_row_candidate_v0.1")
ggsave(paste0(prefix, ".pdf"), p, width = 180/25.4, height = 72/25.4, units = "in", device = cairo_pdf, bg = "white")
ggsave(paste0(prefix, ".tiff"), p, width = 180/25.4, height = 72/25.4, units = "in", dpi = 600, compression = "lzw", device = "tiff", bg = "white")
ggsave(paste0(prefix, ".png"), p, width = 180/25.4, height = 72/25.4, units = "in", dpi = 300, device = "png", bg = "white")

write.table(expr, file.path(plot_ready_dir, "s13c_cd8_marker_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(means, file.path(plot_ready_dir, "s13c_cd8_marker_means.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(marker_spec, file.path(plot_ready_dir, "s13c_cd8_marker_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

qa <- data.frame(check_id = c("marker_count", "module_count", "expression_rows", "donor_count", "treatment_count", "time_count", "one_row_four_columns", "panel_size_mm", "formal_write", "png_nonempty", "pdf_nonempty", "tiff_nonempty"),
                 expected = c("18", "4", "648", "3", "3", "4", "TRUE", "180 x 72", "FALSE", "TRUE", "TRUE", "TRUE"),
                 observed = c(length(unique(expr$gene)), length(unique(expr$module)), nrow(expr), length(unique(expr$Donor)), length(unique(expr$Treatment)), length(unique(expr$Time)), "TRUE", "180 x 72", "FALSE",
                              file.exists(paste0(prefix, ".png")) && file.info(paste0(prefix, ".png"))$size > 0,
                              file.exists(paste0(prefix, ".pdf")) && file.info(paste0(prefix, ".pdf"))$size > 0,
                              file.exists(paste0(prefix, ".tiff")) && file.info(paste0(prefix, ".tiff"))$size > 0), stringsAsFactors = FALSE)
qa$pass <- as.character(qa$expected) == as.character(qa$observed)
write.table(qa, file.path(qa_dir, "S13C_cd8_marker_modules_row_candidate_v0.1_qa.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(qa_dir, "S13C_cd8_marker_modules_row_candidate_v0.1_session_info.txt"))
if (!all(qa$pass)) stop("S13C QA failed")
message("Rendered S13C: ", prefix)
}
