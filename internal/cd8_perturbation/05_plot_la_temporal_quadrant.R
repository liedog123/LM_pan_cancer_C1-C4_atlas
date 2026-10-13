#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 05_plot_la_temporal_quadrant.R --config PATH [--check-config]\n       Rscript 05_plot_la_temporal_quadrant.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("bulk_de_dir","leading_edge_file","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["leading_edge_file"]
de_basenames <- paste0("de_LA_minus_NC_",c("1h","6h","24h","48h"),".tsv")
required_files <- c(unname(input_paths),file.path(unname(paths["bulk_de_dir"]),de_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths["bulk_de_dir"])
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

sha256 <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}
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
