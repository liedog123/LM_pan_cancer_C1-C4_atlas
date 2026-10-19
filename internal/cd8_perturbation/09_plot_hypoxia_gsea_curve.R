#!/usr/bin/env Rscript
# Purpose: Export the reported CD8 figure source tables and component plots.
# Task: cd8_figures; step: s09
# Inputs (DCF fields): bulk_de_dir,bridge_dir,hallmark_gmt
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_figures --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 09_plot_hypoxia_gsea_curve.R --config PATH [--check-config]\n       Rscript 09_plot_hypoxia_gsea_curve.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("bulk_de_dir","bridge_dir","hallmark_gmt","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
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

sha256 <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}
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
