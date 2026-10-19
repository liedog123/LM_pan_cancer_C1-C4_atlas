#!/usr/bin/env Rscript
# Purpose: Export the reported CD8 figure source tables and component plots.
# Task: cd8_figures; step: s06
# Inputs (DCF fields): bulk_de_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_figures --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 06_plot_hcl_temporal_quadrant.R --config PATH [--check-config]\n       Rscript 06_plot_hcl_temporal_quadrant.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("bulk_de_dir","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- character(0L)
de_basenames <- paste0("de_HCl_minus_NC_",c("1h","6h","24h","48h"),".tsv")
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
