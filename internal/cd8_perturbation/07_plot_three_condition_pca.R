#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 07_plot_three_condition_pca.R --config PATH [--check-config]\n       Rscript 07_plot_three_condition_pca.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("pca_coordinates","pca_variance","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("pca_coordinates","pca_variance")]
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
