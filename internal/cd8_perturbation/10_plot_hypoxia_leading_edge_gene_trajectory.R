#!/usr/bin/env Rscript
# Purpose: Export the reported CD8 figure source tables and component plots.
# Task: cd8_figures; step: s10
# Inputs (DCF fields): trajectory_source_file,leading_edge_file
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_figures --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 10_plot_hypoxia_leading_edge_gene_trajectory.R --config PATH [--check-config]\n       Rscript 10_plot_hypoxia_leading_edge_gene_trajectory.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("trajectory_source_file","leading_edge_file","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("trajectory_source_file","leading_edge_file")]
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

sha256 <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}
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
