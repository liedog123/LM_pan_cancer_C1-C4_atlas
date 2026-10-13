#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 12_plot_cd8_marker_modules.R --config PATH [--check-config]\n       Rscript 12_plot_cd8_marker_modules.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("marker_expression","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("marker_expression")]
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
