#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 11_plot_hcl_gene_trajectory.R --config PATH [--check-config]\n       Rscript 11_plot_hcl_gene_trajectory.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("sample_sheet","filtered_logcpm","bulk_de_dir","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("sample_sheet","filtered_logcpm")]
required_files <- c(unname(input_paths),file.path(unname(paths["bulk_de_dir"]),paste0("de_HCl_minus_NC_",c("1h","6h","24h","48h"),".tsv")))

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

rows <- list(); idx <- 0L
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
