#!/usr/bin/env Rscript
# Purpose: Aggregate significant gene-level copy-number differences into descriptive genomic regions.
# Task: cnv_regions; step: s05
# Inputs (DCF fields): manifest,thresholded_cnv,gene_coordinates,cytoband
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cnv_regions --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 05_run_cnv_alteration_regions.R --config PATH [--check-config]\n       Rscript 05_run_cnv_alteration_regions.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("manifest", "thresholded_cnv", "gene_coordinates", "cytoband", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),config_path,package_root)
if(out_dir=="/" || out_dir==normalizePath(path.expand("~")) ||
   any(vapply(protected,function(p) is_within(p,out_dir)||is_within(out_dir,p),logical(1L)))) stop("Output must be independent of inputs/config.")
if(file.exists(out_dir) && !dir.exists(out_dir)) stop("Output is not a directory.")
if(dir.exists(out_dir) && (file.access(out_dir,7L)!=0L || length(list.files(out_dir,all.files=TRUE,no..=TRUE)))) stop("Output must be writable and empty.")
parent <- out_dir
while(!dir.exists(parent)) {
 if(file.exists(parent)) stop("Output has a non-directory ancestor.")
 parent <- dirname(parent)
}
if(file.access(parent,3L)!=0L) stop("No writable searchable output parent.")
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})
# 四项输入与独立输出根参数化；科学和图形主体保持来源实现。
manifest_file <- unname(paths["manifest"])
cnv_file <- unname(paths["thresholded_cnv"])
probe_file <- unname(paths["gene_coordinates"])
cytoband_file <- unname(paths["cytoband"])
source_dir <- file.path(out_dir, "source_data")
export_dir <- file.path(out_dir, "panel_exports")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

MIN_N <- 30L
ALT_FREQ_MIN <- 0.01
FDR_MAX <- 0.05
DIFF_MIN <- 0.15
MAX_GAP <- 10e6
TOP_N <- 15L
STAR_GENES <- c("EGFR", "BRAF", "TP53", "KRAS", "MYC", "PIK3CA", "PTEN",
                "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")

message("[1/7] Read C1/C4 manifest and selected CNV matrix for an unfrozen public rerun")
manifest <- fread(manifest_file)
stopifnot(!anyDuplicated(manifest$Patient_ID), all(manifest$C_subtype %chin% c("C1", "C4")))
coverage <- dcast(manifest[, .N, by = .(Cancer, C_subtype)], Cancer ~ C_subtype,
                  value.var = "N", fill = 0L)
coverage[, formal_eligible := C1 >= MIN_N & C4 >= MIN_N]
eligible_cancers <- coverage[formal_eligible == TRUE, Cancer]
if (!length(eligible_cancers)) stop("No eligible cancers")

thr <- fread(cnv_file, select = c("Sample", manifest$Sample15), showProgress = TRUE)
genes <- as.character(thr[[1L]])
if (anyDuplicated(genes)) stop("Duplicate genes in CNV matrix")
mat <- as.matrix(thr[, -1L])
storage.mode(mat) <- "integer"
if (anyNA(mat) || any(!unique(as.vector(mat)) %in% -2:2)) stop("Invalid thresholded CNV values")
rm(thr); invisible(gc())
manifest[, matrix_col := match(Sample15, colnames(mat))]
if (anyNA(manifest$matrix_col)) stop("Manifest samples missing from matrix")

# 与旧版相同的低变异预筛思想；在当前冻结C1/C4样本宇宙中计算。
alt_freq <- rowMeans(mat != 0)
keep_gene <- alt_freq > ALT_FREQ_MIN
mat <- mat[keep_gene, , drop = FALSE]
genes <- genes[keep_gene]
message("Genes after alt_freq > 0.01: ", length(genes))

fisher_vec <- function(a, b, n1, n4) {
  vapply(seq_along(a), function(i) {
    if (a[i] == 0L && b[i] == 0L) return(1)
    fisher.test(matrix(c(a[i], n1 - a[i], b[i], n4 - b[i]), nrow = 2L))$p.value
  }, numeric(1L))
}

message("[2/7] Run within-cancer C1 vs C4 Fisher tests")
all_stats <- vector("list", length(eligible_cancers))
for (ii in seq_along(eligible_cancers)) {
  cancer <- eligible_cancers[ii]
  idx1 <- manifest[Cancer == cancer & C_subtype == "C1", matrix_col]
  idx4 <- manifest[Cancer == cancer & C_subtype == "C4", matrix_col]
  n1 <- length(idx1); n4 <- length(idx4)
  c1_amp <- rowSums(mat[, idx1, drop = FALSE] > 0)
  c4_amp <- rowSums(mat[, idx4, drop = FALSE] > 0)
  c1_del <- rowSums(mat[, idx1, drop = FALSE] < 0)
  c4_del <- rowSums(mat[, idx4, drop = FALSE] < 0)
  amp <- data.table(Gene = genes, Cancer = cancer, Type = "Amp", n_C1 = n1, n_C4 = n4,
                    event_C1 = c1_amp, event_C4 = c4_amp,
                    Freq_C1 = c1_amp / n1, Freq_C4 = c4_amp / n4)
  del <- data.table(Gene = genes, Cancer = cancer, Type = "Del", n_C1 = n1, n_C4 = n4,
                    event_C1 = c1_del, event_C4 = c4_del,
                    Freq_C1 = c1_del / n1, Freq_C4 = c4_del / n4)
  amp[, `:=`(Diff = Freq_C1 - Freq_C4, P = fisher_vec(event_C1, event_C4, n1, n4))]
  del[, `:=`(Diff = Freq_C1 - Freq_C4, P = fisher_vec(event_C1, event_C4, n1, n4))]
  amp[, FDR := p.adjust(P, "BH")]
  del[, FDR := p.adjust(P, "BH")]
  z <- rbind(amp, del)
  z[, `:=`(Abs_Diff = abs(Diff), Bias = fifelse(Diff > 0, "C1 higher", "C4 higher"))]
  z[, significant := FDR < FDR_MAX & Abs_Diff > DIFF_MIN]
  all_stats[[ii]] <- z
  message(sprintf("  %s (%d/%d): Amp=%d, Del=%d", cancer, ii, length(eligible_cancers),
                  amp[FDR < FDR_MAX & abs(Diff) > DIFF_MIN, .N],
                  del[FDR < FDR_MAX & abs(Diff) > DIFF_MIN, .N]))
}
stats <- rbindlist(all_stats)
fwrite(stats, file.path(source_dir, "fig4_d_cnv_legacy_rule_all_gene_statistics.tsv.gz"),
       sep = "\t", compress = "gzip")
sig <- stats[significant == TRUE]
if (!nrow(sig)) stop("No significant events under legacy-adapted rules")

message("[3/7] Map significant genes to GENCODE v23 coordinates")
probe <- fread(probe_file)
probe <- probe[chrom %chin% paste0("chr", c(1:22, "X", "Y"))]
probe[, symbol_key := sub("\\|.*$", "", gene)]
coords <- probe[, .(Chr = unique(chrom)[1L], Start = min(chromStart), End = max(chromEnd),
                    n_chrom = uniqueN(chrom)), by = symbol_key]
coords <- coords[n_chrom == 1L]
setnames(coords, "symbol_key", "Gene")
mapped <- merge(sig, coords[, .(Gene, Chr, Start, End)], by = "Gene", all = FALSE)
mapping_summary <- data.table(metric = c("significant_gene_events", "mapped_gene_events", "mapping_fraction"),
                              value = c(nrow(sig), nrow(mapped), nrow(mapped) / nrow(sig)))
if (nrow(mapped) / nrow(sig) < 0.90) stop("Coordinate mapping below 90%")

message("[4/7] Merge significant genes into <=10 Mb alteration regions")
setorder(mapped, Cancer, Type, Bias, Chr, Start)
mapped[, new_cluster := as.integer(seq_len(.N) == 1L | Start - shift(End, fill = -1e9) > MAX_GAP),
       by = .(Cancer, Type, Bias, Chr)]
mapped[, cluster_id := cumsum(new_cluster), by = .(Cancer, Type, Bias, Chr)]
regions <- mapped[, {
  stars <- intersect(Gene, STAR_GENES)
  representative <- if (length(stars)) stars[1L] else Gene[which.max(Abs_Diff)]
  .(Start = min(Start), End = max(End), Mean_Diff = mean(Abs_Diff), Max_Diff = max(Abs_Diff),
    Gene_Count = .N, Representative_Gene = representative,
    Star_Genes = paste(stars, collapse = ";"), Is_Star = length(stars) > 0,
    Min_FDR = min(FDR))
}, by = .(Cancer, Type, Bias, Chr, cluster_id)]

message("[5/7] Apply old Top-15 plus star-region display rule")
display_regions <- regions[, {
  stars <- .SD[Is_Star == TRUE]
  ordinary <- .SD[Is_Star == FALSE][order(-Mean_Diff)][seq_len(min(.N, TOP_N))]
  unique(rbind(stars, ordinary, fill = TRUE))
}, by = Cancer]
setorder(display_regions, Cancer, Chr, Start, Type, Bias)

message("[6/7] Build two-column hg38 karyotype")
left_chroms <- paste0("chr", c(1:11, "X"))
right_chroms <- paste0("chr", c(22:12, "Y"))
chrom_levels <- c(left_chroms, right_chroms)
cyto <- fread(cytoband_file,
              col.names = c("chrom", "chromStart", "chromEnd", "band", "gieStain"))
cyto <- cyto[chrom %chin% chrom_levels]
sizes <- cyto[, .(chromLength = max(chromEnd)), by = chrom]
layout <- rbind(
  data.table(chrom = left_chroms, panel_side = "left", row_index = rev(seq_along(left_chroms))),
  data.table(chrom = right_chroms, panel_side = "right", row_index = rev(seq_along(right_chroms)))
)
layout <- merge(layout, sizes, by = "chrom", all.x = TRUE, sort = FALSE)
layout[, row_y := row_index * 0.52]
left_max <- max(layout[panel_side == "left", chromLength])
right_max <- max(layout[panel_side == "right", chromLength])
gap <- 70e6
right_end <- left_max + gap + right_max
to_x <- function(pos, side, len) ifelse(side == "left", pos, right_end - len + pos)
cyto <- merge(cyto, layout, by = "chrom", all.x = TRUE, sort = FALSE)
cyto[, `:=`(x_start = to_x(chromStart, panel_side, chromLength),
            x_end = to_x(chromEnd, panel_side, chromLength))]
outline <- copy(layout)
outline[, `:=`(x_start = fifelse(panel_side == "left", 0, right_end - chromLength),
               x_end = fifelse(panel_side == "left", chromLength, right_end))]
display_regions <- merge(display_regions, layout, by.x = "Chr", by.y = "chrom", all.x = TRUE)
display_regions[, `:=`(mid_pos = (Start + End) / 2,
                       x = to_x((Start + End) / 2, panel_side, chromLength),
                       y = row_y + fifelse(Bias == "C1 higher", 0.115, -0.115))]

labels <- display_regions[Is_Star == TRUE & nzchar(Star_Genes)]
labels[, Label := sub(";.*$", "", Star_Genes)]
setorder(labels, Label, -Mean_Diff)
labels <- labels[, .SD[1L], by = Label]

stain_cols <- c(gneg = "#FFFFFF", gpos25 = "#D9D9D9", gpos50 = "#A6A6A6",
                gpos75 = "#737373", gpos100 = "#111111", acen = "#111111",
                gvar = "#D9D9D9", stalk = "#D9D9D9")
cyto[!gieStain %chin% names(stain_cols), gieStain := "gneg"]
cancers <- sort(unique(display_regions$Cancer))
shape_pool <- c(0:14, 15:18)
if (length(cancers) > length(shape_pool)) stop("Too many cancers for stable shape map")
shape_values <- setNames(shape_pool[seq_along(cancers)], cancers)
axis_breaks <- sort(unique(layout$row_y))
left_axis_labels <- layout[panel_side == "left"][order(row_y), sub("^chr", "", chrom)]
right_axis_labels <- layout[panel_side == "right"][order(row_y), sub("^chr", "", chrom)]

p <- ggplot() +
  geom_segment(data = outline, aes(x = x_start, xend = x_end, y = row_y, yend = row_y),
               linewidth = 2.0, lineend = "round", color = "black") +
  geom_segment(data = cyto, aes(x = x_start, xend = x_end, y = row_y, yend = row_y, color = gieStain),
               linewidth = 1.45, lineend = "butt", show.legend = FALSE) +
  geom_point(data = display_regions, aes(x = x, y = y, shape = Cancer, color = Bias),
             size = 2.15, stroke = 0.62, fill = NA, alpha = 0.94) +
  geom_text_repel(data = labels, aes(x = x, y = y, label = Label),
                  size = 2.15, color = "black", box.padding = 0.24, point.padding = 0.12,
                  min.segment.length = 0, segment.size = 0.18, max.overlaps = Inf,
                  seed = 20260809) +
  scale_color_manual(values = c(stain_cols, "C1 higher" = "#FF2D2D", "C4 higher" = "#1764FF"),
                     breaks = c("C1 higher", "C4 higher"), name = "Enrichment bias") +
  scale_shape_manual(values = shape_values, name = NULL) +
  scale_y_continuous(breaks = axis_breaks, labels = left_axis_labels,
                     sec.axis = dup_axis(breaks = axis_breaks, labels = right_axis_labels, name = NULL),
                     expand = expansion(add = c(0.34, 0.10))) +
  scale_x_continuous(expand = expansion(add = c(8e6, 8e6))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Chromosome") +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(axis.title.y = element_text(size = 7.0, margin = margin(r = 4)),
        axis.text.y = element_text(size = 6.0, color = "black"),
        axis.ticks.y = element_line(linewidth = 0.3), axis.ticks.length.y = unit(1.4, "mm"),
        axis.line.y = element_line(linewidth = 0.3),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.position = "right", legend.box = "vertical", legend.text = element_text(size = 5.8),
        legend.title = element_text(size = 6.0), legend.key.height = unit(7.2, "pt"),
        plot.margin = margin(2.0, 2.5, 2.0, 5.0, unit = "mm")) +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.2)),
         shape = guide_legend(order = 2, override.aes = list(color = "black", size = 2.0)))

message("[7/7] Write unfrozen public rerun source data, QA and exports")
fwrite(display_regions, file.path(source_dir, "fig4_d_cnv_chromosome_map.tsv"), sep = "\t")
fwrite(labels, file.path(source_dir, "fig4_d_cnv_chromosome_labels.tsv"), sep = "\t")
fwrite(layout, file.path(source_dir, "fig4_d_cnv_chromosome_layout.tsv"), sep = "\t")
fwrite(coverage, file.path(source_dir, "fig4_d_cnv_legacy_rule_coverage.tsv"), sep = "\t")
fwrite(mapping_summary, file.path(source_dir, "fig4_d_cnv_legacy_rule_mapping_summary.tsv"), sep = "\t")
qa <- data.table(
  check = c("patient_unique", "only_c1_c4", "min_group_n", "within_cancer_bh",
            "formal_threshold", "mapping_ge90pct", "no_cross_cancer_pooling"),
  pass = c(!anyDuplicated(manifest$Patient_ID), all(manifest$C_subtype %chin% c("C1", "C4")),
           all(stats$n_C1 >= 30 & stats$n_C4 >= 30), TRUE,
           all(sig$FDR < 0.05 & sig$Abs_Diff > 0.15), nrow(mapped) / nrow(sig) >= 0.90, TRUE)
)
fwrite(qa, file.path(source_dir, "fig4_d_cnv_legacy_rule_qa.tsv"), sep = "\t")
if (!all(qa$pass)) stop("QA failed")
manifest_out <- data.table(
  field = c("panel", "method", "comparison", "eligible_cancers", "significant_gene_events",
            "mapped_gene_events", "display_regions", "display_cancers", "label_genes",
            "width_mm", "height_mm", "evidence_boundary"),
  value = c("D", "legacy CNV landscape rules adapted to canonical C1 vs C4", "within-cancer C1 vs C4",
            length(eligible_cancers), nrow(sig), nrow(mapped), nrow(display_regions), length(cancers),
            paste(labels$Label, collapse = ";"), 180, 72,
            "10 Mb alteration clusters; not GISTIC focal peaks or arm-level events")
)
fwrite(manifest_out, file.path(source_dir, "fig4_d_cnv_chromosome_map_manifest.tsv"), sep = "\t")

stem <- file.path(export_dir, "fig4_d_cnv_chromosome_map")
ggsave(paste0(stem, ".pdf"), p, width = 180, height = 72, units = "mm", device = cairo_pdf)
ggsave(paste0(stem, ".tiff"), p, width = 180, height = 72, units = "mm", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(paste0(stem, ".png"), p, width = 180, height = 72, units = "mm", dpi = 300,
       device = "png", bg = "white")
writeLines(capture.output(sessionInfo()), file.path(source_dir, "fig4_d_cnv_legacy_rule_sessionInfo.txt"))
message(sprintf("DONE: %d significant gene-events, %d display regions, %d cancers, %d labels",
                nrow(sig), nrow(display_regions), length(cancers), nrow(labels)))
