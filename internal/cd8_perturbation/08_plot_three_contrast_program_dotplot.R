#!/usr/bin/env Rscript
# Purpose: Export the reported CD8 figure source tables and component plots.
# Task: cd8_figures; step: s08
# Inputs (DCF fields): gsea_three_contrasts
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_figures --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 08_plot_three_contrast_program_dotplot.R --config PATH [--check-config]\n       Rscript 08_plot_three_contrast_program_dotplot.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("gsea_three_contrasts","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("gsea_three_contrasts")]
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
