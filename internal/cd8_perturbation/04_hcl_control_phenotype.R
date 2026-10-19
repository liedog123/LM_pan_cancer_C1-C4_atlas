#!/usr/bin/env Rscript
# Purpose: Fit donor-adjusted CD8 treatment-time responses and fixed programme enrichment.
# Task: cd8_response; step: s04
# Inputs (DCF fields): count_rdata,lm12_membership,msigdb_dir,bulk_de_dir,bridge_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_response --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 04_hcl_control_phenotype.R --config PATH [--check-config]\n       Rscript 04_hcl_control_phenotype.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("count_rdata","lm12_membership","msigdb_dir","bulk_de_dir","bridge_dir","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("count_rdata","lm12_membership")]
gmt_basenames <- c("h.all.v2025.1.Hs.symbols.gmt","c2.cp.reactome.v2025.1.Hs.symbols.gmt","c5.go.bp.v2024.1.Hs.symbols.gmt","c7.all.v2025.1.Hs.symbols.gmt")
de_contrasts <- c("LA_minus_NC","LA_minus_HCl","HCl_minus_NC")
de_basenames <- as.vector(outer(de_contrasts, c("1h","6h","24h","48h"), function(x,y) paste0("de_",x,"_",y,".tsv")))
required_files <- c(unname(input_paths), file.path(unname(paths["msigdb_dir"]),gmt_basenames), file.path(unname(paths["bulk_de_dir"]),de_basenames))
bridge_basenames <- c("fixed_gene_set_manifest.tsv","bulk_all_fixed_preranked_gsea.tsv","bulk_fixed_readout_ssgsea_contrasts.tsv","bulk_fixed_readout_donor_summary.tsv")
required_files <- c(required_files,file.path(unname(paths["bridge_dir"]),bridge_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("msigdb_dir","bulk_de_dir","bridge_dir")])
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

cfg_paths <- paths

# ---- Ordinary-only fixed-resource path routing; no set/content changes ----
relocate_fixed_source <- function(layer, source_path) {
  if (length(layer) != 1L || is.na(layer) || length(source_path) != 1L ||
      is.na(source_path) || !nzchar(source_path)) stop("Invalid fixed-resource identity.")
  if (layer == "LM12_FIG2_FROZEN") return(unname(cfg_paths["lm12_membership"]))
  if (layer != "MSIGDB_FULL_SET") stop("Unknown fixed-resource layer: ", layer)
  file_name <- basename(source_path)
  if (!file_name %in% gmt_basenames) stop("Unregistered fixed GMT basename: ", file_name)
  file.path(unname(cfg_paths["msigdb_dir"]), file_name)
}

# ---- Retained complete source body; documented path adaptations only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(fgsea)
})

analysis_seed <- 20260916L
input_path <- unname(cfg_paths["count_rdata"])
de_root <- unname(cfg_paths["bulk_de_dir"])
bridge_root <- unname(cfg_paths["bridge_dir"])
out_dir <- unname(cfg_paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = TRUE, na = "NA")
stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)
times <- c("1h", "6h", "24h", "48h")
contrasts <- c("LA_minus_NC", "HCl_minus_NC", "LA_minus_HCl")

load_count <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  stop_if(!all(c("count", "sample_anno") %in% ls(e, all.names = TRUE)), "count/sample_anno missing")
  list(count = get("count", e), anno = as.data.table(get("sample_anno", e), keep.rownames = "sample_id"))
}

obj <- load_count(input_path)
count <- obj$count
anno <- obj$anno
stop_if(!identical(colnames(count), anno$sample_id), "count and annotation order mismatch")
anno$sample_id <- as.character(anno$sample_id)
anno$Donor <- as.character(anno$Donor)
anno$Treatment <- as.character(anno$Treatment)
anno$Time <- as.character(anno$Time)
observed_levels <- rbindlist(list(
  data.table(level_type = "Treatment", level = paste(sort(unique(anno$Treatment)), collapse = ";")),
  data.table(level_type = "Time", level = paste(sort(unique(anno$Time)), collapse = ";")),
  data.table(level_type = "Donor", level = paste(sort(unique(anno$Donor)), collapse = ";"))
))
write_tsv(observed_levels, file.path(out_dir, "observed_design_levels.tsv"))

factorial <- anno[Treatment %in% c("NC", "HCl", "LA") & Time %in% times]
factorial$Treatment <- factor(factorial$Treatment, levels = c("NC", "HCl", "LA"))
factorial$Time <- factor(factorial$Time, levels = times)
factorial$Donor <- factor(factorial$Donor)
stop_if(nrow(factorial) != 36L, paste("unexpected factorial sample count:", nrow(factorial)))
stop_if(anyDuplicated(factorial[, .(Donor, Treatment, Time)]) > 0L, "duplicate donor-treatment-time sample")
design <- model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = factorial)
stop_if(qr(design)$rank != ncol(design), "factorial design is rank deficient")
write_tsv(factorial[, .(sample_id, Donor = as.character(Donor), Treatment = as.character(Treatment), Time = as.character(Time))],
          file.path(out_dir, "factorial_sample_manifest.tsv"))

fixed_manifest <- fread(file.path(bridge_root, "fixed_gene_set_manifest.tsv"))
stop_if(!all(c("feature_id", "feature_label", "feature_layer", "source_path") %in% names(fixed_manifest)),
        "fixed gene-set manifest columns missing")
# Ordinary run: map historical manifest resource paths to configured inputs.
fixed_manifest$source_path <- vapply(seq_len(nrow(fixed_manifest)), function(i)
  relocate_fixed_source(as.character(fixed_manifest$feature_layer[i]),
                        as.character(fixed_manifest$source_path[i])), character(1L))
lm12_path <- fixed_manifest[feature_layer == "LM12_FIG2_FROZEN", unique(source_path)][1L]
lm12_tab <- fread(lm12_path)
stop_if(!all(c("feature_id", "gene") %in% names(lm12_tab)), "LM12 membership columns missing")
lm12_ids <- unique(as.character(lm12_tab$feature_id))
lm12_sets <- lapply(lm12_ids, function(id) unique(as.character(lm12_tab[feature_id == id, gene])))
names(lm12_sets) <- paste0("LM12::", lm12_ids)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  toks <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(toks, function(x) x[[1L]], character(1))
  sets <- lapply(toks, function(x) unique(x[-seq_len(min(2L, length(x)))]))
  names(sets) <- ids
  sets[!duplicated(names(sets))]
}
msig_rows <- fixed_manifest[feature_layer == "MSIGDB_FULL_SET"]
gmt_cache <- list()
msig_sets <- list()
for (i in seq_len(nrow(msig_rows))) {
  path <- msig_rows$source_path[i]
  if (is.null(gmt_cache[[path]])) gmt_cache[[path]] <- read_gmt(path)
  raw_id <- msig_rows$feature_id[i]
  stop_if(!raw_id %in% names(gmt_cache[[path]]), paste("MSigDB set missing:", raw_id))
  msig_sets[[paste0("MSIGDB::", raw_id)]] <- gmt_cache[[path]][[raw_id]]
}
all_sets <- c(lm12_sets, msig_sets)
set_meta <- rbind(
  fixed_manifest[feature_layer == "LM12_FIG2_FROZEN", .(
    feature_id = paste0("LM12::", feature_id), feature_label, feature_layer, source_collection,
    n_genes_source, source_path
  )],
  fixed_manifest[feature_layer == "MSIGDB_FULL_SET", .(
    feature_id = paste0("MSIGDB::", feature_id), feature_label, feature_layer, source_collection,
    n_genes_source, source_path
  )]
)
set_meta <- set_meta[match(names(all_sets), feature_id)]
stop_if(anyNA(set_meta$feature_id) || length(all_sets) != nrow(set_meta), "set metadata mapping failed")
write_tsv(set_meta, file.path(out_dir, "fixed_gene_set_manifest_reused.tsv"))

load_de <- function(contrast, tm) {
  path <- file.path(de_root, paste0("de_", contrast, "_", tm, ".tsv"))
  stop_if(!file.exists(path), paste("missing DE input:", path))
  d <- fread(path)
  stop_if(!all(c("gene_id", "logFC", "PValue", "FDR") %in% names(d)), paste("DE columns missing:", path))
  d$rank_score <- sign(d$logFC) * -log10(pmax(d$PValue, .Machine$double.xmin))
  d <- d[is.finite(d$rank_score) & !duplicated(d$gene_id)]
  stats <- d$rank_score
  names(stats) <- d$gene_id
  list(table = d, stats = stats)
}

flatten_gsea <- function(x) {
  x <- as.data.table(x)
  if (nrow(x) > 0L && "leadingEdge" %in% names(x)) {
    x$leadingEdge_genes <- vapply(x$leadingEdge, function(g) paste(as.character(g), collapse = ";"), character(1))
    x$leadingEdge <- NULL
  }
  x
}

run_gsea_one <- function(stats, contrast, tm) {
  present <- lapply(all_sets, function(g) intersect(unique(g), names(stats)))
  n_present <- lengths(present)
  testable <- n_present >= 3L & n_present <= 5000L
  out <- data.table()
  if (any(testable)) {
    set.seed(analysis_seed)
    fg <- fgsea::fgseaMultilevel(pathways = all_sets[testable], stats = stats,
                                  minSize = 3L, maxSize = 5000L, eps = 0)
    out <- flatten_gsea(fg)
    out$feature_id <- as.character(out$pathway)
    out$pathway <- NULL
    out$contrast <- contrast
    out$time <- tm
    out$n_genes_tested <- n_present[match(out$feature_id, names(all_sets))]
    out$status <- "PASS"
  }
  hold_ids <- names(all_sets)[!testable]
  if (length(hold_ids)) {
    hold <- data.table(feature_id = hold_ids, pval = NA_real_, padj = NA_real_, log2err = NA_real_,
                       ES = NA_real_, NES = NA_real_, size = n_present[hold_ids],
                       contrast = contrast, time = tm, n_genes_tested = n_present[hold_ids],
                       status = "HOLD", leadingEdge_genes = NA_character_)
    out <- rbindlist(list(out, hold), fill = TRUE)
  }
  out <- merge(out, set_meta[, .(feature_id, feature_label, feature_layer, source_collection,
                                 n_genes_source)], by = "feature_id", all.x = TRUE, sort = FALSE)
  out[]
}

existing_la <- fread(file.path(bridge_root, "bulk_all_fixed_preranked_gsea.tsv"))
existing_la <- existing_la[context == "bulk_LA_minus_NC" & time %in% times]
existing_la$contrast <- "LA_minus_NC"
existing_la$status <- "REUSED_FROZEN"
existing_la$n_genes_source <- set_meta$n_genes_source[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$feature_label <- set_meta$feature_label[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$feature_layer <- set_meta$feature_layer[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$source_collection <- set_meta$source_collection[match(existing_la$feature_id, set_meta$feature_id)]
existing_la$n_genes_tested <- existing_la$size
existing_la <- existing_la[, .(contrast, time, feature_id, feature_label, feature_layer,
                               source_collection, n_genes_source, n_genes_tested, status,
                               pval, padj, log2err, ES, NES, size, leadingEdge_genes)]

new_gsea <- rbindlist(lapply(contrasts[contrasts != "LA_minus_NC"], function(ct) {
  rbindlist(lapply(times, function(tm) run_gsea_one(load_de(ct, tm)$stats, ct, tm)), fill = TRUE)
}), fill = TRUE)
gsea_three <- rbindlist(list(existing_la, new_gsea), fill = TRUE)
gsea_three <- gsea_three[order(match(contrast, contrasts), match(time, times),
                               match(feature_id, set_meta$feature_id))]
write_tsv(gsea_three, file.path(out_dir, "gsea_three_contrasts.tsv"))
write_tsv(new_gsea, file.path(out_dir, "gsea_new_control_contrasts.tsv"))

reuse_ref <- fread(file.path(bridge_root, "bulk_all_fixed_preranked_gsea.tsv"))[
  context == "bulk_LA_minus_NC" & time %in% times,
  .(feature_id, time, source_NES = NES, source_FDR = padj, source_ES = ES)
]
reuse_check <- merge(existing_la[, .(feature_id, time, reused_NES = NES, reused_FDR = padj, reused_ES = ES)],
                     reuse_ref, by = c("feature_id", "time"), all = TRUE)
reuse_check$NES_delta <- reuse_check$reused_NES - reuse_check$source_NES
reuse_check$FDR_delta <- reuse_check$reused_FDR - reuse_check$source_FDR
reuse_check$ES_delta <- reuse_check$reused_ES - reuse_check$source_ES
reuse_check$exact_reuse <- (is.na(reuse_check$reused_NES) == is.na(reuse_check$source_NES)) &
  (is.na(reuse_check$reused_NES) | abs(reuse_check$NES_delta) < 1e-12) &
  (is.na(reuse_check$reused_FDR) | abs(reuse_check$FDR_delta) < 1e-12) &
  (is.na(reuse_check$reused_ES) | abs(reuse_check$ES_delta) < 1e-12)
write_tsv(reuse_check, file.path(out_dir, "la_nc_reuse_check.tsv"))

score <- fread(file.path(bridge_root, "bulk_fixed_readout_ssgsea_contrasts.tsv"))[
  contrast %in% contrasts & time %in% times]
donor <- fread(file.path(bridge_root, "bulk_fixed_readout_donor_summary.tsv"))[
  contrast %in% contrasts & time %in% times]
stop_if(nrow(score) != 21L * 4L * 3L, paste("unexpected score rows:", nrow(score)))
stop_if(nrow(donor) != 21L * 4L * 3L, paste("unexpected donor score rows:", nrow(donor)))
write_tsv(score, file.path(out_dir, "program_score_three_contrasts.tsv"))
write_tsv(donor, file.path(out_dir, "program_donor_effects_three_contrasts.tsv"))

key_labels <- c("Hypoxia", "Glycolysis", "Unfolded protein response", "TCR signaling",
                "T-cell mediated cytotoxicity", "Effector-enriched (Goldrath effector vs memory, UP)",
                "Memory-enriched (Goldrath effector vs memory, DN)",
                "Antiviral mechanism by IFN-stimulated genes",
                "Class I MHC mediated antigen processing & presentation")
key_score <- score[feature_label %in% key_labels]
key_donor <- donor[feature_id %in% unique(key_score$feature_id)]
write_tsv(key_score, file.path(out_dir, "key_program_score_summary.tsv"))
write_tsv(key_donor, file.path(out_dir, "key_program_donor_summary.tsv"))

marker_groups <- list(
  activation = c("CD69", "IL2RA"),
  cytotoxic_effector = c("GZMB", "PRF1", "NKG7", "GZMA", "IFNG"),
  checkpoint_exhaustion_associated = c("PDCD1", "TIGIT", "LAG3", "CTLA4", "HAVCR2", "TOX"),
  memory_associated = c("CCR7", "IL7R", "TCF7", "SELL", "LEF1")
)
marker_manifest <- rbindlist(lapply(names(marker_groups), function(g) data.table(
  marker_group = g, gene = marker_groups[[g]], prespecified = TRUE
)))
marker_manifest$in_count <- marker_manifest$gene %in% rownames(count)
marker_manifest$n_samples_detected <- vapply(marker_manifest$gene, function(g) {
  if (!g %in% rownames(count)) return(0L)
  sum(count[g, ] > 0)
}, integer(1))
marker_manifest$detection_status <- ifelse(!marker_manifest$in_count, "NOT_IN_COUNT",
                                            ifelse(marker_manifest$n_samples_detected == 0L, "NO_COUNTS", "PRESENT"))
write_tsv(marker_manifest, file.path(out_dir, "cd8_marker_manifest.tsv"))

fac_ids <- factorial$sample_id
dge <- edgeR::DGEList(counts = count[, fac_ids, drop = FALSE], samples = factorial)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
expr <- rbindlist(lapply(marker_manifest$gene, function(g) {
  if (!g %in% rownames(logcpm)) return(data.table(gene = g))
  data.table(gene = g, sample_id = fac_ids, logCPM = as.numeric(logcpm[g, fac_ids]))
}), fill = TRUE)
expr <- merge(expr, factorial[, .(sample_id, Donor = as.character(Donor),
                                  Treatment = as.character(Treatment), Time = as.character(Time))],
              by = "sample_id", all.x = TRUE, sort = FALSE)
expr <- merge(expr, marker_manifest[, .(marker_group, gene, detection_status)],
              by = "gene", all.x = TRUE, sort = FALSE)
write_tsv(expr, file.path(out_dir, "cd8_marker_expression.tsv"))

marker_de <- rbindlist(lapply(contrasts, function(ct) {
  rbindlist(lapply(times, function(tm) {
    d <- load_de(ct, tm)$table
    x <- data.table(gene = marker_manifest$gene)
    x <- merge(x, d[, .(gene = gene_id, logFC, PValue, FDR)], by = "gene", all.x = TRUE, sort = FALSE)
    x$contrast <- ct
    x$time <- tm
    x[]
  }), fill = TRUE)
}), fill = TRUE)
marker_de <- merge(marker_de, marker_manifest[, .(marker_group, gene, detection_status, n_samples_detected)],
                   by = "gene", all.x = TRUE, sort = FALSE)
marker_de$de_status <- ifelse(marker_de$detection_status != "PRESENT", marker_de$detection_status,
                              ifelse(is.na(marker_de$logFC), "FILTERED_FROM_DE", "EVALUABLE"))
marker_de <- marker_de[order(match(marker_group, names(marker_groups)), match(gene, marker_manifest$gene),
                             match(contrast, contrasts), match(time, times))]
write_tsv(marker_de, file.path(out_dir, "cd8_marker_de.tsv"))

marker_donor_rows <- rbindlist(lapply(contrasts, function(ct) {
  pair <- switch(ct,
                 LA_minus_NC = c("LA", "NC"),
                 HCl_minus_NC = c("HCl", "NC"),
                 LA_minus_HCl = c("LA", "HCl"))
  rbindlist(lapply(times, function(tm) {
    left <- expr[Treatment == pair[1] & Time == tm & !is.na(logCPM),
                 .(Donor, gene, left_logCPM = logCPM)]
    right <- expr[Treatment == pair[2] & Time == tm & !is.na(logCPM),
                  .(Donor, gene, right_logCPM = logCPM)]
    x <- merge(left, right, by = c("Donor", "gene"), all = FALSE)
    x$contrast <- ct
    x$time <- tm
    x$donor_delta_logCPM <- x$left_logCPM - x$right_logCPM
    x[, .(contrast, time, Donor, gene, left_logCPM, right_logCPM, donor_delta_logCPM)]
  }), fill = TRUE)
}), fill = TRUE)
marker_donor_rows <- merge(marker_donor_rows,
                            marker_manifest[, .(marker_group, gene, detection_status)],
                            by = "gene", all.x = TRUE, sort = FALSE)
write_tsv(marker_donor_rows, file.path(out_dir, "cd8_marker_donor_effects.tsv"))
marker_donor_summary <- marker_donor_rows[, .(
  n_donors = uniqueN(Donor),
  mean_delta_logCPM = mean(donor_delta_logCPM),
  median_delta_logCPM = median(donor_delta_logCPM),
  sd_delta_logCPM = sd(donor_delta_logCPM),
  n_positive = sum(donor_delta_logCPM > 0),
  n_negative = sum(donor_delta_logCPM < 0),
  n_zero = sum(donor_delta_logCPM == 0)
), by = .(contrast, time, gene, marker_group, detection_status)]
marker_donor_summary$all_donors_same_direction <-
  marker_donor_summary$n_positive == marker_donor_summary$n_donors |
  marker_donor_summary$n_negative == marker_donor_summary$n_donors
write_tsv(marker_donor_summary, file.path(out_dir, "cd8_marker_donor_summary.tsv"))

qa <- data.table(
  check_id = c("observed_treatment_levels", "factorial_samples", "design_full_rank",
               "three_contrasts", "fixed_set_count", "new_gsea_rows", "la_nc_exact_reuse",
               "score_rows", "donor_rows", "marker_rows", "marker_expression_rows",
               "marker_donor_rows", "marker_donor_summary_rows",
               "no_unexpected_fourth_condition"),
  status = c(
    ifelse(identical(sort(unique(anno$Treatment)), sort(c("HCl", "LA", "NC", "Unsti"))), "PASS", "HOLD"),
    ifelse(nrow(factorial) == 36L, "PASS", "HOLD"),
    ifelse(qr(design)$rank == ncol(design), "PASS", "HOLD"),
    ifelse(all(contrasts %in% unique(score$contrast)), "PASS", "HOLD"),
    ifelse(length(all_sets) == 21L, "PASS", "HOLD"),
    ifelse(nrow(new_gsea) == length(all_sets) * 4L * 2L, "PASS", "HOLD"),
    ifelse(all(reuse_check$exact_reuse), "PASS", "HOLD"),
    ifelse(nrow(score) == 252L, "PASS", "HOLD"),
    ifelse(nrow(donor) == 252L && all(donor$n_donors == 3L), "PASS", "HOLD"),
    ifelse(nrow(marker_de) == nrow(marker_manifest) * 12L, "PASS", "HOLD"),
    ifelse(nrow(expr) > 0L, "PASS", "HOLD"),
    ifelse(nrow(marker_donor_rows) == nrow(marker_manifest) * 12L * 3L, "PASS", "HOLD"),
    ifelse(nrow(marker_donor_summary) == nrow(marker_manifest) * 12L, "PASS", "HOLD"),
    ifelse(identical(sort(unique(anno$Treatment)), sort(c("HCl", "LA", "NC", "Unsti"))), "PASS", "HOLD")
  ),
  detail = c(
    paste(sort(unique(anno$Treatment)), collapse = ";"), nrow(factorial),
    paste(qr(design)$rank, "of", ncol(design)), paste(contrasts, collapse = ";"),
    length(all_sets), nrow(new_gsea), sum(reuse_check$exact_reuse), nrow(score),
    nrow(donor), nrow(marker_de), nrow(expr), nrow(marker_donor_rows),
    nrow(marker_donor_summary), paste(sort(unique(anno$Treatment)), collapse = ";")
  )
)
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
overall <- if (all(qa$status == "PASS")) "PASS_CANDIDATE" else "HOLD_CANDIDATE"
writeLines(c(
  paste0("status=", overall),
  paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)),
  paste0("analysis_seed=", analysis_seed),
  paste0("factorial_samples=", nrow(factorial), ";donors=", uniqueN(factorial$Donor)),
  paste0("fixed_sets=", length(all_sets), ";new_control_gsea_rows=", nrow(new_gsea)),
  paste0("marker_genes=", uniqueN(marker_manifest$gene), ";marker_de_rows=", nrow(marker_de))
), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("STATUS\t", overall, "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
