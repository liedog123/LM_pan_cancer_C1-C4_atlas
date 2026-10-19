#!/usr/bin/env Rscript
# Purpose: Fit donor-adjusted CD8 treatment-time responses and fixed programme enrichment.
# Task: cd8_response; step: s03
# Inputs (DCF fields): count_rdata,lm12_membership,msigdb_dir,melanoma_counts,melanoma_deg,bulk_de_dir
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_response --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 03_fixed_programme_gsea.R --config PATH [--check-config]\n       Rscript 03_fixed_programme_gsea.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("count_rdata","lm12_membership","msigdb_dir","melanoma_counts","melanoma_deg","bulk_de_dir","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("count_rdata","lm12_membership","melanoma_counts","melanoma_deg")]
gmt_basenames <- c("h.all.v2025.1.Hs.symbols.gmt","c2.cp.reactome.v2025.1.Hs.symbols.gmt","c5.go.bp.v2024.1.Hs.symbols.gmt","c7.all.v2025.1.Hs.symbols.gmt")
de_contrasts <- c("LA_minus_NC")
de_basenames <- as.vector(outer(de_contrasts, c("1h","6h","24h","48h"), function(x,y) paste0("de_",x,"_",y,".tsv")))
required_files <- c(unname(input_paths), file.path(unname(paths["msigdb_dir"]),gmt_basenames), file.path(unname(paths["bulk_de_dir"]),de_basenames))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- unname(paths[c("msigdb_dir","bulk_de_dir")])
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

# ---- Retained complete source body; documented path adaptations only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(GSVA)
  library(fgsea)
})

# v0.2 专门补足：基因集 ID 映射、donor 一致性/区间、leading-edge 拆解和 GSEA 可复现性。
analysis_seed <- 20260911L
boundary_seed <- 20260912L
set.seed(analysis_seed)

out_dir <- unname(cfg_paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  bulk_input = unname(cfg_paths["count_rdata"]),
  lm12_long = unname(cfg_paths["lm12_membership"]),
  msigdb_root = unname(cfg_paths["msigdb_dir"]),
  melanoma_counts = unname(cfg_paths["melanoma_counts"]),
  melanoma_deg = unname(cfg_paths["melanoma_deg"]),
  bulk_de_root = unname(cfg_paths["bulk_de_dir"])
)

msigdb_files <- c(
  hallmark = file.path(paths$msigdb_root, "h.all.v2025.1.Hs.symbols.gmt"),
  reactome = file.path(paths$msigdb_root, "c2.cp.reactome.v2025.1.Hs.symbols.gmt"),
  go_bp = file.path(paths$msigdb_root, "c5.go.bp.v2024.1.Hs.symbols.gmt"),
  c7 = file.path(paths$msigdb_root, "c7.all.v2025.1.Hs.symbols.gmt")
)

fixed_msigdb <- data.table(
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
  feature_label = c(
    "CD8 activation",
    "TCR signaling",
    "T-cell mediated cytotoxicity",
    "IFN-gamma response",
    "Memory-enriched (Goldrath effector vs memory, DN)",
    "Effector-enriched (Goldrath effector vs memory, UP)",
    "Glycolysis",
    "Hypoxia",
    "Unfolded protein response"
  ),
  source_collection = c("C5_GO_BP", "C2_REACTOME", "C5_GO_BP", "HALLMARK", "C7_IMMUNO", "C7_IMMUNO", "HALLMARK", "HALLMARK", "HALLMARK"),
  source_file_key = c("go_bp", "reactome", "go_bp", "hallmark", "c7", "c7", "hallmark", "hallmark", "hallmark"),
  direction = c("activity", "activity", "activity", "activity", "memory_enriched", "effector_enriched", "activity", "activity", "activity")
)

write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = TRUE, na = "NA")
stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)

for (f in c(unlist(paths[c("bulk_input", "lm12_long", "melanoma_counts", "melanoma_deg", "bulk_de_root")]), unname(msigdb_files))) {
  stop_if(!file.exists(f), paste("missing input:", f))
}

input_paths <- c(unlist(paths[c("bulk_input", "lm12_long", "melanoma_counts", "melanoma_deg")], use.names = TRUE), msigdb_files)
file_manifest <- data.table(
  role = c(c("bulk_input", "lm12_long", "melanoma_counts", "melanoma_deg"), paste0("msigdb_", names(msigdb_files))),
  path = unname(input_paths)
)
file_manifest[, `:=`(exists = file.exists(path), bytes = as.numeric(file.info(path)$size), sha256 = unname(tools::sha256sum(path)))]
write_tsv(file_manifest, file.path(out_dir, "input_file_manifest.tsv"))

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  toks <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(toks, function(x) x[[1L]], character(1))
  sets <- lapply(toks, function(x) unique(x[-seq_len(min(2L, length(x)))]))
  names(sets) <- ids
  sets <- sets[!duplicated(names(sets))]
  list(sets = sets, n_source = lengths(sets))
}

gmt_cache <- lapply(msigdb_files, read_gmt)
names(gmt_cache) <- names(msigdb_files)
msigdb_sets <- lapply(seq_len(nrow(fixed_msigdb)), function(i) {
  key <- fixed_msigdb$source_file_key[i]
  id <- fixed_msigdb$feature_id[i]
  stop_if(!id %in% names(gmt_cache[[key]]$sets), paste("fixed MSigDB set missing:", id, key))
  gmt_cache[[key]]$sets[[id]]
})
names(msigdb_sets) <- fixed_msigdb$feature_id

lm12_tab <- fread(paths$lm12_long, sep = "\t", quote = "\"", data.table = TRUE)
stop_if(!all(c("feature_id", "gene", "pathway") %in% names(lm12_tab)), "LM12 membership columns incomplete")
lm12_ids <- unique(as.character(lm12_tab$feature_id))
lm12_sets <- setNames(lapply(lm12_ids, function(id) unique(as.character(lm12_tab[feature_id == id, gene]))), lm12_ids)
lm12_meta <- data.table(feature_id = lm12_ids, feature_label = vapply(lm12_ids, function(id) unique(as.character(lm12_tab[feature_id == id, pathway]))[1L], character(1)))
lm12_meta[, `:=`(feature_layer = "LM12_FIG2_FROZEN", source_collection = "FIG2_FROZEN_LM12", source_file_key = "lm12_long", direction = "activity")]

msigdb_meta <- copy(fixed_msigdb)
msigdb_meta[, `:=`(feature_layer = "MSIGDB_FULL_SET", source_path = unname(msigdb_files[source_file_key]))]
msigdb_meta[, n_genes_source := lengths(msigdb_sets[feature_id])]
msigdb_meta[, set_sha256 := vapply(msigdb_sets[feature_id], function(g) digest::digest(sort(unique(g)), algo = "sha256"), character(1))]

all_sets <- c(setNames(lm12_sets, paste0("LM12::", names(lm12_sets))), setNames(msigdb_sets, paste0("MSIGDB::", names(msigdb_sets))))
feature_meta <- rbind(
  lm12_meta[, .(feature_id = paste0("LM12::", feature_id), feature_label, feature_layer, source_collection, direction)],
  msigdb_meta[, .(feature_id = paste0("MSIGDB::", feature_id), feature_label, feature_layer, source_collection, direction)]
)
feature_meta[, feature_id_raw := sub("^[^:]+::", "", feature_id)]
set_label <- setNames(feature_meta$feature_label, feature_meta$feature_id)
set_layer <- setNames(feature_meta$feature_layer, feature_meta$feature_id)
source_size <- data.table(feature_id = names(all_sets), n_genes_source = lengths(all_sets))
source_path <- c(setNames(rep(paths$lm12_long, length(lm12_sets)), names(all_sets)[seq_along(lm12_sets)]),
                 setNames(unname(msigdb_files[msigdb_meta$source_file_key]), names(all_sets)[(length(lm12_sets) + 1L):length(all_sets)]))

write_tsv(rbind(
  lm12_meta[, .(feature_id, feature_label, feature_layer, source_collection, source_path = paths$lm12_long, direction, n_genes_source = unname(lengths(lm12_sets)[feature_id]), set_sha256 = vapply(lm12_sets[feature_id], function(g) digest::digest(sort(unique(g)), algo = "sha256"), character(1)))],
  msigdb_meta[, .(feature_id, feature_label, feature_layer, source_collection, source_path, direction, n_genes_source, set_sha256)]
), file.path(out_dir, "fixed_gene_set_manifest.tsv"))

load_count <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  stop_if(!all(c("count", "sample_anno") %in% ls(e, all.names = TRUE)), "bulk RData lacks count/sample_anno")
  list(count = get("count", e), sample_anno = get("sample_anno", e))
}

score_ssgsea <- function(logcpm, sets, dataset) {
  present <- lapply(sets, function(g) intersect(unique(g), rownames(logcpm)))
  names(present) <- names(sets)
  status <- data.table(dataset = dataset, feature_id = names(present), n_genes_input = lengths(sets), n_genes_present = lengths(present), coverage = lengths(present) / lengths(sets), status = ifelse(lengths(present) >= 3L, "PASS", "HOLD"))
  stop_if(any(status$status == "HOLD"), paste("set coverage below 3 genes:", paste(status[status == "HOLD", feature_id], collapse = ",")))
  set.seed(analysis_seed)
  score <- GSVA::gsva(GSVA::ssgseaParam(exprData = logcpm, geneSets = present, minSize = 1L, alpha = 0.25, normalize = TRUE), verbose = FALSE)
  score <- as.matrix(score)
  returned_ids <- rownames(score)
  stop_if(is.null(returned_ids) || !setequal(returned_ids, names(present)), paste(dataset, "GSVA returned IDs do not match fixed gene-set IDs"))
  returned_index <- match(names(present), returned_ids)
  stop_if(anyNA(returned_index), paste(dataset, "GSVA ID matching failed"))
  score <- score[returned_index, , drop = FALSE]
  rownames(score) <- names(present)
  mapping <- data.table(dataset = dataset, expected_id = names(present), returned_id = returned_ids[returned_index], returned_index = returned_index, mapping_ok = names(present) == returned_ids[returned_index])
  stop_if(!all(mapping$mapping_ok), paste(dataset, "GSVA ID mapping is not exact"))
  list(score = score, present = present, coverage = status, mapping = mapping)
}

make_long <- function(score, sample_sheet, score_layer) {
  long <- melt(as.data.table(t(score), keep.rownames = "sample_id"), id.vars = "sample_id", variable.name = "feature_id", value.name = "ssgsea_score")
  long[, feature_id := as.character(feature_id)]
  long <- merge(long, as.data.table(sample_sheet), by = "sample_id", all.x = TRUE, sort = FALSE)
  long[, `:=`(score_layer = score_layer, feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]))]
  stop_if(anyNA(long$feature_label) || anyNA(long$feature_layer), paste(score_layer, "feature label mapping failed"))
  setcolorder(long, c("sample_id", "feature_id", "feature_label", "feature_layer", "score_layer"))
  long[]
}

model_score_contrasts <- function(score, design, contrast_specs, feature_meta_local) {
  fit <- limma::lmFit(score, design)
  coef_names <- colnames(design)
  find_coef <- function(term) {
    hit <- which(coef_names == term)
    stop_if(length(hit) != 1L, paste("coefficient not unique:", term))
    hit
  }
  out <- lapply(seq_len(nrow(contrast_specs)), function(i) {
    cf <- limma::contrasts.fit(fit, contrast_specs$vector[[i]])
    cf <- limma::eBayes(cf, robust = TRUE)
    tt <- as.data.table(limma::topTable(cf, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature_id")
    stop_if(!setequal(tt$feature_id, rownames(score)), "score contrast lost feature IDs")
    tt <- tt[match(rownames(score), feature_id)]
    se <- as.numeric(cf$stdev.unscaled[, 1L] * cf$sigma)
    names(se) <- rownames(cf$stdev.unscaled)
    se <- se[match(tt$feature_id, names(se))]
    df_resid <- max(cf$df.residual, na.rm = TRUE)
    tt[, `:=`(contrast = contrast_specs$contrast[i], time = contrast_specs$time[i], SE = se, df = df_resid)]
    tt[, `:=`(CI_low = logFC - qt(0.975, df_resid) * SE, CI_high = logFC + qt(0.975, df_resid) * SE)]
    tt[, .(feature_id, contrast, time, logFC, AveExpr, t, P.Value, adj.P.Val, B, SE, df, CI_low, CI_high)]
  })
  merge(rbindlist(out, fill = TRUE), feature_meta_local, by = "feature_id", all.x = TRUE, sort = FALSE)
}

bulk_de_stats <- function(time) {
  f <- file.path(paths$bulk_de_root, paste0("de_LA_minus_NC_", time, ".tsv"))
  stop_if(!file.exists(f), paste("missing existing bulk DE:", f))
  d <- fread(f, sep = "\t")
  d[, rank_score := sign(logFC) * -log10(pmax(PValue, .Machine$double.xmin))]
  d[!duplicated(gene_id), .(gene = gene_id, rank_score, logFC, PValue, FDR)]
}

run_fgsea <- function(stats, sets, context, time = NA_character_, donor = NA_character_, seed = analysis_seed) {
  stats <- stats[is.finite(stats)]
  stats <- stats[!duplicated(names(stats))]
  set.seed(seed)
  res <- fgsea::fgseaMultilevel(pathways = sets, stats = stats, minSize = 3L, maxSize = 5000L, eps = 0)
  res <- as.data.table(res)
  res[, pathway := as.character(pathway)]
  res[, `:=`(context = context, time = time, donor = donor, feature_id = pathway, seed = seed)]
  res[, pathway := NULL]
  res[]
}

build_bulk <- function() {
  obj <- load_count(paths$bulk_input)
  count <- obj$count
  ann <- as.data.table(obj$sample_anno, keep.rownames = "sample_id")
  stop_if(!identical(colnames(count), ann$sample_id), "bulk count/sample annotation order mismatch")
  ann <- ann[Treatment %in% c("NC", "HCl", "LA")]
  ann[, `:=`(Donor = factor(Donor), Treatment = factor(Treatment, levels = c("NC", "HCl", "LA")), Time = factor(Time, levels = c("1h", "6h", "24h", "48h")))]
  stop_if(anyNA(ann$Donor) || anyNA(ann$Treatment) || anyNA(ann$Time), "bulk factorial metadata contains NA")
  design <- model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = ann)
  stop_if(qr(design)$rank != ncol(design), "bulk design not full rank")
  dge <- edgeR::DGEList(counts = count[, ann$sample_id, drop = FALSE], samples = ann)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  scored <- score_ssgsea(logcpm, all_sets, "bulk")
  sample_sheet <- as.data.frame(ann)
  score_long <- make_long(scored$score, sample_sheet, "bulk_ssGSEA")
  coef_names <- colnames(design)
  find_coef <- function(term) {
    hit <- which(coef_names == term)
    stop_if(length(hit) != 1L, paste("coefficient not unique:", term))
    hit
  }
  contrast_vector <- function(left, right, time) {
    v <- numeric(ncol(design))
    add <- function(treatment, sign) {
      if (treatment != "NC") {
        v[find_coef(paste0("Treatment", treatment))] <<- v[find_coef(paste0("Treatment", treatment))] + sign
        if (time != "1h") {
          hit <- which(coef_names %in% c(paste0("Time", time, ":Treatment", treatment), paste0("Treatment", treatment, ":Time", time)))
          stop_if(length(hit) != 1L, "bulk interaction coefficient missing")
          v[hit] <<- v[hit] + sign
        }
      }
    }
    add(left, 1); add(right, -1); v
  }
  specs <- CJ(contrast = c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"), time = c("1h", "6h", "24h", "48h"), unique = TRUE)
  specs[, `:=`(left = c("LA", "LA", "HCl")[match(contrast, c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))], right = c("NC", "HCl", "NC")[match(contrast, c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))])]
  specs[, vector := lapply(seq_len(.N), function(i) contrast_vector(left[i], right[i], time[i]))]
  contrasts <- model_score_contrasts(scored$score, design, specs, feature_meta)
  donor_rows <- list()
  for (i in seq_len(nrow(specs))) {
    sp <- specs[i]
    sub <- score_long[Time == sp$time & Treatment %in% c(sp$left, sp$right)]
    left <- sub[Treatment == sp$left, .(Donor = as.character(Donor), feature_id, left_score = ssgsea_score)]
    right <- sub[Treatment == sp$right, .(Donor = as.character(Donor), feature_id, right_score = ssgsea_score)]
    dd <- merge(left, right, by = c("Donor", "feature_id"), all = FALSE)
    dd[, `:=`(contrast = sp$contrast, time = sp$time, donor_delta = left_score - right_score)]
    donor_rows[[length(donor_rows) + 1L]] <- dd[, .(contrast, time, feature_id, donor = Donor, donor_delta)]
  }
  donor_delta <- rbindlist(donor_rows, fill = TRUE)
  donor_summary <- donor_delta[, .(n_donors = uniqueN(donor), mean_delta = mean(donor_delta), median_delta = median(donor_delta), sd_delta = sd(donor_delta), n_positive = sum(donor_delta > 0), n_negative = sum(donor_delta < 0), n_zero = sum(donor_delta == 0)), by = .(contrast, time, feature_id)]
  donor_summary[, se_delta := sd_delta / sqrt(n_donors)]
  donor_summary[, `:=`(CI_low = mean_delta - qt(0.975, pmax(n_donors - 1, 1)) * se_delta, CI_high = mean_delta + qt(0.975, pmax(n_donors - 1, 1)) * se_delta, donor_all_same_direction = n_positive == n_donors | n_negative == n_donors)]
  list(score_long = score_long, score_contrasts = contrasts, donor_delta = donor_delta, donor_summary = donor_summary, sample_sheet = sample_sheet, ann = ann, logcpm = logcpm, present = scored$present, coverage = scored$coverage, mapping = scored$mapping, design = design, design_rank = qr(design)$rank)
}

build_melanoma <- function() {
  obj <- readRDS(paths$melanoma_counts)
  stop_if(!all(c("counts", "sample_sheet") %in% names(obj)), "melanoma CD8 pseudobulk bundle incomplete")
  count <- obj$counts; ss_all <- as.data.table(obj$sample_sheet)
  req <- c("unit_id", "sample_ID_long", "patient", "sample_C_subtype", "pre_post", "eligible_ge20", "n_cells", "library_size")
  stop_if(!all(req %in% names(ss_all)), "melanoma sample sheet incomplete")
  ss_all <- ss_all[eligible_ge20 == TRUE & sample_C_subtype %in% c("C1", "C4")]
  stop_if(anyDuplicated(ss_all$unit_id) > 0L, "duplicate melanoma CD8 unit IDs")
  ss <- copy(ss_all); selection_mode <- "all_patients_unique"
  if (anyDuplicated(ss$patient) > 0L) {
    ss[, pre_priority := fifelse(pre_post == "Pre", 0L, 1L)]
    setorder(ss, patient, pre_priority, -n_cells, -library_size, sample_ID_long)
    ss <- ss[, .SD[1L], by = patient]; ss[, pre_priority := NULL]
    selection_mode <- "independent_patient_fallback_pre_priority"
  }
  stop_if(anyDuplicated(ss$patient) > 0L, "duplicate melanoma patients after selection")
  selection_manifest <- copy(ss_all); selection_manifest[, retained := unit_id %in% ss$unit_id]; selection_manifest[, selection_reason := fifelse(retained, selection_mode, "dropped_duplicate_patient")]
  count <- count[, ss$unit_id, drop = FALSE]
  group <- factor(ss$sample_C_subtype, levels = c("C4", "C1")); ss[, group := group]; ss[, pre_post := factor(pre_post)]
  design <- model.matrix(~ 0 + group + pre_post, data = ss)
  stop_if(qr(design)$rank != ncol(design), "melanoma CD8 design not full rank")
  dge <- edgeR::DGEList(counts = count, samples = ss); dge <- edgeR::calcNormFactors(dge, method = "TMM"); logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  scored <- score_ssgsea(logcpm, all_sets, "melanoma_CD8")
  sample_sheet <- data.frame(sample_id = ss$unit_id, sample_ID_long = ss$sample_ID_long, patient = ss$patient, state = ss$sample_C_subtype, pre_post = as.character(ss$pre_post), cell_type = "CD8 T", n_cells = ss$n_cells, stringsAsFactors = FALSE)
  score_long <- make_long(scored$score, sample_sheet, "melanoma_CD8_ssGSEA")
  fit <- limma::lmFit(scored$score, design); cf <- limma::contrasts.fit(fit, c(-1, 1, rep(0, ncol(design) - 2L))); cf <- limma::eBayes(cf, robust = TRUE)
  tt <- as.data.table(limma::topTable(cf, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature_id")
  stop_if(!setequal(tt$feature_id, rownames(scored$score)), "melanoma score contrast lost feature IDs")
  tt <- tt[match(rownames(scored$score), feature_id)]
  se <- as.numeric(cf$stdev.unscaled[, 1L] * cf$sigma); names(se) <- rownames(cf$stdev.unscaled); se <- se[match(tt$feature_id, names(se))]
  df_resid <- max(cf$df.residual, na.rm = TRUE)
  tt[, `:=`(contrast = "C1_minus_C4", SE = se, df = df_resid)]
  tt[, `:=`(CI_low = logFC - qt(0.975, df_resid) * SE, CI_high = logFC + qt(0.975, df_resid) * SE)]
  tt <- tt[, .(feature_id, contrast, logFC, AveExpr, t, P.Value, adj.P.Val, B, SE, df, CI_low, CI_high)]
  tt <- merge(tt, feature_meta, by = "feature_id", all.x = TRUE, sort = FALSE)
  group_summary <- score_long[, .(mean_score = mean(ssgsea_score), sd_score = sd(ssgsea_score), n_samples = .N, n_patients = uniqueN(patient)), by = .(feature_id, feature_label, feature_layer, state)]
  deg <- fread(cmd = paste("gzip -cd", shQuote(paths$melanoma_deg)), sep = "\t", quote = "\"", data.table = TRUE); deg <- deg[cell_type == "CD8 T" & comparison == "C1_vs_C4"]
  union_genes <- unique(unlist(all_sets)); deg_fixed <- deg[gene %in% union_genes]
  list(score_long = score_long, score_contrasts = tt, group_summary = group_summary, sample_sheet = sample_sheet, present = scored$present, coverage = scored$coverage, mapping = scored$mapping, design = design, design_rank = qr(design)$rank, n_samples = nrow(ss), n_patients = uniqueN(ss$patient), deg_fixed = deg_fixed, selection_manifest = selection_manifest, selection_mode = selection_mode, logcpm = logcpm)
}

if (FALSE) {
# v0.1 tail retained below only as an auditable historical code path; it is not executed.
bulk <- build_bulk(); mel <- build_melanoma()
membership <- merge(feature_meta, merge(bulk$coverage[, .(feature_id, n_genes_bulk = n_genes_present, bulk_coverage = coverage)], mel$coverage[, .(feature_id, n_genes_melanoma = n_genes_present, melanoma_coverage = coverage)], by = "feature_id", all = TRUE), by = "feature_id", all = TRUE, sort = FALSE)
membership[, `:=`(n_genes_source = ifelse(grepl("^LM12::", feature_id), lengths(lm12_sets)[sub("^LM12::", "", feature_id)], msigdb_meta$n_genes_source[match(sub("^MSIGDB::", "", feature_id), msigdb_meta$feature_id)]), source_path = ifelse(grepl("^LM12::", feature_id), paths$lm12_long, unname(msigdb_files[msigdb_meta$source_file_key[match(sub("^MSIGDB::", "", feature_id), msigdb_meta$feature_id)] ]))) ]
membership[, `:=`(bulk_score_status = ifelse(n_genes_bulk >= 3L, "PASS", "HOLD"), melanoma_score_status = ifelse(n_genes_melanoma >= 3L, "PASS", "HOLD"))]

write_tsv(membership, file.path(out_dir, "fixed_gene_set_coverage_manifest.tsv"))
write_tsv(rbind(bulk$mapping, mel$mapping, fill = TRUE), file.path(out_dir, "score_id_mapping_qa.tsv"))
write_tsv(bulk$score_long, file.path(out_dir, "bulk_fixed_readout_ssgsea_long.tsv"))
write_tsv(bulk$score_contrasts, file.path(out_dir, "bulk_fixed_readout_ssgsea_contrasts.tsv"))
write_tsv(bulk$donor_delta, file.path(out_dir, "bulk_fixed_readout_donor_deltas.tsv"))
write_tsv(bulk$donor_summary, file.path(out_dir, "bulk_fixed_readout_donor_summary.tsv"))
write_tsv(mel$score_long, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_long.tsv"))
write_tsv(mel$score_contrasts, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_c1_vs_c4.tsv"))
write_tsv(mel$group_summary, file.path(out_dir, "melanoma_cd8_fixed_readout_group_summary.tsv"))
write_tsv(mel$selection_manifest, file.path(out_dir, "melanoma_cd8_patient_selection_manifest.tsv"))
write_tsv(mel$deg_fixed, file.path(out_dir, "melanoma_cd8_fixed_readout_gene_deg.tsv"))
write_tsv(bridge, file.path(out_dir, "melanoma_cd8_bulk_bridge_by_time.tsv"))
write_tsv(bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_bridge_summary.tsv"))
write_tsv(gene_bridge, file.path(out_dir, "melanoma_cd8_bulk_gene_bridge.tsv"))
write_tsv(flatten_gsea(bulk_fgsea_raw), file.path(out_dir, "bulk_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(mel_fgsea_raw), file.path(out_dir, "melanoma_cd8_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(bulk_donor_gsea_raw), file.path(out_dir, "bulk_donor_directional_preranked_gsea.tsv"))
write_tsv(gsea_bridge, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_by_time.tsv"))
write_tsv(gsea_bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_summary.tsv"))
write_tsv(gsea_repro, file.path(out_dir, "gsea_reproducibility_boundary.tsv"))
write_tsv(method_direction_audit, file.path(out_dir, "ssgsea_gsea_direction_audit.tsv"))
write_tsv(memory_direction_audit, file.path(out_dir, "memory_direction_audit.tsv"))
write_tsv(diag_specs, file.path(out_dir, "diagnostic_gene_set_manifest.tsv"))
write_tsv(edge_long, file.path(out_dir, "gsea_leading_edge_long.tsv"))
write_tsv(edge_cat, file.path(out_dir, "gsea_leading_edge_category_long.tsv"))
write_tsv(edge_summary, file.path(out_dir, "gsea_leading_edge_category_summary.tsv"))
write_tsv(gsea_repro[boundary_a | boundary_b | boundary_changed], file.path(out_dir, "gsea_boundary_focus.tsv"))
write_tsv(data.table(
  object = c("bulk_count", "melanoma_cd8_pseudobulk"),
  analysis_unit = c("donor x condition x time sample", "sample x native CD8 T pseudobulk; patient biological unit"),
  primary_comparison = c("time-matched LA_minus_NC", "C1_minus_C4"),
  n_independent_units = c(length(unique(bulk$sample_sheet$Donor)), mel$n_patients),
  interpretation = c("CD8 bulk aggregate response", "melanoma CD8 state association")
), file.path(out_dir, "analysis_unit_manifest.tsv"))

mapping_ok <- all(bulk$mapping$mapping_ok) && all(mel$mapping$mapping_ok)
gsea_rows_ok <- nrow(bulk_fgsea_raw) == length(all_sets) * 4L && nrow(mel_fgsea_raw) == length(all_sets)
donor_rows_ok <- nrow(bulk$donor_delta[contrast == "LA_minus_NC"]) == length(all_sets) * 4L * 3L
repro_ok <- nrow(gsea_repro) > 0L && all(is.finite(gsea_repro$abs_NES_delta)) && all(is.finite(gsea_repro$abs_FDR_delta))
qa <- data.table(
  check_id = c("bulk_factorial_samples", "bulk_design_rank", "melanoma_c1c4_samples", "melanoma_c1c4_patients", "fixed_readout_count", "score_id_mapping_exact", "bulk_coverage", "melanoma_coverage", "all_ssgsea_finite", "bulk_primary_donor_rows", "gsea_rows_all_fixed", "donor_gsea_rows", "gsea_seed_reproducible_run", "leading_edge_table"),
  status = c(
    ifelse(nrow(bulk$sample_sheet) == 36L, "PASS", "HOLD"),
    ifelse(bulk$design_rank == ncol(bulk$design), "PASS", "HOLD"),
    ifelse(mel$n_samples >= 6L, "PASS", "HOLD"),
    ifelse(mel$n_patients >= 6L, "PASS", "HOLD"),
    ifelse(length(all_sets) == 21L, "PASS", "HOLD"),
    ifelse(mapping_ok, "PASS", "HOLD"),
    ifelse(all(membership$bulk_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(membership$melanoma_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(is.finite(bulk$score_long$ssgsea_score)) && all(is.finite(mel$score_long$ssgsea_score)), "PASS", "HOLD"),
    ifelse(donor_rows_ok, "PASS", "HOLD"),
    ifelse(gsea_rows_ok, "PASS", "HOLD"),
    ifelse(nrow(bulk_donor_gsea_raw) == length(all_sets) * 4L * 3L, "PASS", "HOLD"),
    ifelse(repro_ok, "PASS", "HOLD"),
    ifelse(nrow(edge_long) > 0L, "PASS", "HOLD")
  ),
  detail = c(nrow(bulk$sample_sheet), paste(bulk$design_rank, "of", ncol(bulk$design)), mel$n_samples, mel$n_patients, length(all_sets), paste(all(bulk$mapping$mapping_ok), all(mel$mapping$mapping_ok), sep = ";"), sprintf("min %.3f", min(membership$bulk_coverage)), sprintf("min %.3f", min(membership$melanoma_coverage)), paste(nrow(bulk$score_long), nrow(mel$score_long), sep = " + "), nrow(bulk$donor_delta[contrast == "LA_minus_NC"]), paste(nrow(bulk_fgsea_raw), nrow(mel_fgsea_raw), sep = " + "), nrow(bulk_donor_gsea_raw), sprintf("max |delta NES| %.4g; max |delta FDR| %.4g", max(gsea_repro$abs_NES_delta), max(gsea_repro$abs_FDR_delta)), nrow(edge_long))
)
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
overall_status <- if (all(qa$status == "PASS")) "PASS_CANDIDATE" else "HOLD_CANDIDATE"
writeLines(c(paste0("status=", overall_status), paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)), paste0("analysis_seed=", analysis_seed), paste0("boundary_seed=", boundary_seed), paste0("fixed_readouts=", length(all_sets)), paste0("bulk_samples=", nrow(bulk$sample_sheet), ";bulk_donors=", length(unique(bulk$sample_sheet$Donor))), paste0("melanoma_cd8_samples=", mel$n_samples, ";melanoma_patients=", mel$n_patients, ";selection=", mel$selection_mode), paste0("gsea_rows=", nrow(bulk_fgsea_raw), "+", nrow(mel_fgsea_raw)), paste0("leading_edge_rows=", nrow(edge_long))), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("STATUS\t", overall_status, "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")

bulk_primary <- bulk$score_contrasts[contrast == "LA_minus_NC"]
mel_primary <- mel$score_contrasts
bridge <- merge(mel_primary[, .(feature_id, melanoma_logFC = logFC, melanoma_p = P.Value, melanoma_fdr = adj.P.Val)], bulk_primary[, .(feature_id, time, bulk_logFC = logFC, bulk_p = P.Value, bulk_fdr = adj.P.Val)], by = "feature_id", allow.cartesian = TRUE)
donor_support <- merge(bulk$donor_delta[contrast == "LA_minus_NC"], mel_primary[, .(feature_id, melanoma_logFC = logFC)], by = "feature_id", all.x = TRUE)
donor_support[, donor_same_direction := sign(donor_delta) == sign(melanoma_logFC)]
support_summary <- donor_support[, .(n_donors = uniqueN(donor), n_donors_same_direction = sum(donor_same_direction, na.rm = TRUE), donor_same_fraction = mean(donor_same_direction, na.rm = TRUE)), by = .(feature_id, time)]
bridge <- merge(bridge, support_summary, by = c("feature_id", "time"), all.x = TRUE)
bridge[, `:=`(feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]), melanoma_direction = fifelse(melanoma_logFC > 0, "positive", fifelse(melanoma_logFC < 0, "negative", "zero")), bulk_direction = fifelse(bulk_logFC > 0, "positive", fifelse(bulk_logFC < 0, "negative", "zero")))]
bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
bridge[, evidence_class := ifelse(melanoma_fdr < 0.05 & bulk_fdr < 0.05 & donor_same_fraction >= 2/3, ifelse(direction_relation == "same_direction", "supported_concordance", "supported_divergence"), ifelse(direction_relation == "same_direction" & (donor_same_fraction >= 1 | melanoma_p < 0.05 | bulk_p < 0.05), "suggestive_concordance", "insufficient_or_mixed"))]
bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(bridge, layer_order, feature_id, time)
bridge[, layer_order := NULL]
bridge_summary <- bridge[, .(bridge_class = if (all(evidence_class == "supported_concordance")) "supported_concordance" else if (any(evidence_class == "supported_divergence")) "supported_divergence_or_mixed" else if (any(evidence_class == "suggestive_concordance")) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), mean_abs_bulk_logFC = mean(abs(bulk_logFC), na.rm = TRUE), melanoma_logFC = first(melanoma_logFC), melanoma_fdr = first(melanoma_fdr)), by = .(feature_id, feature_label, feature_layer)]

gene_bridge <- merge(mel$deg_fixed[, .(gene, melanoma_logFC = log2fc_c1_minus_ref, melanoma_p = p_value, melanoma_fdr = fdr)], rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) { d <- bulk_de_stats(tm); d[, time := tm]; d[gene %in% unique(unlist(all_sets))] })), by = "gene", allow.cartesian = TRUE)
gene_bridge[, direction_relation := fifelse(sign(melanoma_logFC) == sign(logFC), "same_direction", "opposite_direction")]
gene_bridge[, set_membership := vapply(gene, function(g) paste(names(all_sets)[vapply(all_sets, function(s) g %in% s, logical(1))], collapse = ";"), character(1))]

bulk_fgsea <- rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) {
  d <- bulk_de_stats(tm); stats <- d$rank_score; names(stats) <- d$gene
  run_fgsea(stats, msigdb_sets, "bulk_LA_minus_NC", tm)
}), fill = TRUE)
mel_deg <- mel$deg_fixed
mel_all_deg <- fread(cmd = paste("gzip -cd", shQuote(paths$melanoma_deg)), sep = "\t", quote = "\"", data.table = TRUE)[cell_type == "CD8 T" & comparison == "C1_vs_C4"]
mel_stats <- mel_all_deg$moderated_t; names(mel_stats) <- mel_all_deg$gene
mel_fgsea <- run_fgsea(mel_stats, msigdb_sets, "melanoma_CD8_C1_minus_C4")
gsea_bridge <- merge(
  mel_fgsea[, .(feature_id, melanoma_NES = NES, melanoma_p = pval, melanoma_fdr = padj, melanoma_size = size)],
  bulk_fgsea[, .(feature_id, time, bulk_NES = NES, bulk_p = pval, bulk_fdr = padj, bulk_size = size)],
  by = "feature_id", allow.cartesian = TRUE
)
gsea_bridge[, `:=`(
  feature_label = unname(set_label[feature_id]),
  feature_layer = unname(set_layer[feature_id]),
  melanoma_direction = fifelse(melanoma_NES > 0, "positive", fifelse(melanoma_NES < 0, "negative", "zero")),
  bulk_direction = fifelse(bulk_NES > 0, "positive", fifelse(bulk_NES < 0, "negative", "zero"))
)]
gsea_bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
gsea_bridge[, evidence_class := ifelse(melanoma_fdr < 0.05 & bulk_fdr < 0.05, ifelse(direction_relation == "same_direction", "supported_concordance", "supported_divergence"), ifelse(direction_relation == "same_direction" & (melanoma_p < 0.05 | bulk_p < 0.05), "suggestive_concordance", "insufficient_or_mixed"))]
gsea_bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(gsea_bridge, layer_order, feature_id, time)
gsea_bridge[, layer_order := NULL]
gsea_bridge_summary <- gsea_bridge[, .(bridge_class = if (all(evidence_class == "supported_concordance")) "supported_concordance" else if (any(evidence_class == "supported_divergence")) "supported_divergence_or_mixed" else if (any(evidence_class == "suggestive_concordance")) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), melanoma_NES = first(melanoma_NES), melanoma_fdr = first(melanoma_fdr), mean_bulk_abs_NES = mean(abs(bulk_NES), na.rm = TRUE)), by = .(feature_id, feature_label, feature_layer)]

write_tsv(membership, file.path(out_dir, "fixed_gene_set_coverage_manifest.tsv"))
write_tsv(bulk$score_long, file.path(out_dir, "bulk_fixed_readout_ssgsea_long.tsv"))
write_tsv(bulk$score_contrasts, file.path(out_dir, "bulk_fixed_readout_ssgsea_contrasts.tsv"))
write_tsv(bulk$donor_delta, file.path(out_dir, "bulk_fixed_readout_donor_deltas.tsv"))
write_tsv(mel$score_long, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_long.tsv"))
write_tsv(mel$score_contrasts, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_c1_vs_c4.tsv"))
write_tsv(mel$group_summary, file.path(out_dir, "melanoma_cd8_fixed_readout_group_summary.tsv"))
write_tsv(mel$selection_manifest, file.path(out_dir, "melanoma_cd8_patient_selection_manifest.tsv"))
write_tsv(mel$deg_fixed, file.path(out_dir, "melanoma_cd8_fixed_readout_gene_deg.tsv"))
write_tsv(bridge, file.path(out_dir, "melanoma_cd8_bulk_bridge_by_time.tsv"))
write_tsv(bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_bridge_summary.tsv"))
write_tsv(gene_bridge, file.path(out_dir, "melanoma_cd8_bulk_gene_bridge.tsv"))
write_tsv(bulk_fgsea, file.path(out_dir, "bulk_msigdb_preranked_gsea.tsv"))
write_tsv(mel_fgsea, file.path(out_dir, "melanoma_cd8_msigdb_preranked_gsea.tsv"))
write_tsv(gsea_bridge, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_by_time.tsv"))
write_tsv(gsea_bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_summary.tsv"))
write_tsv(rbind(bulk$coverage[, `:=`(dataset = "bulk")], mel$coverage[, `:=`(dataset = "melanoma_CD8")]), file.path(out_dir, "score_coverage_qa.tsv"))
write_tsv(data.table(
  check_id = c("bulk_factorial_samples", "bulk_design_rank", "melanoma_c1c4_samples", "melanoma_c1c4_patients", "melanoma_unique_patient_units", "fixed_readout_count", "bulk_coverage", "melanoma_coverage", "all_ssgsea_finite", "bridge_rows", "bulk_fgsea_rows", "melanoma_fgsea_rows"),
  status = c(ifelse(nrow(bulk$sample_sheet) == 36L, "PASS", "HOLD"), ifelse(bulk$design_rank == ncol(bulk$design), "PASS", "HOLD"), ifelse(mel$n_samples >= 6L, "PASS", "HOLD"), ifelse(mel$n_patients >= 6L, "PASS", "HOLD"), ifelse(!anyDuplicated(mel$sample_sheet$patient), "PASS", "HOLD"), ifelse(length(all_sets) == 21L, "PASS", "HOLD"), ifelse(all(membership$bulk_score_status == "PASS"), "PASS", "HOLD"), ifelse(all(membership$melanoma_score_status == "PASS"), "PASS", "HOLD"), ifelse(all(is.finite(bulk$score_long$ssgsea_score)) && all(is.finite(mel$score_long$ssgsea_score)), "PASS", "HOLD"), ifelse(nrow(bridge) == length(all_sets) * 4L, "PASS", "HOLD"), ifelse(nrow(bulk_fgsea) == length(msigdb_sets) * 4L, "PASS", "HOLD"), ifelse(nrow(mel_fgsea) == length(msigdb_sets), "PASS", "HOLD")),
  detail = c(nrow(bulk$sample_sheet), paste(bulk$design_rank, "of", ncol(bulk$design)), mel$n_samples, mel$n_patients, uniqueN(mel$sample_sheet$patient), length(all_sets), sprintf("min %.3f", min(membership$bulk_coverage)), sprintf("min %.3f", min(membership$melanoma_coverage)), paste(nrow(bulk$score_long), nrow(mel$score_long), sep = " + "), nrow(bridge), nrow(bulk_fgsea), nrow(mel_fgsea))
), file.path(out_dir, "qa_summary.tsv"))
write_tsv(data.table(object = c("bulk_count", "melanoma_cd8_pseudobulk"), analysis_unit = c("donor x condition x time sample", "sample x native CD8 T pseudobulk; patient biological unit"), primary_comparison = c("time-matched LA_minus_NC", "C1_minus_C4"), n_independent_units = c(length(unique(bulk$sample_sheet$Donor)), mel$n_patients), interpretation = c("CD8 bulk aggregate response", "melanoma CD8 state association")), file.path(out_dir, "analysis_unit_manifest.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
writeLines(c("status=PASS_CANDIDATE", paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)), paste0("fixed_readouts=", length(all_sets)), paste0("bulk_samples=", nrow(bulk$sample_sheet), ";bulk_donors=", length(unique(bulk$sample_sheet$Donor))), paste0("melanoma_cd8_samples=", mel$n_samples, ";melanoma_patients=", mel$n_patients, ";selection=", mel$selection_mode), paste0("bridge_rows=", nrow(bridge))), file.path(out_dir, "run_summary.txt"))
cat("STATUS\tPASS_CANDIDATE\n")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
}

# ----------------------------- v0.2 execution -----------------------------
bulk <- build_bulk()
mel <- build_melanoma()

# donor-specific rank enrichment is descriptive and supplies a direction gate;
# it does not replace the formal limma rank used for the bulk GSEA.
build_bulk_donor_gsea <- function(bulk_obj) {
  out <- list(); k <- 0L
  for (donor in unique(as.character(bulk_obj$ann$Donor))) for (tm in c("1h", "6h", "24h", "48h")) {
    la <- bulk_obj$ann$sample_id[as.character(bulk_obj$ann$Donor) == donor & as.character(bulk_obj$ann$Time) == tm & as.character(bulk_obj$ann$Treatment) == "LA"]
    nc <- bulk_obj$ann$sample_id[as.character(bulk_obj$ann$Donor) == donor & as.character(bulk_obj$ann$Time) == tm & as.character(bulk_obj$ann$Treatment) == "NC"]
    stop_if(length(la) != 1L || length(nc) != 1L, paste("donor rank pair missing", donor, tm))
    stats <- bulk_obj$logcpm[, la] - bulk_obj$logcpm[, nc]
    names(stats) <- rownames(bulk_obj$logcpm)
    k <- k + 1L
    out[[k]] <- run_fgsea(stats, all_sets, "bulk_donor_LA_minus_NC", tm, donor, analysis_seed)
  }
  rbindlist(out, fill = TRUE)
}

bulk_donor_gsea_raw <- build_bulk_donor_gsea(bulk)
bulk_fgsea_raw <- rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) {
  d <- bulk_de_stats(tm); stats <- d$rank_score; names(stats) <- d$gene
  run_fgsea(stats, all_sets, "bulk_LA_minus_NC", tm, seed = analysis_seed)
}), fill = TRUE)
mel_all_deg <- fread(cmd = paste("gzip -cd", shQuote(paths$melanoma_deg)), sep = "\t", quote = "\"", data.table = TRUE)
mel_all_deg <- mel_all_deg[cell_type == "CD8 T" & comparison == "C1_vs_C4"]
mel_stats <- mel_all_deg$moderated_t; names(mel_stats) <- mel_all_deg$gene
mel_fgsea_raw <- run_fgsea(mel_stats, all_sets, "melanoma_CD8_C1_minus_C4", seed = analysis_seed)

flatten_gsea <- function(x) {
  x <- copy(as.data.table(x))
  x[, leadingEdge_genes := vapply(seq_len(.N), function(i) paste(as.character(leadingEdge[[i]]), collapse = ";"), character(1))]
  x[, leadingEdge := NULL]
  x[]
}

gsea_repro_one <- function(stats, context, time = NA_character_) {
  a <- run_fgsea(stats, all_sets, context, time, seed = analysis_seed)
  b <- run_fgsea(stats, all_sets, context, time, seed = boundary_seed)
  out <- merge(a[, .(feature_id, NES_a = NES, FDR_a = padj, p_a = pval)], b[, .(feature_id, NES_b = NES, FDR_b = padj, p_b = pval)], by = "feature_id", all = TRUE)
  out[, `:=`(context = context, time = time, seed_a = analysis_seed, seed_b = boundary_seed, abs_NES_delta = abs(NES_a - NES_b), abs_FDR_delta = abs(FDR_a - FDR_b), boundary_a = FDR_a >= 0.02 & FDR_a <= 0.10, boundary_b = FDR_b >= 0.02 & FDR_b <= 0.10, boundary_changed = (FDR_a < 0.05) != (FDR_b < 0.05))]
  out[]
}

gsea_repro <- rbindlist(c(
  list(gsea_repro_one(mel_stats, "melanoma_CD8_C1_minus_C4")),
  lapply(c("1h", "6h", "24h", "48h"), function(tm) { d <- bulk_de_stats(tm); stats <- d$rank_score; names(stats) <- d$gene; gsea_repro_one(stats, "bulk_LA_minus_NC", tm) })
), fill = TRUE)

# Program-level ssGSEA bridge with explicit donor summaries and t-based intervals.
bulk_primary <- bulk$score_contrasts[contrast == "LA_minus_NC"]
mel_primary <- mel$score_contrasts
bridge <- merge(
  mel_primary[, .(feature_id, melanoma_logFC = logFC, melanoma_p = P.Value, melanoma_fdr = adj.P.Val, melanoma_CI_low = CI_low, melanoma_CI_high = CI_high)],
  bulk_primary[, .(feature_id, time, bulk_logFC = logFC, bulk_p = P.Value, bulk_fdr = adj.P.Val, bulk_CI_low = CI_low, bulk_CI_high = CI_high)],
  by = "feature_id", allow.cartesian = TRUE
)
donor_support <- merge(bulk$donor_delta[contrast == "LA_minus_NC"], mel_primary[, .(feature_id, melanoma_logFC = logFC)], by = "feature_id", all.x = TRUE)
donor_support[, donor_relation := fifelse(donor_delta == 0 | melanoma_logFC == 0, "zero", fifelse(sign(donor_delta) == sign(melanoma_logFC), "same", "opposite"))]
support_summary <- donor_support[, .(n_donors = uniqueN(donor), n_donors_same_direction = sum(donor_relation == "same"), n_donors_opposite_direction = sum(donor_relation == "opposite"), n_donors_zero = sum(donor_relation == "zero")), by = .(feature_id, time)]
support_summary[, `:=`(donor_same_fraction = n_donors_same_direction / n_donors, donor_relation_consistency = fifelse(n_donors_same_direction == n_donors, "same_3of3", fifelse(n_donors_opposite_direction == n_donors, "opposite_3of3", fifelse(pmax(n_donors_same_direction, n_donors_opposite_direction) >= 2, "2of3_consistent", "mixed"))))]
bridge <- merge(bridge, support_summary, by = c("feature_id", "time"), all.x = TRUE)
bridge[, `:=`(feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]), melanoma_direction = fifelse(melanoma_logFC > 0, "positive", fifelse(melanoma_logFC < 0, "negative", "zero")), bulk_direction = fifelse(bulk_logFC > 0, "positive", fifelse(bulk_logFC < 0, "negative", "zero")))]
bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
classify_bridge <- function(mel_fdr, bulk_fdr, mel_p, bulk_p, relation, n_same, n_opp) {
  if (is.na(relation) || relation == "insufficient_evidence") return("insufficient_evidence")
  if (mel_fdr < 0.05 && bulk_fdr < 0.05 && relation == "same_direction" && n_same == 3) return("supported_concordance")
  if (mel_fdr < 0.05 && bulk_fdr < 0.05 && relation == "opposite_direction" && n_opp == 3) return("supported_divergence")
  if (relation == "same_direction" && (mel_p < 0.05 || bulk_p < 0.05) && n_same >= 2) return("suggestive_concordance")
  if (relation == "opposite_direction" && (mel_p < 0.05 || bulk_p < 0.05) && n_opp >= 2) return("suggestive_divergence")
  if (max(n_same, n_opp) >= 2) return("donor_consistent_direction_only")
  "time_dependent_or_mixed"
}
bridge[, evidence_class := mapply(classify_bridge, melanoma_fdr, bulk_fdr, melanoma_p, bulk_p, direction_relation, n_donors_same_direction, n_donors_opposite_direction)]
bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(bridge, layer_order, feature_id, time)
bridge[, layer_order := NULL]
bridge_summary <- bridge[, .(bridge_class = if (any(evidence_class == "supported_concordance")) "supported_concordance_present" else if (any(evidence_class == "supported_divergence")) "supported_divergence_present" else if (any(evidence_class %in% c("suggestive_concordance", "suggestive_divergence"))) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), mean_abs_bulk_logFC = mean(abs(bulk_logFC), na.rm = TRUE), melanoma_logFC = first(melanoma_logFC), melanoma_fdr = first(melanoma_fdr)), by = .(feature_id, feature_label, feature_layer)]

# GSEA bridge uses formal bulk/melanoma ranks plus donor-specific directional consistency.
gsea_bridge <- merge(
  mel_fgsea_raw[, .(feature_id, melanoma_NES = NES, melanoma_p = pval, melanoma_fdr = padj, melanoma_size = size)],
  bulk_fgsea_raw[, .(feature_id, time, bulk_NES = NES, bulk_p = pval, bulk_fdr = padj, bulk_size = size)],
  by = "feature_id", allow.cartesian = TRUE
)
gsea_donor <- merge(bulk_donor_gsea_raw[, .(feature_id, time, donor, donor_NES = NES)], mel_fgsea_raw[, .(feature_id, melanoma_NES = NES)], by = "feature_id", all.x = TRUE)
gsea_donor[, donor_relation := fifelse(donor_NES == 0 | melanoma_NES == 0, "zero", fifelse(sign(donor_NES) == sign(melanoma_NES), "same", "opposite"))]
gsea_donor_summary <- gsea_donor[, .(n_donors = uniqueN(donor), n_donors_same_direction = sum(donor_relation == "same"), n_donors_opposite_direction = sum(donor_relation == "opposite"), n_donors_zero = sum(donor_relation == "zero")), by = .(feature_id, time)]
gsea_donor_summary[, `:=`(donor_same_fraction = n_donors_same_direction / n_donors, donor_relation_consistency = fifelse(n_donors_same_direction == n_donors, "same_3of3", fifelse(n_donors_opposite_direction == n_donors, "opposite_3of3", fifelse(pmax(n_donors_same_direction, n_donors_opposite_direction) >= 2, "2of3_consistent", "mixed"))))]
gsea_bridge <- merge(gsea_bridge, gsea_donor_summary, by = c("feature_id", "time"), all.x = TRUE)
gsea_bridge[, `:=`(feature_label = unname(set_label[feature_id]), feature_layer = unname(set_layer[feature_id]), melanoma_direction = fifelse(melanoma_NES > 0, "positive", fifelse(melanoma_NES < 0, "negative", "zero")), bulk_direction = fifelse(bulk_NES > 0, "positive", fifelse(bulk_NES < 0, "negative", "zero")))]
gsea_bridge[, direction_relation := fifelse(melanoma_direction == "zero" | bulk_direction == "zero", "insufficient_evidence", fifelse(melanoma_direction == bulk_direction, "same_direction", "opposite_direction"))]
gsea_bridge[, evidence_class := mapply(classify_bridge, melanoma_fdr, bulk_fdr, melanoma_p, bulk_p, direction_relation, n_donors_same_direction, n_donors_opposite_direction)]
gsea_bridge[, layer_order := match(feature_layer, c("LM12_FIG2_FROZEN", "MSIGDB_FULL_SET"))]
setorder(gsea_bridge, layer_order, feature_id, time)
gsea_bridge[, layer_order := NULL]
gsea_bridge_summary <- gsea_bridge[, .(bridge_class = if (any(evidence_class == "supported_concordance")) "supported_concordance_present" else if (any(evidence_class == "supported_divergence")) "supported_divergence_present" else if (any(evidence_class %in% c("suggestive_concordance", "suggestive_divergence"))) "suggestive_or_mixed" else "insufficient_or_mixed", n_timepoints = .N, n_same = sum(direction_relation == "same_direction"), n_opposite = sum(direction_relation == "opposite_direction"), melanoma_NES = first(melanoma_NES), melanoma_fdr = first(melanoma_fdr), mean_bulk_abs_NES = mean(abs(bulk_NES), na.rm = TRUE)), by = .(feature_id, feature_label, feature_layer)]

# 方法方向审计：同一 feature 的 ssGSEA 与 GSEA 不强行视为同一结果。
method_direction_audit <- merge(
  mel_primary[, .(feature_id, context = "melanoma_CD8", time = NA_character_, ssgsea_effect = logFC, ssgsea_direction = sign(logFC))],
  mel_fgsea_raw[, .(feature_id, gsea_effect = NES, gsea_direction = sign(NES))], by = "feature_id", all = TRUE
)
method_direction_audit[, `:=`(method_direction_discordant = !is.na(ssgsea_direction) & !is.na(gsea_direction) & ssgsea_direction != 0 & gsea_direction != 0 & ssgsea_direction != gsea_direction, interpretation = "ssGSEA is sample-level score; GSEA is rank enrichment; discordance requires review")]
bulk_method <- merge(bulk_primary[, .(feature_id, time, ssgsea_effect = logFC, ssgsea_direction = sign(logFC))], bulk_fgsea_raw[, .(feature_id, time, gsea_effect = NES, gsea_direction = sign(NES))], by = c("feature_id", "time"), all = TRUE)
bulk_method[, `:=`(context = "bulk_LA_minus_NC", method_direction_discordant = !is.na(ssgsea_direction) & !is.na(gsea_direction) & ssgsea_direction != 0 & gsea_direction != 0 & ssgsea_direction != gsea_direction, interpretation = "ssGSEA is sample-level score; GSEA is rank enrichment; discordance requires review")]
method_direction_audit <- rbind(method_direction_audit, bulk_method, fill = TRUE)
memory_direction_audit <- method_direction_audit[feature_id %in% c("MSIGDB::GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN", "MSIGDB::GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_UP")]
memory_direction_audit[, semantic_note := fifelse(grepl("_DN$", feature_id), "DN set is lower in effector than memory; interpret as the memory-enriched side of the paired signature", "UP set is higher in effector than memory; interpret as the effector-enriched side of the paired signature")]

# Leading-edge 拆解使用本地完整 MSigDB 诊断集合，不另造 CD8 marker 集合。
diag_specs <- data.table(
  category = c("cell_cycle_E2F", "cell_cycle_G2M", "IFN_gamma", "T_cell_cytotoxicity", "TCR_signaling", "hypoxia", "UPR", "glycolysis"),
  source_file_key = c("hallmark", "hallmark", "hallmark", "go_bp", "reactome", "hallmark", "hallmark", "hallmark"),
  feature_id = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "GOBP_T_CELL_MEDIATED_CYTOTOXICITY", "REACTOME_TCR_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_GLYCOLYSIS")
)
diag_specs[, exists := mapply(function(key, id) id %in% names(gmt_cache[[key]]$sets), source_file_key, feature_id)]
diag_specs[, n_genes_source := mapply(function(key, id) if (id %in% names(gmt_cache[[key]]$sets)) length(gmt_cache[[key]]$sets[[id]]) else NA_integer_, source_file_key, feature_id)]
diag_specs[, source_path := unname(msigdb_files[source_file_key])]
diag_existing <- diag_specs[exists == TRUE]
diag_sets <- setNames(lapply(seq_len(nrow(diag_existing)), function(i) gmt_cache[[diag_existing$source_file_key[i]]]$sets[[diag_existing$feature_id[i]]]), diag_existing$category)
gene_to_cat <- list()
for (cat in names(diag_sets)) for (g in diag_sets[[cat]]) gene_to_cat[[g]] <- unique(c(gene_to_cat[[g]], cat))

gsea_primary_raw <- rbind(mel_fgsea_raw, bulk_fgsea_raw, fill = TRUE)
edge_rows <- list(); k <- 0L
for (i in seq_len(nrow(gsea_primary_raw))) {
  genes <- as.character(gsea_primary_raw$leadingEdge[[i]])
  if (!length(genes)) next
  k <- k + 1L
  edge_rows[[k]] <- data.table(context = gsea_primary_raw$context[i], time = gsea_primary_raw$time[i], feature_id = gsea_primary_raw$feature_id[i], gene = genes, rank_in_leading_edge = seq_along(genes))
}
edge_long <- rbindlist(edge_rows, fill = TRUE)
if (nrow(edge_long)) {
  edge_cat <- rbindlist(lapply(seq_len(nrow(edge_long)), function(i) {
    cats <- gene_to_cat[[edge_long$gene[i]]]; if (is.null(cats) || !length(cats)) cats <- "unclassified"
    data.table(context = edge_long$context[i], time = edge_long$time[i], feature_id = edge_long$feature_id[i], gene = edge_long$gene[i], rank_in_leading_edge = edge_long$rank_in_leading_edge[i], category = cats)
  }), fill = TRUE)
  edge_summary <- edge_cat[, .(n_genes = uniqueN(gene), mean_rank = mean(rank_in_leading_edge), genes = paste(head(unique(gene), 25L), collapse = ";")), by = .(context, time, feature_id, category)]
} else {
  edge_cat <- data.table(); edge_summary <- data.table()
}

# 基因级桥接只用于候选优先级；不把名义方向直接称为验证。
mel_gene <- mel$deg_fixed[order(fdr, p_value)][!duplicated(gene)]
gene_bridge <- merge(mel_gene[, .(gene, melanoma_logFC = log2fc_c1_minus_ref, melanoma_p = p_value, melanoma_fdr = fdr)], rbindlist(lapply(c("1h", "6h", "24h", "48h"), function(tm) { d <- bulk_de_stats(tm); d[, time := tm]; d[gene %in% unique(unlist(all_sets))] })), by = "gene", allow.cartesian = TRUE)
gene_bridge[, direction_relation := fifelse(sign(melanoma_logFC) == sign(logFC), "same_direction", "opposite_direction")]
gene_bridge[, set_membership := vapply(gene, function(g) paste(names(all_sets)[vapply(all_sets, function(s) g %in% s, logical(1))], collapse = ";"), character(1))]
gene_bridge[, candidate_gate := fifelse(melanoma_fdr < 0.05 & FDR < 0.05 & direction_relation == "same_direction", "strict_two_sided_gate", "not_passed")]

membership <- merge(feature_meta, source_size, by = "feature_id", all = TRUE, sort = FALSE)
membership[, source_path := unname(source_path[feature_id])]
membership <- merge(membership, bulk$coverage[, .(feature_id, n_genes_bulk = n_genes_present, bulk_coverage = coverage)], by = "feature_id", all = TRUE)
membership <- merge(membership, mel$coverage[, .(feature_id, n_genes_melanoma = n_genes_present, melanoma_coverage = coverage)], by = "feature_id", all = TRUE)
membership[, `:=`(bulk_score_status = ifelse(n_genes_bulk >= 3L, "PASS", "HOLD"), melanoma_score_status = ifelse(n_genes_melanoma >= 3L, "PASS", "HOLD"))]

write_tsv(membership, file.path(out_dir, "fixed_gene_set_coverage_manifest.tsv"))
write_tsv(rbind(bulk$mapping, mel$mapping, fill = TRUE), file.path(out_dir, "score_id_mapping_qa.tsv"))
write_tsv(bulk$score_long, file.path(out_dir, "bulk_fixed_readout_ssgsea_long.tsv"))
write_tsv(bulk$score_contrasts, file.path(out_dir, "bulk_fixed_readout_ssgsea_contrasts.tsv"))
write_tsv(bulk$donor_delta, file.path(out_dir, "bulk_fixed_readout_donor_deltas.tsv"))
write_tsv(bulk$donor_summary, file.path(out_dir, "bulk_fixed_readout_donor_summary.tsv"))
write_tsv(mel$score_long, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_long.tsv"))
write_tsv(mel$score_contrasts, file.path(out_dir, "melanoma_cd8_fixed_readout_ssgsea_c1_vs_c4.tsv"))
write_tsv(mel$group_summary, file.path(out_dir, "melanoma_cd8_fixed_readout_group_summary.tsv"))
write_tsv(mel$selection_manifest, file.path(out_dir, "melanoma_cd8_patient_selection_manifest.tsv"))
write_tsv(mel$deg_fixed, file.path(out_dir, "melanoma_cd8_fixed_readout_gene_deg.tsv"))
write_tsv(bridge, file.path(out_dir, "melanoma_cd8_bulk_bridge_by_time.tsv"))
write_tsv(bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_bridge_summary.tsv"))
write_tsv(gene_bridge, file.path(out_dir, "melanoma_cd8_bulk_gene_bridge.tsv"))
write_tsv(flatten_gsea(bulk_fgsea_raw), file.path(out_dir, "bulk_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(mel_fgsea_raw), file.path(out_dir, "melanoma_cd8_all_fixed_preranked_gsea.tsv"))
write_tsv(flatten_gsea(bulk_donor_gsea_raw), file.path(out_dir, "bulk_donor_directional_preranked_gsea.tsv"))
write_tsv(gsea_bridge, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_by_time.tsv"))
write_tsv(gsea_bridge_summary, file.path(out_dir, "melanoma_cd8_bulk_gsea_bridge_summary.tsv"))
write_tsv(gsea_repro, file.path(out_dir, "gsea_reproducibility_boundary.tsv"))
write_tsv(method_direction_audit, file.path(out_dir, "ssgsea_gsea_direction_audit.tsv"))
write_tsv(memory_direction_audit, file.path(out_dir, "memory_direction_audit.tsv"))
write_tsv(diag_specs, file.path(out_dir, "diagnostic_gene_set_manifest.tsv"))
write_tsv(edge_long, file.path(out_dir, "gsea_leading_edge_long.tsv"))
write_tsv(edge_cat, file.path(out_dir, "gsea_leading_edge_category_long.tsv"))
write_tsv(edge_summary, file.path(out_dir, "gsea_leading_edge_category_summary.tsv"))
write_tsv(gsea_repro[boundary_a | boundary_b | boundary_changed], file.path(out_dir, "gsea_boundary_focus.tsv"))
write_tsv(data.table(
  object = c("bulk_count", "melanoma_cd8_pseudobulk"),
  analysis_unit = c("donor x condition x time sample", "sample x native CD8 T pseudobulk; patient biological unit"),
  primary_comparison = c("time-matched LA_minus_NC", "C1_minus_C4"),
  n_independent_units = c(length(unique(bulk$sample_sheet$Donor)), mel$n_patients),
  interpretation = c("CD8 bulk aggregate response", "melanoma CD8 state association")
), file.path(out_dir, "analysis_unit_manifest.tsv"))

mapping_ok <- all(bulk$mapping$mapping_ok) && all(mel$mapping$mapping_ok)
gsea_rows_ok <- nrow(bulk_fgsea_raw) == length(all_sets) * 4L && nrow(mel_fgsea_raw) == length(all_sets)
donor_rows_ok <- nrow(bulk$donor_delta[contrast == "LA_minus_NC"]) == length(all_sets) * 4L * 3L
repro_ok <- nrow(gsea_repro) > 0L && all(is.finite(gsea_repro$abs_NES_delta)) && all(is.finite(gsea_repro$abs_FDR_delta))
qa <- data.table(
  check_id = c("bulk_factorial_samples", "bulk_design_rank", "melanoma_c1c4_samples", "melanoma_c1c4_patients", "fixed_readout_count", "score_id_mapping_exact", "bulk_coverage", "melanoma_coverage", "all_ssgsea_finite", "bulk_primary_donor_rows", "gsea_rows_all_fixed", "donor_gsea_rows", "gsea_seed_reproducible_run", "leading_edge_table"),
  status = c(
    ifelse(nrow(bulk$sample_sheet) == 36L, "PASS", "HOLD"),
    ifelse(bulk$design_rank == ncol(bulk$design), "PASS", "HOLD"),
    ifelse(mel$n_samples >= 6L, "PASS", "HOLD"),
    ifelse(mel$n_patients >= 6L, "PASS", "HOLD"),
    ifelse(length(all_sets) == 21L, "PASS", "HOLD"),
    ifelse(mapping_ok, "PASS", "HOLD"),
    ifelse(all(membership$bulk_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(membership$melanoma_score_status == "PASS"), "PASS", "HOLD"),
    ifelse(all(is.finite(bulk$score_long$ssgsea_score)) && all(is.finite(mel$score_long$ssgsea_score)), "PASS", "HOLD"),
    ifelse(donor_rows_ok, "PASS", "HOLD"),
    ifelse(gsea_rows_ok, "PASS", "HOLD"),
    ifelse(nrow(bulk_donor_gsea_raw) == length(all_sets) * 4L * 3L, "PASS", "HOLD"),
    ifelse(repro_ok, "PASS", "HOLD"),
    ifelse(nrow(edge_long) > 0L, "PASS", "HOLD")
  ),
  detail = c(nrow(bulk$sample_sheet), paste(bulk$design_rank, "of", ncol(bulk$design)), mel$n_samples, mel$n_patients, length(all_sets), paste(all(bulk$mapping$mapping_ok), all(mel$mapping$mapping_ok), sep = ";"), sprintf("min %.3f", min(membership$bulk_coverage)), sprintf("min %.3f", min(membership$melanoma_coverage)), paste(nrow(bulk$score_long), nrow(mel$score_long), sep = " + "), nrow(bulk$donor_delta[contrast == "LA_minus_NC"]), paste(nrow(bulk_fgsea_raw), nrow(mel_fgsea_raw), sep = " + "), nrow(bulk_donor_gsea_raw), sprintf("max |delta NES| %.4g; max |delta FDR| %.4g", max(gsea_repro$abs_NES_delta), max(gsea_repro$abs_FDR_delta)), nrow(edge_long))
)
write_tsv(qa, file.path(out_dir, "qa_summary.tsv"))
overall_status <- if (all(qa$status == "PASS")) "PASS_CANDIDATE" else "HOLD_CANDIDATE"
writeLines(c(paste0("status=", overall_status), paste0("output_dir=", normalizePath(out_dir, mustWork = FALSE)), paste0("analysis_seed=", analysis_seed), paste0("boundary_seed=", boundary_seed), paste0("fixed_readouts=", length(all_sets)), paste0("bulk_samples=", nrow(bulk$sample_sheet), ";bulk_donors=", length(unique(bulk$sample_sheet$Donor))), paste0("melanoma_cd8_samples=", mel$n_samples, ";melanoma_patients=", mel$n_patients, ";selection=", mel$selection_mode), paste0("gsea_rows=", nrow(bulk_fgsea_raw), "+", nrow(mel_fgsea_raw)), paste0("leading_edge_rows=", nrow(edge_long))), file.path(out_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("STATUS\t", overall_status, "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
