#!/usr/bin/env Rscript
# Purpose: Fit donor-adjusted CD8 treatment-time responses and fixed programme enrichment.
# Task: cd8_response; step: s02
# Inputs (DCF fields): count_rdata
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task cd8_response --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 02_differential_response.R --config PATH [--check-config]\n       Rscript 02_differential_response.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("count_rdata","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths["count_rdata"]
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

# ---- Retained complete source body; CLI/path substitutions only ----
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

input_file <- unname(paths["count_rdata"])
out_dir <- unname(paths["output_dir"])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path, row.names = FALSE) {
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE, row.names = row.names,
    na = "NA", col.names = TRUE
  )
}
stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

stop_if(!requireNamespace("edgeR", quietly = TRUE), "R 包 edgeR 未安装。")

e <- new.env(parent = emptyenv())
stop_if(!file.exists(input_file), paste("输入文件不存在:", input_file))
load(input_file, envir = e)
stop_if(!all(c("count", "sample_anno") %in% ls(e, all.names = TRUE)),
        "RData 必须同时包含 count 和 sample_anno。")
count <- get("count", envir = e)
sample_anno <- get("sample_anno", envir = e)
stop_if(!identical(colnames(count), rownames(sample_anno)),
        "count 列名与 sample_anno 行名不一致。")

anno <- sample_anno[colnames(count), c("Donor", "Treatment", "Time"), drop = FALSE]
anno$SampleID <- rownames(anno)
factorial <- as.character(anno$Treatment) %in% c("HCl", "LA", "NC")
fa <- anno[factorial, , drop = FALSE]
fa$Donor <- factor(as.character(fa$Donor))
fa$Treatment <- factor(as.character(fa$Treatment), levels = c("NC", "HCl", "LA"))
time_levels <- c("1h", "6h", "24h", "48h")
stop_if(!all(as.character(fa$Time) %in% time_levels),
        "发现未登记的 factorial Time 水平。")
fa$Time <- factor(as.character(fa$Time), levels = time_levels)
stop_if(anyNA(fa$Donor) || anyNA(fa$Treatment) || anyNA(fa$Time),
        "factorial 设计因子包含 NA。")

design <- stats::model.matrix(~ Donor + Time + Treatment + Time:Treatment, data = fa)
colnames(design) <- make.unique(colnames(design))
design_rank <- qr(design)$rank
stop_if(design_rank != ncol(design), "差异模型设计矩阵不满秩。")

dge <- edgeR::DGEList(counts = count[, fa$SampleID, drop = FALSE], samples = fa)
keep <- edgeR::filterByExpr(dge, design = design)
stop_if(sum(keep) < 1000L, "过滤后基因数过少。")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")
dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
fit <- edgeR::glmQLFit(dge, design, robust = TRUE)

coef_names <- colnames(design)
find_coef <- function(term) {
  hit <- which(coef_names == term)
  if (length(hit) != 1L) {
    hit <- which(coef_names %in% c(term, sub(":", ":", term, fixed = TRUE)))
  }
  stop_if(length(hit) != 1L, paste("找不到唯一模型系数:", term))
  hit
}

contrast_vector <- function(left, right, time) {
  v <- numeric(ncol(design))
  add_mean <- function(treatment, sign) {
    base_term <- paste0("Treatment", treatment)
    if (treatment != "NC") v[find_coef(base_term)] <<- v[find_coef(base_term)] + sign
    if (time != "1h" && treatment != "NC") {
      interaction_terms <- c(
        paste0("Time", time, ":Treatment", treatment),
        paste0("Treatment", treatment, ":Time", time)
      )
      hit <- which(coef_names %in% interaction_terms)
      stop_if(length(hit) != 1L, paste("找不到时间交互项:", paste(interaction_terms, collapse = "/")))
      v[hit] <<- v[hit] + sign
    }
  }
  add_mean(left, 1)
  add_mean(right, -1)
  v
}

contrast_specs <- expand.grid(
  contrast = c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"),
  time = time_levels,
  stringsAsFactors = FALSE
)
contrast_specs$left <- c("LA", "LA", "HCl")[match(contrast_specs$contrast,
                                                    c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))]
contrast_specs$right <- c("NC", "HCl", "NC")[match(contrast_specs$contrast,
                                                       c("LA_minus_NC", "LA_minus_HCl", "HCl_minus_NC"))]
contrast_specs$contrast_vector <- lapply(seq_len(nrow(contrast_specs)), function(i) {
  contrast_vector(contrast_specs$left[[i]], contrast_specs$right[[i]], contrast_specs$time[[i]])
})

format_de <- function(test, spec) {
  tab <- as.data.frame(edgeR::topTags(test, n = Inf, sort.by = "none")$table,
                       stringsAsFactors = FALSE)
  tab$gene_id <- rownames(tab)
  tab <- tab[, c("gene_id", setdiff(colnames(tab), "gene_id")), drop = FALSE]
  tab$contrast <- spec$contrast
  tab$time <- spec$time
  tab$left <- spec$left
  tab$right <- spec$right
  tab
}

contrast_manifest <- do.call(rbind, lapply(seq_len(nrow(contrast_specs)), function(i) {
  v <- contrast_specs$contrast_vector[[i]]
  data.frame(
    contrast_id = paste(contrast_specs$contrast[[i]], contrast_specs$time[[i]], sep = "__"),
    contrast = contrast_specs$contrast[[i]], time = contrast_specs$time[[i]],
    left = contrast_specs$left[[i]], right = contrast_specs$right[[i]],
    n_left = sum(fa$Treatment == contrast_specs$left[[i]] & fa$Time == contrast_specs$time[[i]]),
    n_right = sum(fa$Treatment == contrast_specs$right[[i]] & fa$Time == contrast_specs$time[[i]]),
    contrast_nonzero = paste(paste(coef_names[which(v != 0)], v[v != 0], sep = "="), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

for (i in seq_len(nrow(contrast_specs))) {
  spec <- contrast_specs[i, , drop = FALSE]
  test <- edgeR::glmQLFTest(fit, contrast = contrast_specs$contrast_vector[[i]])
  result <- format_de(test, spec)
  out_name <- paste0("de_", spec$contrast, "_", spec$time, ".tsv")
  write_tsv(result, file.path(out_dir, out_name))
}

interaction_terms <- grep("^Time.*:Treatment|^Treatment.*:Time", coef_names)
stop_if(length(interaction_terms) != 6L, "时间×处理交互项数量不是预期的 6。")
interaction_test <- edgeR::glmQLFTest(fit, coef = interaction_terms)
interaction_result <- as.data.frame(edgeR::topTags(interaction_test, n = Inf, sort.by = "none")$table,
                                    stringsAsFactors = FALSE)
interaction_result$gene_id <- rownames(interaction_result)
interaction_result <- interaction_result[, c("gene_id", setdiff(colnames(interaction_result), "gene_id")), drop = FALSE]
write_tsv(interaction_result, file.path(out_dir, "time_treatment_interaction.tsv"))

sample_design <- data.frame(fa, stringsAsFactors = FALSE)
write_tsv(sample_design, file.path(out_dir, "factorial_sample_design.tsv"))
write_tsv(data.frame(coefficient = coef_names, index = seq_along(coef_names), stringsAsFactors = FALSE),
          file.path(out_dir, "model_coefficients.tsv"))
write_tsv(contrast_manifest, file.path(out_dir, "contrast_manifest.tsv"))
write_tsv(data.frame(
  input_file = normalizePath(input_file, mustWork = TRUE),
  factorial_samples = nrow(fa), genes_before_filter = nrow(count),
  genes_after_filter = sum(keep), design_rank = design_rank,
  design_columns = ncol(design), n_interaction_terms = length(interaction_terms),
  stringsAsFactors = FALSE
), file.path(out_dir, "run_summary.tsv"))
write_tsv(data.frame(
  check_id = c("factorial_samples", "design_full_rank", "genes_after_filter",
               "all_time_contrasts", "interaction_test"),
  status = c(
    ifelse(nrow(fa) == 36L, "PASS", "HOLD"),
    ifelse(design_rank == ncol(design), "PASS", "HOLD"),
    ifelse(sum(keep) >= 1000L, "PASS", "HOLD"),
    ifelse(nrow(contrast_specs) == 12L, "PASS", "HOLD"),
    ifelse(length(interaction_terms) == 6L, "PASS", "HOLD")
  ),
  detail = c(nrow(fa), paste(design_rank, "of", ncol(design)), sum(keep),
             nrow(contrast_specs), length(interaction_terms)),
  stringsAsFactors = FALSE
), file.path(out_dir, "qa_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

cat("DE_STATUS\tPASS\n")
cat("GENES_AFTER_FILTER\t", sum(keep), "\n", sep = "")
cat("CONTRAST_FILES\t", nrow(contrast_specs), "\n", sep = "")
cat("OUTPUT_DIR\t", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")

