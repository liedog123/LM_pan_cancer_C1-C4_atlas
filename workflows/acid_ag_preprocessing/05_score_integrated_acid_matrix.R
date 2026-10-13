#!/usr/bin/env Rscript
# Fig. 1C核心计算；路径检查不加载研究输入。

options(stringsAsFactors = FALSE)
usage <- function() {
  cat("Usage: Rscript 05_score_integrated_acid_matrix.R --config PATH [--check-config]\n",
      "       Rscript 05_score_integrated_acid_matrix.R --help\n\n",
      "DCF fields: clinical_groups, final_genes, integrated_expression, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "Scores the full integrated expression matrix using existing acid groups; no clustering or plotting.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (!(length(args) %in% c(2L, 3L)) || args[1L] != "--config" || !nzchar(args[2L]) ||
    (length(args) == 3L && args[3L] != "--check-config")) {
  usage()
  stop("Expected --config PATH, optionally followed by --check-config.", call. = FALSE)
}
flags <- if (length(args) == 3L) args[3L] else character()

# Resolve paths relative to the DCF directory, including existing ancestors and symlinks.
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
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
config_path <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path, mode = 4L) != 0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size <= 0) {
  stop("Configuration must be a readable nonempty regular file.", call. = FALSE)
}
cfg <- read.dcf(config_path)
required_fields <- c("clinical_groups", "final_genes", "integrated_expression", "output_dir")
if (nrow(cfg) != 1L || !setequal(colnames(cfg), required_fields) ||
    ncol(cfg) != length(required_fields)) {
  stop("Configuration must contain one DCF record with exactly: ",
       paste(required_fields, collapse = ", "), call. = FALSE)
}
cfg <- setNames(as.character(cfg[1L, required_fields]), required_fields)
if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME", cfg, fixed = TRUE))) {
  stop("Replace all example placeholders with explicit paths.", call. = FALSE)
}
paths <- vapply(cfg, resolve_path, character(1L), base_dir = dirname(config_path))
required_files <- unname(paths[setdiff(names(paths), "output_dir")])
out_dir <- unname(paths["output_dir"])


# Protect source inputs and existing results; check mode creates no directories.
protected_paths <- c(dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
if ("--check-config" %in% flags) {
  cat("Configuration paths passed metadata checks. No inputs loaded; no analysis or outputs created.\n")
  quit(status = 0L)
}

# Isolate each serialized input and retrieve only the source's named object.
load_named_object <- function(path, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(path, envir = input_env)
  if (!exists(object_name, envir = input_env, inherits = FALSE)) {
    stop("Required object not found: ", object_name, " in ", path, call. = FALSE)
  }
  get(object_name, envir = input_env, inherits = FALSE)
}
if (!dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) && !dir.exists(out_dir)) {
  stop("Cannot create output_dir.", call. = FALSE)
}
suppressPackageStartupMessages({library(data.table); library(dplyr); library(GSVA)})
clinical_data_cluster_list <- load_named_object(paths["clinical_groups"], "clinical_data_cluster_list")
final_selected_genes <- load_named_object(paths["final_genes"], "final_selected_genes")
expr_file <- unname(paths["integrated_expression"])
output_rdata <- file.path(out_dir, "pancancer_integrated_ssgsea_acid_scores.rdata")
output_csv <- file.path(out_dir, "pancancer_integrated_ssgsea_acid_scores_by_sample.csv")

all_clinical_data <- bind_rows(clinical_data_cluster_list) %>%
  select(sample, Group, `cancer type abbreviation`) %>%
  rename(Cancer = `cancer type abbreviation`)

message("Loading integrated pan-cancer matrix...")
expr <- fread(expr_file, sep = "\t", header = TRUE, data.table = FALSE)
if (sum(duplicated(expr$sample)) > 0) {
  expr <- expr[!duplicated(expr$sample), ]
}
rownames(expr) <- expr$sample
expr <- expr[, -1, drop = FALSE]

common_samples <- intersect(all_clinical_data$sample, colnames(expr))
all_clinical_data <- all_clinical_data %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, .keep_all = TRUE)

expr_matrix <- as.matrix(expr[, all_clinical_data$sample, drop = FALSE])

gene_list <- as.character(unique(final_selected_genes$Gene))
alias_map <- c(PKM = "PKM2")
mapped_gene_list <- vapply(
  gene_list,
  function(gene) {
    if (gene %in% rownames(expr_matrix)) {
      gene
    } else if (!is.na(alias_map[gene]) && alias_map[gene] %in% rownames(expr_matrix)) {
      unname(alias_map[gene])
    } else {
      NA_character_
    }
  },
  character(1)
)

missing_genes <- gene_list[is.na(mapped_gene_list)]
if (length(missing_genes) > 0) {
  stop("Integrated matrix is missing genes: ", paste(missing_genes, collapse = ", "))
}

if (any(mapped_gene_list != gene_list)) {
  message(
    "Gene alias mapping applied: ",
    paste(gene_list[mapped_gene_list != gene_list], mapped_gene_list[mapped_gene_list != gene_list], sep = " -> ", collapse = "; ")
  )
}

message("Running ssGSEA on integrated matrix...")
ssgsea_param <- ssgseaParam(
  exprData = expr_matrix,
  geneSets = list(Acid_Score = unname(mapped_gene_list)),
  normalize = TRUE
)
ssgsea_scores <- gsva(ssgsea_param, verbose = FALSE)

sample_scores_df <- data.frame(
  sample = colnames(ssgsea_scores),
  Acid_Score = as.numeric(ssgsea_scores[1, ]),
  stringsAsFactors = FALSE
) %>%
  left_join(all_clinical_data, by = "sample")

stopifnot(!any(is.na(sample_scores_df$Acid_Score)))
stopifnot(!any(is.na(sample_scores_df$Cancer)))
stopifnot(!any(is.na(sample_scores_df$Group)))

pancancer_integrated_ssgsea_acid_scores <- list(
  score_matrix = ssgsea_scores,
  sample_scores = sample_scores_df,
  metadata = list(
    method = "ssGSEA",
    normalize = TRUE,
    input_matrix = expr_file,
    input_scale = "EB++AdjustPANCAN integrated matrix used directly",
    gene_count = length(gene_list),
    gene_mapping = stats::setNames(as.list(mapped_gene_list), gene_list),
    sample_count = nrow(sample_scores_df),
    cancer_count = dplyr::n_distinct(sample_scores_df$Cancer),
    mode = "full"
  )
)

save(pancancer_integrated_ssgsea_acid_scores, file = output_rdata)
write.csv(sample_scores_df, output_csv, row.names = FALSE)

summary_df <- sample_scores_df %>%
  group_by(Cancer) %>%
  summarise(
    N = n(),
    Median_Acid = median(Acid_Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Median_Acid))

message("Done.")
message(
  "Range: ",
  sprintf("%.4f", min(sample_scores_df$Acid_Score, na.rm = TRUE)),
  " to ",
  sprintf("%.4f", max(sample_scores_df$Acid_Score, na.rm = TRUE))
)
print(summary_df)

# 沿用panel-source的癌种及High/Low中位摘要；不添加纳入门槛。
top_order_df <- sample_scores_df %>%
  group_by(Cancer) %>%
  summarise(
    Tumor_Median = median(Acid_Score, na.rm = TRUE),
    High_Median = median(Acid_Score[Group == "High"], na.rm = TRUE),
    Low_Median = median(Acid_Score[Group == "Low"], na.rm = TRUE),
    Tumor_N = n(),
    High_N = sum(Group == "High"),
    Low_N = sum(Group == "Low"),
    .groups = "drop"
  ) %>%
  arrange(desc(Tumor_Median)) %>%
  mutate(Order = row_number())

write.csv(summary_df, file.path(out_dir, "pancancer_integrated_ssgsea_acid_summary.csv"), row.names = FALSE)
write.csv(top_order_df, file.path(out_dir, "Fig1C_top_cancer_group_medians.csv"), row.names = FALSE)
