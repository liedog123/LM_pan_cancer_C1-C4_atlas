#!/usr/bin/env Rscript
# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Task: ag_survival; step: s03
# Inputs (DCF fields): survival_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task ag_survival --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 03_export_formal_ag_km.R --config PATH [--check-config]\n",
      "       Rscript 03_export_formal_ag_km.R --help\n\n",
      "DCF fields: survival_data, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation executes the documented AG continuation.\n",
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

# Resolve DCF-relative paths through existing ancestors and symlinks.
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
required_fields <- c("survival_data", "output_dir")
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
input_rdata <- unname(paths["survival_data"])
out_dir <- unname(paths["output_dir"])
required_files <- unname(paths[setdiff(required_fields, "output_dir")])
# Resolve each input's final target before checking output independence.
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)

# Protect inputs and existing results; checks create no directories.
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

suppressPackageStartupMessages({ library(survival) })
load_named <- function(path, name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) stop("Missing named object: ", name)
  get(name, envir = env, inherits = FALSE)
}
survival_data_list <- load_named(input_rdata, "survival_data_list")
if (!is.list(survival_data_list)) stop("survival_data_list is not a list")
qa_dir <- file.path(out_dir, "qa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

expected_cancers <- c("BLCA", "BRCA", "CESC", "HNSC", "KIRC", "KIRP",
                      "LGG", "LIHC", "LUAD", "SARC", "SKCM")
expected_n <- c(BLCA = 248L, BRCA = 418L, CESC = 173L, HNSC = 353L,
                KIRC = 286L, KIRP = 118L, LGG = 336L, LIHC = 263L,
                LUAD = 287L, SARC = 178L, SKCM = 344L)
expected_low <- c(BLCA = 172L, BRCA = 213L, CESC = 95L, HNSC = 91L,
                  KIRC = 129L, KIRP = 67L, LGG = 181L, LIHC = 214L,
                  LUAD = 246L, SARC = 117L, SKCM = 220L)
expected_high <- c(BLCA = 76L, BRCA = 205L, CESC = 78L, HNSC = 262L,
                   KIRC = 157L, KIRP = 51L, LGG = 155L, LIHC = 49L,
                   LUAD = 41L, SARC = 61L, SKCM = 124L)
expected_events <- c(BLCA = 103L, BRCA = 65L, CESC = 39L, HNSC = 162L,
                     KIRC = 107L, KIRP = 22L, LGG = 100L, LIHC = 85L,
                     LUAD = 92L, SARC = 71L, SKCM = 159L)
expected_p <- c(BLCA = 0.0231412, BRCA = 0.000311019, CESC = 0.000247714,
                HNSC = 0.000851538, KIRC = 0.0196141, KIRP = 0.0175173,
                LGG = 1.28265e-06, LIHC = 8.43172e-05, LUAD = 4.66134e-05,
                SARC = 0.00816214, SKCM = 0.00604262)

resolve_name <- function(cancer) {
  candidates <- c(cancer, paste0("TCGA-", cancer))
  hit <- candidates[candidates %in% names(survival_data_list)]
  if (length(hit) != 1L) {
    stop("Cannot uniquely resolve cancer name: ", cancer,
         "; available matches: ", paste(hit, collapse = ","))
  }
  hit[[1L]]
}

as_numeric_safe <- function(x, field, label) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(out)) && any(!is.na(x))) {
    stop("Cannot coerce ", field, " to numeric for ", label)
  }
  out
}

prepare_cancer <- function(cancer) {
  source_name <- resolve_name(cancer)
  d <- survival_data_list[[source_name]]
  required <- c("Group", "OS", "OS.time")
  missing <- setdiff(required, colnames(d))
  if (length(missing) > 0L) stop("Missing fields for ", source_name, ": ", paste(missing, collapse = ","))
  d <- d[d$Group %in% c("High", "Low"), , drop = FALSE]
  d$Group <- factor(as.character(d$Group), levels = c("Low", "High"))
  d$OS <- as_numeric_safe(d$OS, "OS", source_name)
  d$OS.time <- as_numeric_safe(d$OS.time, "OS.time", source_name)
  d <- d[complete.cases(d[, c("Group", "OS", "OS.time")]), , drop = FALSE]
  d$Cancer <- cancer
  d$panel_id <- sprintf("KM%02d_%s", match(cancer, expected_cancers), cancer)
  d$panel_order <- match(cancer, expected_cancers)
  d$display_label <- cancer
  d$source_scope <- "within_cancer"
  d$sample_id <- rownames(d)
  if (is.null(d$sample_id) || anyNA(d$sample_id) || any(d$sample_id == "")) {
    d$sample_id <- paste0(d$Cancer, "_row_", seq_len(nrow(d)))
  }
  d[, c("panel_order", "panel_id", "display_label", "Cancer", "source_scope",
        "sample_id", "Group", "OS", "OS.time")]
}

get_logrank <- function(d) {
  fit <- survdiff(Surv(OS.time, OS) ~ Group, data = d)
  1 - pchisq(fit$chisq, df = length(fit$n) - 1L)
}

within_data <- do.call(rbind, lapply(expected_cancers, prepare_cancer))
rownames(within_data) <- NULL

# 严格复用 KM_total_samples.R：七个生存字段整体 complete-case 后再拟合 OS。
all_parts <- lapply(names(survival_data_list), function(source_name) {
  d <- survival_data_list[[source_name]]
  required <- c("Group", "OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time")
  missing <- setdiff(required, colnames(d))
  if (length(missing) > 0L) stop("Missing pan-cancer fields for ", source_name, ": ", paste(missing, collapse = ","))
  keep <- d$Group %in% c("High", "Low")
  d <- d[keep, required, drop = FALSE]
  d$Group <- factor(as.character(d$Group), levels = c("Low", "High"))
  for (field in setdiff(required, "Group")) d[[field]] <- as_numeric_safe(d[[field]], field, source_name)
  d$Cancer <- source_name
  d$sample_id <- rownames(d)
  if (is.null(d$sample_id) || anyNA(d$sample_id) || any(d$sample_id == "")) {
    d$sample_id <- paste0(source_name, "_row_", seq_len(nrow(d)))
  }
  d
})
all_survival <- do.call(rbind, all_parts)
all_survival <- all_survival[complete.cases(all_survival[, c("Group", "OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time")]), , drop = FALSE]
pan_data <- all_survival[, c("Group", "OS", "OS.time", "Cancer", "sample_id"), drop = FALSE]
pan_data$panel_order <- 12L
pan_data$panel_id <- "KM12_Pan-cancer"
pan_data$display_label <- "Pan-cancer"
pan_data$source_scope <- "pan_cancer_complete_case_all_endpoints"
pan_data <- pan_data[, c("panel_order", "panel_id", "display_label", "Cancer", "source_scope",
                         "sample_id", "Group", "OS", "OS.time")]

plot_data <- rbind(within_data, pan_data)
plot_data$panel_order <- as.integer(plot_data$panel_order)
plot_data <- plot_data[order(plot_data$panel_order, plot_data$sample_id), , drop = FALSE]
rownames(plot_data) <- NULL

summary_rows <- lapply(seq_len(12L), function(i) {
  panel <- plot_data[plot_data$panel_order == i, , drop = FALSE]
  p <- get_logrank(panel)
  data.frame(
    panel_order = i,
    panel_id = unique(panel$panel_id),
    display_label = unique(panel$display_label),
    source_scope = unique(panel$source_scope),
    n = nrow(panel),
    n_low = sum(panel$Group == "Low"),
    n_high = sum(panel$Group == "High"),
    events = sum(panel$OS == 1),
    logrank_p = p,
    significance_rule = if (i < 12L) "original KM.R log-rank p < 0.05" else "descriptive pan-cancer OS log-rank",
    significant = if (i < 12L) p < 0.05 else NA,
    complete_case_rule = if (i < 12L) "Group + OS + OS.time complete" else "Group + OS + OS.time + PFI + PFI.time + DSS + DSS.time complete",
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)

write.table(plot_data, file.path(path.expand(out_dir), "km_12panel_plot_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_df, file.path(path.expand(out_dir), "km_12panel_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_df[, c("panel_order", "panel_id", "display_label")],
            file.path(path.expand(out_dir), "km_12panel_panel_order.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(path.expand(out_dir), "sessionInfo.txt"))

qa <- data.frame(check_id = character(), value = character(), pass = logical(), stringsAsFactors = FALSE)
add_qa <- function(check_id, value, pass) {
  qa <<- rbind(qa, data.frame(check_id = check_id, value = as.character(value), pass = isTRUE(pass), stringsAsFactors = FALSE))
}

add_qa("input_exists", file.exists(path.expand(input_rdata)), file.exists(path.expand(input_rdata)))
add_qa("panel_count", paste(sort(unique(plot_data$panel_order)), collapse = ","), identical(sort(unique(plot_data$panel_order)), 1:12))
add_qa("within_cancer_order", paste(unique(plot_data$display_label[plot_data$panel_order < 12]), collapse = ","), identical(unique(plot_data$display_label[plot_data$panel_order < 12]), expected_cancers))
add_qa("all_groups_valid", paste(sort(unique(as.character(plot_data$Group))), collapse = ","), all(as.character(plot_data$Group) %in% c("Low", "High")))
add_qa("plot_data_no_missing", sum(!complete.cases(plot_data[, c("panel_order", "panel_id", "Group", "OS", "OS.time")])), all(complete.cases(plot_data[, c("panel_order", "panel_id", "Group", "OS", "OS.time")])))

within_summary <- summary_df[summary_df$panel_order < 12, , drop = FALSE]
within_match <- all(within_summary$n == unname(expected_n[within_summary$display_label])) &&
  all(within_summary$n_low == unname(expected_low[within_summary$display_label])) &&
  all(within_summary$n_high == unname(expected_high[within_summary$display_label])) &&
  all(within_summary$events == unname(expected_events[within_summary$display_label]))
add_qa("within_counts", paste(within_summary$n, collapse = ","), within_match)
add_qa("within_significant_count", sum(within_summary$significant), sum(within_summary$significant) == 11L)
add_qa("within_logrank_p_reproducible", paste(format(within_summary$logrank_p, scientific = TRUE, digits = 6), collapse = ","), all(abs(within_summary$logrank_p - unname(expected_p[within_summary$display_label])) < 1e-07))

pan_summary <- summary_df[summary_df$panel_order == 12, , drop = FALSE]
pan_match <- pan_summary$n == 5188L && pan_summary$n_low == 2639L && pan_summary$n_high == 2549L &&
  pan_summary$events == 1508L && abs(pan_summary$logrank_p - 1.165430382e-04) < 1e-10
add_qa("pan_cancer_counts_and_p", paste(pan_summary$n, pan_summary$n_low, pan_summary$n_high, pan_summary$events, format(pan_summary$logrank_p, scientific = TRUE), sep = ";"), pan_match)

out_files <- c("km_12panel_plot_data.tsv", "km_12panel_summary.tsv", "km_12panel_panel_order.tsv", "sessionInfo.txt")
nonempty <- all(file.info(file.path(path.expand(out_dir), out_files))$size > 0)
add_qa("output_files_nonempty", paste(out_files, collapse = ","), nonempty)
add_qa("all_pass", "see pass column", all(qa$pass))
write.table(qa, file.path(path.expand(qa_dir), "km_12panel_machine_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sprintf("all_qa_pass=%s\noutput_dir=%s\nqa_dir=%s", all(qa$pass), path.expand(out_dir), path.expand(qa_dir)),
           file.path(path.expand(qa_dir), "km_12panel_qa_summary.txt"))

manifest <- data.frame(
  file = out_files,
  bytes = file.info(file.path(path.expand(out_dir), out_files))$size,
  md5 = unname(tools::md5sum(file.path(path.expand(out_dir), out_files))),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(path.expand(qa_dir), "km_12panel_output_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!all(qa$pass)) stop("KM 12-panel machine QA failed; inspect ", path.expand(qa_dir))
message("KM 12-panel figure-ready export and QA completed: ", path.expand(out_dir))
