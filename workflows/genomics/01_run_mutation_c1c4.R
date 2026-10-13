#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 01_run_mutation_c1c4.R --config PATH [--check-config]\n       Rscript 01_run_mutation_c1c4.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("callable_dir", "mutation_events", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
callable_dir <- unname(paths["callable_dir"])
if(!dir.exists(callable_dir) || file.access(callable_dir,5L)!=0L) stop("Callable input must be a readable searchable directory.")
callable_names <- c("mutation_callable_decision_v0.1.tsv", "mutation_callable_recovery_qa.tsv",
 "mutation_callable_patient_audit.tsv", "mutation_callable_formal_eligibility_by_cancer.tsv")
required_files <- c(file.path(callable_dir, callable_names), unname(paths["mutation_events"]))
required_files <- vapply(required_files,resolve_path,character(1L),base_dir=dirname(config_path))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(callable_dir,dirname(required_files),config_path,package_root)
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

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))
xena_path <- unname(paths["mutation_events"])

# 仅普通执行进入保留的分析主体；输入与新输出隔离。
log_msg <- function(level, step, ...) {
  detail <- paste(..., collapse = " ")
  cat(sprintf(
    "%s\t%s\t%s\t%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level, step, detail
  ))
  flush.console()
}

assert_nonempty <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }
  size <- file.info(path)$size
  if (is.na(size) || size <= 0) {
    stop("Empty input: ", path)
  }
  invisible(size)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

save_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), path)
}

nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Translation_Start_Site",
  "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"
)

fisher_gene <- function(a, n1, c, n4) {
  ft <- fisher.test(
    matrix(c(a, n1 - a, c, n4 - c), nrow = 2L, byrow = TRUE),
    alternative = "two.sided"
  )
  c(
    p_value = unname(ft$p.value),
    odds_ratio_C1_vs_C4 = if (length(ft$estimate)) unname(ft$estimate) else NA_real_,
    odds_ratio_ci_low = unname(ft$conf.int[1L]),
    odds_ratio_ci_high = unname(ft$conf.int[2L])
  )
}

tryCatch({
  log_msg("START", "mutation_c1c4", "Formal within-cancer C1 vs C4 mutation analysis")
  required <- c(
    file.path(callable_dir, "mutation_callable_decision_v0.1.tsv"),
    file.path(callable_dir, "mutation_callable_recovery_qa.tsv"),
    file.path(callable_dir, "mutation_callable_patient_audit.tsv"),
    file.path(callable_dir, "mutation_callable_formal_eligibility_by_cancer.tsv"),
    xena_path
  )
  invisible(lapply(required, assert_nonempty))

  decision <- fread(required[1L])
  recovery_qa <- fread(required[2L])
  if (!all(recovery_qa$pass)) stop("Callable recovery QA did not pass")
  if (!isTRUE(decision$accepted_for_formal_denominator[1L])) {
    stop("Callable universe was not accepted for formal mutation analysis")
  }

  callable_audit <- fread(required[3L])
  callable <- callable_audit[
    mc3_callable_exact_sample15 == TRUE,
    .(Sample15 = selected_canonical_sample15, Patient_ID, Cancer, C_subtype)
  ]
  if (anyDuplicated(callable$Patient_ID)) stop("Duplicate callable patients")
  if (any(!callable$C_subtype %chin% c("C1", "C4"))) stop("Non-C1/C4 subtype")

  write_tsv(callable, file.path(out_dir, "mutation_callable_patient_manifest.tsv"))

  log_msg("CHECKPOINT", "event_read", "Reading frozen Xena PASS mutation events")
  events <- fread(xena_path, select = c("sample", "gene", "effect"))
  events[, Sample15 := substr(as.character(sample), 1L, 15L)]
  events <- events[
    effect %chin% nonsilent_classes &
      !is.na(gene) & gene != "" &
      Sample15 %chin% callable$Sample15,
    .(Sample15, Gene = as.character(gene), Variant_Classification = effect)
  ]
  event_class_counts <- events[, .(
    event_rows = .N,
    unique_sample_gene_events = uniqueN(paste(Sample15, Gene, sep = "|")),
    unique_samples = uniqueN(Sample15),
    unique_genes = uniqueN(Gene)
  ), by = Variant_Classification][order(-event_rows)]
  write_tsv(
    event_class_counts,
    file.path(out_dir, "mutation_nonsilent_event_class_counts.tsv")
  )

  sparse_binary <- unique(events[, .(Sample15, Gene)])
  sparse_binary[, mutated := 1L]
  sparse_binary <- merge(
    sparse_binary,
    callable,
    by = "Sample15",
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(sparse_binary$Patient_ID))) {
    stop("Mutation event outside accepted callable manifest")
  }
  setorder(sparse_binary, Cancer, Sample15, Gene)
  write_tsv(
    sparse_binary[, .(
      Sample15, Patient_ID, Cancer, C_subtype, Gene, mutated
    )],
    file.path(out_dir, "mutation_nonsilent_patient_gene_sparse_binary.tsv")
  )
  rm(events)
  gc()

  cancer_coverage <- callable[, .N, by = .(Cancer, C_subtype)]
  cancer_coverage <- dcast(
    cancer_coverage,
    Cancer ~ C_subtype,
    value.var = "N",
    fill = 0L
  )
  if (!"C1" %chin% names(cancer_coverage)) cancer_coverage[, C1 := 0L]
  if (!"C4" %chin% names(cancer_coverage)) cancer_coverage[, C4 := 0L]
  cancer_coverage[, `:=`(
    min_group_n = 30L,
    formal_eligible = C1 >= 30L & C4 >= 30L
  )]
  setorder(cancer_coverage, Cancer)
  write_tsv(
    cancer_coverage,
    file.path(out_dir, "mutation_formal_coverage_by_cancer.tsv")
  )

  eligible_cancers <- cancer_coverage[formal_eligible == TRUE, Cancer]
  result_list <- vector("list", length(eligible_cancers))

  for (i in seq_along(eligible_cancers)) {
    cancer <- eligible_cancers[i]
    group_dat <- callable[Cancer == cancer]
    n1 <- group_dat[C_subtype == "C1", .N]
    n4 <- group_dat[C_subtype == "C4", .N]

    event_dat <- sparse_binary[Cancer == cancer]
    counts <- event_dat[, .(
      mutated_patients = uniqueN(Sample15)
    ), by = .(Gene, C_subtype)]
    counts_wide <- dcast(
      counts,
      Gene ~ C_subtype,
      value.var = "mutated_patients",
      fill = 0L
    )
    if (!"C1" %chin% names(counts_wide)) counts_wide[, C1 := 0L]
    if (!"C4" %chin% names(counts_wide)) counts_wide[, C4 := 0L]
    counts_wide <- counts_wide[C1 + C4 >= 5L]

    if (nrow(counts_wide) == 0L) {
      log_msg("CHECKPOINT", "cancer", paste(cancer, "no gene passed event threshold"))
      next
    }

    test_values <- t(vapply(
      seq_len(nrow(counts_wide)),
      function(j) fisher_gene(counts_wide$C1[j], n1, counts_wide$C4[j], n4),
      numeric(4L)
    ))
    results <- cbind(
      data.table(
        Cancer = cancer,
        Gene = counts_wide$Gene,
        n_C1 = n1,
        n_C4 = n4,
        mutated_C1 = counts_wide$C1,
        mutated_C4 = counts_wide$C4
      ),
      as.data.table(test_values)
    )
    results[, `:=`(
      frequency_C1 = mutated_C1 / n_C1,
      frequency_C4 = mutated_C4 / n_C4
    )]
    results[, frequency_difference_C1_minus_C4 := frequency_C1 - frequency_C4]
    results[, abs_frequency_difference := abs(frequency_difference_C1_minus_C4)]
    results[, FDR_within_cancer := p.adjust(p_value, method = "BH")]
    results[, single_cancer_significant :=
      FDR_within_cancer < 0.05 &
        abs(frequency_difference_C1_minus_C4) >= 0.10
    ]
    results[, direction := fifelse(
      frequency_difference_C1_minus_C4 > 0,
      "C1_higher",
      fifelse(frequency_difference_C1_minus_C4 < 0, "C4_higher", "tie")
    )]
    setorder(results, FDR_within_cancer, -abs_frequency_difference)
    result_list[[i]] <- results
    log_msg(
      "CHECKPOINT", "cancer",
      paste(
        cancer,
        "n_C1:", n1,
        "n_C4:", n4,
        "tested genes:", nrow(results),
        "significant:", sum(results$single_cancer_significant)
      )
    )
  }

  mutation_results <- rbindlist(result_list, fill = TRUE)
  if (nrow(mutation_results) == 0L) stop("No mutation gene tests were generated")
  write_tsv(
    mutation_results,
    file.path(out_dir, "mutation_gene_c1_vs_c4_by_cancer.tsv")
  )

  cross_cancer <- mutation_results[, .(
    n_tested_cancers = uniqueN(Cancer),
    n_significant_cancers = sum(single_cancer_significant),
    n_significant_C1_higher = sum(
      single_cancer_significant & direction == "C1_higher"
    ),
    n_significant_C4_higher = sum(
      single_cancer_significant & direction == "C4_higher"
    ),
    significant_cancers_C1_higher = paste(
      sort(Cancer[single_cancer_significant & direction == "C1_higher"]),
      collapse = ";"
    ),
    significant_cancers_C4_higher = paste(
      sort(Cancer[single_cancer_significant & direction == "C4_higher"]),
      collapse = ";"
    ),
    max_abs_frequency_difference = max(
      abs(frequency_difference_C1_minus_C4),
      na.rm = TRUE
    ),
    min_within_cancer_FDR = min(FDR_within_cancer, na.rm = TRUE)
  ), by = Gene]
  cross_cancer[, context_dependent := n_significant_C1_higher > 0L &
    n_significant_C4_higher > 0L
  ]
  cross_cancer[, dominant_significant_count := pmax(
    n_significant_C1_higher,
    n_significant_C4_higher
  )]
  setorder(
    cross_cancer,
    -n_significant_cancers,
    -dominant_significant_count,
    min_within_cancer_FDR,
    Gene
  )
  write_tsv(
    cross_cancer,
    file.path(out_dir, "mutation_gene_cross_cancer_significant_counts.tsv")
  )
  write_tsv(
    cross_cancer[n_significant_cancers > 0L],
    file.path(out_dir, "mutation_single_cancer_significant_gene_summary.tsv")
  )

  significant_units <- mutation_results[single_cancer_significant == TRUE]
  write_tsv(
    significant_units,
    file.path(out_dir, "mutation_significant_gene_cancer_units.tsv")
  )

  endpoint_summary <- data.table(
    endpoint = "nonsilent_gene_mutation",
    n_callable_patients = nrow(callable),
    n_formal_cancers = length(eligible_cancers),
    n_gene_cancer_tests = nrow(mutation_results),
    n_significant_gene_cancer_units = nrow(significant_units),
    n_genes_significant_in_at_least_one_cancer =
      cross_cancer[n_significant_cancers > 0L, .N],
    n_genes_recurrent_same_direction_ge2_cancers = cross_cancer[
      pmax(n_significant_C1_higher, n_significant_C4_higher) >= 2L,
      .N
    ],
    n_context_dependent_genes = cross_cancer[context_dependent == TRUE, .N]
  )
  write_tsv(
    endpoint_summary,
    file.path(out_dir, "mutation_endpoint_summary.tsv")
  )

  qa <- data.table(
    check = c(
      "callable_recovery_qa_pass",
      "only_exact_callable_samples_used",
      "callable_patient_unique",
      "only_c1_c4",
      "formal_cancers_meet_group_n_30",
      "all_tested_genes_meet_total_mutated_n_5",
      "bh_is_within_cancer",
      "effect_threshold_applied_to_significance",
      "no_cross_cancer_bh",
      "no_figure_output_generated"
    ),
    pass = c(
      all(recovery_qa$pass),
      all(sparse_binary$Sample15 %chin% callable$Sample15),
      !anyDuplicated(callable$Patient_ID),
      all(callable$C_subtype %chin% c("C1", "C4")),
      cancer_coverage[formal_eligible == TRUE, all(C1 >= 30L & C4 >= 30L)],
      mutation_results[, all(mutated_C1 + mutated_C4 >= 5L)],
      mutation_results[, all(
        abs(
          FDR_within_cancer -
            ave(p_value, Cancer, FUN = function(z) p.adjust(z, method = "BH"))
        ) < 1e-12
      )],
      mutation_results[single_cancer_significant == TRUE, all(
        FDR_within_cancer < 0.05 &
          abs(frequency_difference_C1_minus_C4) >= 0.10
      )],
      !any(c("cross_cancer_FDR", "meta_FDR") %chin% names(cross_cancer)),
      !dir.exists(file.path(out_dir, "figures/01_mutation"))
    )
  )
  write_tsv(qa, file.path(out_dir, "mutation_analysis_qa.tsv"))
  if (!all(qa$pass)) {
    stop(
      "Mutation analysis QA failed: ",
      paste(qa[pass == FALSE, check], collapse = ", ")
    )
  }

  save_session_info(file.path(out_dir, "mutation_analysis_sessionInfo.txt"))
  log_msg(
    "DONE", "mutation_c1c4",
    paste(
      "Callable patients:", nrow(callable),
      "formal cancers:", length(eligible_cancers),
      "gene-cancer tests:", nrow(mutation_results),
      "significant units:", nrow(significant_units)
    )
  )
}, error = function(e) {
  log_msg("ERROR", "mutation_c1c4", conditionMessage(e))
  quit(status = 1L)
})
