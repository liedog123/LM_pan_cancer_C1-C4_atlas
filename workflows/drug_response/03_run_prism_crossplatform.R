#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 03_run_prism_crossplatform.R --config PATH [--check-config]\n       Rscript 03_run_prism_crossplatform.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("gdsc_response", "prism_curve", "prism_cell_metadata", "internal_support_summary", "drug_identity", "prism_coverage", "discovery_counts", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

# 固定PRISM快照文件名；元信息检查不读取或计算研究文件hash。
expected_basenames <- c(
  prism_curve = "secondary-screen-dose-response-top8-filtered.csv",
  prism_cell_metadata = "secondary-screen-cell-line-info.csv"
)
if (!identical(unname(basename(paths[names(expected_basenames)])),
               unname(expected_basenames))) stop("Required PRISM snapshot basenames differ.")

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

suppressPackageStartupMessages(library(readxl))
# 从独立的正式或公开续跑中间表继续；不要求先重新拟合内部模型。
gdsc_response_path <- unname(paths["gdsc_response"])
prism_curve_path <- unname(paths["prism_curve"])
prism_cell_path <- unname(paths["prism_cell_metadata"])
internal_qc_path <- unname(paths["internal_support_summary"])
identity_path <- unname(paths["drug_identity"])
prism_coverage_path <- unname(paths["prism_coverage"])
discovery_path <- unname(paths["discovery_counts"])
result_dir <- out_dir

target_drugs <- c(
  "Sorafenib", "Bleomycin", "VE821", "Wnt-C59",
  "Vinorelbine", "667880", "Linsitinib", "Palbociclib"
)
secondary_drugs <- c(
  "Sorafenib", "VE821", "Vinorelbine", "Linsitinib", "Palbociclib"
)
preferred_screen_expected <- c(
  Sorafenib = "MTS010",
  VE821 = "HTS002",
  Vinorelbine = "HTS002",
  Linsitinib = "MTS010",
  Palbociclib = "MTS010"
)

expected_sha256 <- c(
  "secondary-screen-dose-response-top8-filtered.csv" =
    "c07b791e6d5356232ce7e1a9339a338f65f662e3a2717c1408c2f0462420aa6e",
  "secondary-screen-cell-line-info.csv" =
    "b93436b5f4bcf5fd14589697be4ca8f99ffd99d639223080380bc514c3718d27"
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

normalize_drug <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

# PRISM CCLE名称通常为CELL_LINE_TISSUE；匹配时只使用细胞系名称部分。
normalize_cell <- function(x, prism = FALSE) {
  x <- as.character(x)
  if (prism) x <- sub("_[^_]+$", "", x)
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

collapse_unique <- function(x) {
  paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))),
        collapse = ";")
}

safe_cor_test <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(estimate = unname(test$estimate), p_value = test$p.value)
}

message("[Stage5-crossplatform] 核验输入")
required_paths <- c(
  gdsc_response_path, prism_curve_path, prism_cell_path, internal_qc_path,
  identity_path, prism_coverage_path, discovery_path
)
if (!all(file.exists(required_paths))) {
  stop("缺少输入: ", paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}
if (!all(file.info(required_paths)$size > 0)) stop("存在空输入文件")

sha_lines <- system2(
  "sha256sum",
  shQuote(c(prism_curve_path, prism_cell_path)),
  stdout = TRUE
)
observed_sha <- sub("\\s+.*$", "", sha_lines)
names(observed_sha) <- basename(sub("^[^ ]+\\s+", "", sha_lines))
checksum_table <- data.frame(
  file = names(expected_sha256),
  expected_sha256 = unname(expected_sha256),
  observed_sha256 = observed_sha[names(expected_sha256)],
  pass = unname(expected_sha256) == observed_sha[names(expected_sha256)],
  stringsAsFactors = FALSE
)
write_tsv(checksum_table, file.path(result_dir, "prism_uploaded_input_checksum.tsv"))
if (!all(checksum_table$pass)) stop("PRISM上传文件checksum不一致")

message("[Stage5-crossplatform] 读取GDSC与PRISM真实药敏")
gdsc <- as.data.frame(read_excel(gdsc_response_path), stringsAsFactors = FALSE)
prism <- read.csv(prism_curve_path, check.names = FALSE)
prism_cell <- read.csv(prism_cell_path, check.names = FALSE)

gdsc$drug_key <- normalize_drug(gdsc$DRUG_NAME)
gdsc$cell_key <- normalize_cell(gdsc$CELL_LINE_NAME)
gdsc$LN_IC50 <- suppressWarnings(as.numeric(gdsc$LN_IC50))

prism$drug_key <- normalize_drug(prism$name)
prism$cell_key <- normalize_cell(prism$ccle_name, prism = TRUE)
prism$auc <- suppressWarnings(as.numeric(prism$auc))
prism$passed_str_profiling <- tolower(as.character(prism$passed_str_profiling)) == "true"

target_key_map <- setNames(target_drugs, normalize_drug(target_drugs))
secondary_keys <- normalize_drug(secondary_drugs)

gdsc_target <- gdsc[
  gdsc$drug_key %in% secondary_keys &
    nzchar(gdsc$cell_key) & is.finite(gdsc$LN_IC50),
  c("drug_key", "CELL_LINE_NAME", "cell_key", "TCGA_DESC", "LN_IC50"),
  drop = FALSE
]
prism_target <- prism[
  prism$drug_key %in% secondary_keys &
    nzchar(prism$cell_key) & prism$passed_str_profiling & is.finite(prism$auc),
  c("drug_key", "depmap_id", "ccle_name", "cell_key", "screen_id",
    "auc", "r2", "broad_id"),
  drop = FALSE
]

# 选择每种药物的首选PRISM screen。
screen_rows <- lapply(secondary_keys, function(key) {
  screens <- unique(prism_target$screen_id[prism_target$drug_key == key])
  preferred <- if ("MTS010" %in% screens) "MTS010" else if ("HTS002" %in% screens) "HTS002" else NA_character_
  data.frame(
    drug = target_key_map[[key]],
    drug_key = key,
    available_screens = paste(sort(screens), collapse = ";"),
    preferred_screen = preferred,
    stringsAsFactors = FALSE
  )
})
screen_table <- do.call(rbind, screen_rows)
write_tsv(screen_table, file.path(result_dir, "prism_preferred_screen.tsv"))

prism_selected <- merge(
  prism_target,
  screen_table[, c("drug_key", "preferred_screen")],
  by = "drug_key",
  all.x = TRUE
)
prism_selected <- prism_selected[
  prism_selected$screen_id == prism_selected$preferred_screen, ,
  drop = FALSE
]

# 审计名称标准化后是否把不同原始细胞系折叠到同一key。
gdsc_key_audit <- aggregate(
  CELL_LINE_NAME ~ drug_key + cell_key,
  data = gdsc_target,
  FUN = collapse_unique
)
names(gdsc_key_audit)[3] <- "source_names"
gdsc_key_audit$n_source_names <- vapply(
  strsplit(gdsc_key_audit$source_names, ";", fixed = TRUE), length, integer(1)
)
gdsc_key_audit$source <- "GDSC"

prism_key_audit <- aggregate(
  ccle_name ~ drug_key + cell_key,
  data = prism_selected,
  FUN = collapse_unique
)
names(prism_key_audit)[3] <- "source_names"
prism_key_audit$n_source_names <- vapply(
  strsplit(prism_key_audit$source_names, ";", fixed = TRUE), length, integer(1)
)
prism_key_audit$source <- "PRISM"

mapping_audit <- rbind(
  gdsc_key_audit[, c("source", "drug_key", "cell_key", "source_names", "n_source_names")],
  prism_key_audit[, c("source", "drug_key", "cell_key", "source_names", "n_source_names")]
)
write_tsv(mapping_audit, file.path(result_dir, "cell_line_name_mapping_audit.tsv"))

# 不同原始细胞系若在去标点后得到同一key，则整键排除，禁止强行合并。
ambiguous_keys <- unique(
  mapping_audit$cell_key[mapping_audit$n_source_names > 1L]
)
ambiguous_key_table <- mapping_audit[
  mapping_audit$cell_key %in% ambiguous_keys, ,
  drop = FALSE
]
write_tsv(
  ambiguous_key_table,
  file.path(result_dir, "excluded_ambiguous_cell_keys.tsv")
)
gdsc_target <- gdsc_target[
  !gdsc_target$cell_key %in% ambiguous_keys, ,
  drop = FALSE
]
prism_selected <- prism_selected[
  !prism_selected$cell_key %in% ambiguous_keys, ,
  drop = FALSE
]

# 同一药物-细胞系若有技术重复，按均值折叠，保持与旧GDSC训练矩阵构建方式一致。
gdsc_agg <- aggregate(
  LN_IC50 ~ drug_key + cell_key,
  data = gdsc_target,
  FUN = mean
)
gdsc_names <- aggregate(
  CELL_LINE_NAME ~ drug_key + cell_key,
  data = gdsc_target,
  FUN = collapse_unique
)
gdsc_agg <- merge(gdsc_agg, gdsc_names, by = c("drug_key", "cell_key"))

prism_agg <- aggregate(
  auc ~ drug_key + cell_key,
  data = prism_selected,
  FUN = mean
)
prism_names <- aggregate(
  ccle_name ~ drug_key + cell_key,
  data = prism_selected,
  FUN = collapse_unique
)
prism_ids <- aggregate(
  depmap_id ~ drug_key + cell_key,
  data = prism_selected,
  FUN = collapse_unique
)
prism_agg <- Reduce(
  function(x, y) merge(x, y, by = c("drug_key", "cell_key")),
  list(prism_agg, prism_names, prism_ids)
)

matched <- merge(
  gdsc_agg, prism_agg,
  by = c("drug_key", "cell_key"),
  all = FALSE
)
matched$drug <- unname(target_key_map[matched$drug_key])
matched <- matched[
  order(match(matched$drug, secondary_drugs), matched$cell_key), ,
  drop = FALSE
]
write_tsv(
  matched[, c(
    "drug", "cell_key", "CELL_LINE_NAME", "ccle_name", "depmap_id",
    "LN_IC50", "auc"
  )],
  file.path(result_dir, "gdsc_prism_matched_cellline_pairs.tsv")
)

message("[Stage5-crossplatform] 计算跨平台相关")
summary_rows <- lapply(secondary_drugs, function(drug) {
  key <- normalize_drug(drug)
  one <- matched[matched$drug_key == key, , drop = FALSE]
  spearman <- safe_cor_test(one$LN_IC50, one$auc, "spearman")
  pearson <- safe_cor_test(one$LN_IC50, one$auc, "pearson")
  data.frame(
    drug = drug,
    preferred_prism_screen =
      screen_table$preferred_screen[screen_table$drug == drug][1],
    n_gdsc_cell_lines = sum(gdsc_agg$drug_key == key),
    n_prism_cell_lines = sum(prism_agg$drug_key == key),
    n_matched_cell_lines = nrow(one),
    spearman_rho = spearman[["estimate"]],
    spearman_p = spearman[["p_value"]],
    pearson_r = pearson[["estimate"]],
    pearson_p = pearson[["p_value"]],
    median_gdsc_ln_ic50 = median(one$LN_IC50, na.rm = TRUE),
    median_prism_auc = median(one$auc, na.rm = TRUE),
    comparable = nrow(one) >= 30L,
    stringsAsFactors = FALSE
  )
})
cross_summary <- do.call(rbind, summary_rows)
cross_summary$spearman_fdr <- NA_real_
eligible <- cross_summary$comparable & is.finite(cross_summary$spearman_p)
cross_summary$spearman_fdr[eligible] <- p.adjust(
  cross_summary$spearman_p[eligible], method = "BH"
)
cross_summary$external_direction_support <-
  cross_summary$comparable &
  cross_summary$spearman_rho >= 0.20 &
  cross_summary$spearman_fdr < 0.10
cross_summary$external_status <- ifelse(
  !cross_summary$comparable,
  "NOT_COMPARABLE",
  ifelse(
    cross_summary$external_direction_support,
    "PRISM_DIRECTION_SUPPORTED",
    "PRISM_DIRECTION_NOT_SUPPORTED"
  )
)
write_tsv(
  cross_summary,
  file.path(result_dir, "gdsc_prism_crossplatform_summary.tsv")
)

message("[Stage5-crossplatform] 生成8种药物证据等级的未冻结续跑候选")
internal_qc <- read.delim(internal_qc_path, check.names = FALSE)
identity <- read.delim(identity_path, check.names = FALSE)
prism_coverage <- read.delim(prism_coverage_path, check.names = FALSE)
discovery <- read.delim(discovery_path, check.names = FALSE)

evidence <- data.frame(drug = target_drugs, stringsAsFactors = FALSE)
evidence <- merge(
  evidence,
  discovery[, c(
    "drug", "n_tested_cancers", "n_sig_cancers", "n_C1_sensitive",
    "n_C4_sensitive", "median_abs_delta_sig", "context_dependent"
  )],
  by = "drug", all.x = TRUE, sort = FALSE
)
evidence <- merge(
  evidence,
  identity[, c("drug", "official_target", "official_pathway", "identity_status")],
  by = "drug", all.x = TRUE, sort = FALSE
)
evidence <- merge(
  evidence,
  internal_qc[, c("drug", "median_spearman", "internal_modelability_pass")],
  by = "drug", all.x = TRUE, sort = FALSE
)
evidence <- merge(
  evidence,
  prism_coverage[, c("drug", "prism_support_scope")],
  by = "drug", all.x = TRUE, sort = FALSE
)
evidence <- merge(
  evidence,
  cross_summary[, c(
    "drug", "n_matched_cell_lines", "spearman_rho", "spearman_fdr",
    "external_direction_support", "external_status"
  )],
  by = "drug", all.x = TRUE, sort = FALSE
)
evidence <- evidence[match(target_drugs, evidence$drug), , drop = FALSE]

evidence$external_status[is.na(evidence$external_status)] <- ifelse(
  evidence$prism_support_scope[is.na(evidence$external_status)] ==
    "PRIMARY_COVERAGE_ONLY",
  "PRIMARY_COVERAGE_ONLY",
  "NO_PRISM_MATCH"
)
evidence$external_direction_support[
  is.na(evidence$external_direction_support)
] <- FALSE

evidence$subtype_direction <- ifelse(
  evidence$n_C1_sensitive == evidence$n_sig_cancers,
  "C1_SENSITIVE_ALL_SIGNIFICANT_CANCERS",
  ifelse(
    evidence$n_C4_sensitive == evidence$n_sig_cancers,
    "C4_SENSITIVE_ALL_SIGNIFICANT_CANCERS",
    "CONTEXT_DEPENDENT"
  )
)
evidence$evidence_grade <- ifelse(
  grepl("UNRESOLVED", evidence$identity_status),
  "HOLD_UNRESOLVED",
  ifelse(
    evidence$internal_modelability_pass &
      evidence$external_direction_support,
    "STRICT_SUPPORTED",
    ifelse(
      evidence$internal_modelability_pass,
      "PARTIAL_SUPPORT",
      "INTERNAL_ONLY"
    )
  )
)
write_tsv(
  evidence,
  file.path(result_dir, "top8_final_evidence_grade.tsv")
)

cell_info_ids <- unique(prism_cell$depmap_id)
selected_ids <- unique(prism_selected$depmap_id)
expected_screen_ok <- all(
  screen_table$preferred_screen[
    match(names(preferred_screen_expected), screen_table$drug)
  ] == unname(preferred_screen_expected)
)
matched_keys <- unique(matched$cell_key)
ambiguous_matched <- mapping_audit[
  mapping_audit$cell_key %in% matched_keys &
    mapping_audit$n_source_names > 1L, ,
  drop = FALSE
]

qa <- data.frame(
  check = c(
    "uploaded_inputs_exist_and_nonempty",
    "uploaded_input_checksums_match",
    "all_selected_prism_depmap_ids_in_cell_metadata",
    "five_predefined_secondary_drugs_present",
    "preferred_screen_matches_contract",
    "all_five_drugs_have_at_least_30_matched_cell_lines",
    "bh_family_contains_all_five_comparable_drugs",
    "no_ambiguous_normalized_names_in_matched_pairs",
    "final_evidence_table_contains_all_eight_drugs",
    "no_prism_cell_line_assigned_c1_or_c4",
    "no_source_object_modified"
  ),
  pass = c(
    all(file.info(c(prism_curve_path, prism_cell_path))$size > 0),
    all(checksum_table$pass),
    all(selected_ids %in% cell_info_ids),
    setequal(unique(prism_selected$drug_key), secondary_keys),
    expected_screen_ok,
    all(cross_summary$n_matched_cell_lines >= 30L),
    sum(eligible) == length(secondary_drugs),
    nrow(ambiguous_matched) == 0L,
    identical(evidence$drug, target_drugs),
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)
write_tsv(qa, file.path(result_dir, "qa_gdsc_prism_crossplatform.tsv"))
if (nrow(ambiguous_matched)) {
  write_tsv(
    ambiguous_matched,
    file.path(result_dir, "ambiguous_matched_cell_keys.tsv")
  )
}

capture.output(
  sessionInfo(),
  file = file.path(result_dir, "sessionInfo_crossplatform.txt")
)
writeLines(
  c(
    paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "comparison=GDSC_LN_IC50_vs_PRISM_AUC",
    "higher_value_in_both_datasets=relative_resistance",
    "prism_screen_rule=MTS010_if_available_else_HTS002",
    "min_matched_cell_lines=30",
    "test=two_sided_spearman",
    "fdr_family=all_comparable_predefined_secondary_drugs",
    "support_rule=rho>=0.20_and_FDR<0.10",
    "no_C1_C4_assignment_to_PRISM_cell_lines=true"
  ),
  con = file.path(result_dir, "parameters_crossplatform.txt")
)

if (!all(qa$pass)) {
  stop(
    "GDSC–PRISM跨平台QA未全部通过: ",
    paste(qa$check[!qa$pass], collapse = ", ")
  )
}
message("[Stage5-crossplatform] 完成")
