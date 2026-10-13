#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 05_integrate_three_cohort.R --config PATH [--check-config]\n       Rscript 05_integrate_three_cohort.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("riaz_gide_preparation_dir","riaz_gide_analysis_dir","braun_preparation_dir","braun_analysis_dir","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
cfg_paths <- paths
# 固定文件名路由只读取各阶段输入目录；不读取表内容。
icb_input <- function(group, name, optional = FALSE) {
  if (!(length(group)==1L && group %in% c("riaz_gide","braun"))) stop("Unknown input group.")
  prep <- if(group=="riaz_gide") c("B_baseline_patient_inclusion_and_scores.tsv","B_member_coverage.tsv") else c("B_Braun_patient_inclusion_and_scores.tsv","B_Braun_member_coverage.tsv")
  analysis <- if(group=="riaz_gide") c("B_four_state_response_rates.tsv","B_primary_exact_comparisons.tsv","B_continuous_Firth_auxiliary.tsv","B_C1_internal_programme_exploration.tsv","B_CRPR_vs_PD_sensitivity.tsv","B_clinical_Firth_sensitivity.tsv") else c("B_Braun_four_state_response_rates.tsv","B_Braun_primary_exact_comparisons.tsv","B_Braun_continuous_Firth_auxiliary.tsv","B_Braun_C1_internal_programme_exploration.tsv","B_Braun_CRPR_vs_PD_sensitivity.tsv","B_Braun_clinical_Firth_sensitivity.tsv")
  extras <- if(group=="braun") c("B_Braun_expression_scale_and_parameters.tsv","B_Braun_duplicate_gene_handling.tsv","B_Braun_source_resolution_summary.tsv") else character()
  allowed <- if(optional) extras else c(prep,analysis)
  if(anyNA(name) || any(!name %in% allowed)) stop("Unknown fixed input basename.")
  key <- paste0(group, if(optional || all(name %in% prep)) "_preparation_dir" else "_analysis_dir")
  vapply(name, function(nm) {
    selected_key <- if(optional || nm %in% prep) paste0(group,"_preparation_dir") else paste0(group,"_analysis_dir")
    base <- unname(cfg_paths[selected_key])
    path <- resolve_path(file.path(base,nm),base)
    if(!is_within(path,base)) stop("Input symlink escapes configured directory.")
    if(file.exists(path) && !file_test("-f",path)) stop("Input is not a regular file.")
    if(!optional && (!file_test("-f",path) || file.access(path,4L)!=0L || is.na(file.info(path)$size) || file.info(path)$size<=0)) stop("Required input must be readable nonempty regular file.")
    if(optional && file.exists(path) && file.access(path,4L)!=0L) stop("Optional input is unreadable.")
    path
  }, character(1L), USE.NAMES=FALSE)
}

input_dirs <- unname(paths[c("riaz_gide_preparation_dir","riaz_gide_analysis_dir","braun_preparation_dir","braun_analysis_dir")])
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
required_files <- c(icb_input("riaz_gide",c("B_baseline_patient_inclusion_and_scores.tsv","B_member_coverage.tsv","B_four_state_response_rates.tsv","B_primary_exact_comparisons.tsv","B_continuous_Firth_auxiliary.tsv","B_C1_internal_programme_exploration.tsv","B_CRPR_vs_PD_sensitivity.tsv","B_clinical_Firth_sensitivity.tsv")),icb_input("braun",c("B_Braun_patient_inclusion_and_scores.tsv","B_Braun_member_coverage.tsv","B_Braun_four_state_response_rates.tsv","B_Braun_primary_exact_comparisons.tsv","B_Braun_continuous_Firth_auxiliary.tsv","B_Braun_C1_internal_programme_exploration.tsv","B_Braun_CRPR_vs_PD_sensitivity.tsv","B_Braun_clinical_Firth_sensitivity.tsv")))
invisible(icb_input("braun",c("B_Braun_expression_scale_and_parameters.tsv","B_Braun_duplicate_gene_handling.tsv","B_Braun_source_resolution_summary.tsv"),optional=TRUE))

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)

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
if(file.exists(cfg_out_dir) || dir.exists(cfg_out_dir)) stop("Integration output must be absent.")
if(!dir.exists(dirname(cfg_out_dir))) stop("Integration output parent must already exist.")
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}

cfg_paths <- paths

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))

# Private project root replaced by explicit DCF inputs.
# Riaz/Gide inputs route to their separate preparation/analysis directories.
# Braun inputs route to their separate preparation/analysis directories.
out_dir <- unname(cfg_paths["output_dir"])
if (dir.exists(out_dir)) stop("Versioned output directory already exists: ", out_dir)

required <- c(
  "B_Braun_patient_inclusion_and_scores.tsv",
  "B_Braun_four_state_response_rates.tsv",
  "B_Braun_primary_exact_comparisons.tsv",
  "B_Braun_continuous_Firth_auxiliary.tsv",
  "B_Braun_C1_internal_programme_exploration.tsv",
  "B_Braun_member_coverage.tsv",
  "B_Braun_CRPR_vs_PD_sensitivity.tsv",
  "B_Braun_clinical_Firth_sensitivity.tsv"
)
missing <- required[!file.exists(icb_input("braun", required))]
if (length(missing)) stop("Braun outputs are not ready: ", paste(missing, collapse=", "))
if (any(file.info(icb_input("braun", required))$size <= 0)) stop("At least one required Braun output is empty.")

read_table <- function(path) fread(path, na.strings=c("NA",""))
assert_braun <- function(x, label) {
  if (!"cohort" %in% names(x) || any(x$cohort != "Braun")) stop(label, " must contain Braun rows only")
  invisible(x)
}
combine_replace_braun <- function(old_name, braun_name, out_name) {
  old <- read_table(icb_input("riaz_gide", old_name))
  new <- read_table(icb_input("braun", braun_name))
  assert_braun(new, braun_name)
  out <- rbindlist(list(old[cohort != "Braun"], new), use.names=TRUE, fill=TRUE)
  if ("contrast" %in% names(out)) {
    out[, cohort_order:=match(cohort,c("Riaz","Gide","Braun"))]
    out[, contrast_order:=match(contrast,c("C1/C4","C1/C3"))]
    setorderv(out,c("cohort_order","contrast_order"),na.last=TRUE)
    out[,c("cohort_order","contrast_order"):=NULL]
  } else if ("C_state" %in% names(out)) {
    out[, cohort_order:=match(cohort,c("Riaz","Gide","Braun"))]
    out[, state_order:=match(C_state,c("C1","C2","C3","C4"))]
    setorderv(out,c("cohort_order","state_order"),na.last=TRUE)
    out[,c("cohort_order","state_order"):=NULL]
  } else if ("axis" %in% names(out)) {
    out[, cohort_order:=match(cohort,c("Riaz","Gide","Braun"))]
    out[, axis_order:=match(axis,c("Acid_Score","LM12_inflammatory_axis"))]
    setorderv(out,c("cohort_order","axis_order"),na.last=TRUE)
    out[,c("cohort_order","axis_order"):=NULL]
  } else if ("program" %in% names(out)) {
    out[, cohort_order:=match(cohort,c("Riaz","Gide","Braun"))]
    out[, program_order:=match(program,c("Hypoxia","ECM","Cytotoxic"))]
    setorderv(out,c("cohort_order","program_order"),na.last=TRUE)
    out[,c("cohort_order","program_order"):=NULL]
  } else {
    setorder(out,cohort)
  }
  fwrite(out,file.path(out_dir,out_name),sep="\t",na="NA")
  out
}
check_grid <- function(d, fields, expected_n, label) {
  if (nrow(d) != expected_n || anyDuplicated(d[, ..fields])) stop(label, " does not contain the fixed complete slot grid")
}
recompute_bh <- function(d, pcol, qcol, family_n, label) {
  if (!pcol %in% names(d)) stop(label, " is missing ", pcol)
  d[, (qcol):=NA_real_]
  ok <- which(is.finite(d[[pcol]]))
  if (length(ok)) d[[qcol]][ok] <- p.adjust(d[[pcol]][ok],method="BH",n=family_n)
  d
}

if (!dir.create(out_dir, recursive=FALSE)) stop("Could not create versioned output directory")
# 合并患者级结果和覆盖记录，供完整复核与作图读取。
patients_old <- read_table(icb_input("riaz_gide","B_baseline_patient_inclusion_and_scores.tsv"))
patients_new <- read_table(icb_input("braun","B_Braun_patient_inclusion_and_scores.tsv"))
assert_braun(patients_new,"B_Braun_patient_inclusion_and_scores.tsv")
patients <- rbindlist(list(patients_old[cohort != "Braun"],patients_new),use.names=TRUE,fill=TRUE)
if (anyDuplicated(patients[,.(cohort,sample_id)])) stop("Duplicate cohort/sample_id in combined patient table")
patients[,cohort_order:=match(cohort,c("Riaz","Gide","Braun"))]
setorderv(patients,c("cohort_order","sample_id"),na.last=TRUE)
patients[,cohort_order:=NULL]
fwrite(patients,file.path(out_dir,"B_baseline_patient_inclusion_and_scores.tsv"),sep="\t",na="NA")

coverage_old <- read_table(icb_input("riaz_gide","B_member_coverage.tsv"))
coverage_new <- read_table(icb_input("braun","B_Braun_member_coverage.tsv"))
assert_braun(coverage_new,"B_Braun_member_coverage.tsv")
coverage <- rbindlist(list(coverage_old[cohort != "Braun"],coverage_new),use.names=TRUE,fill=TRUE)
fwrite(coverage,file.path(out_dir,"B_member_coverage.tsv"),sep="\t",na="NA")

rates <- combine_replace_braun("B_four_state_response_rates.tsv","B_Braun_four_state_response_rates.tsv","B_four_state_response_rates.tsv")
check_grid(rates,c("cohort","C_state"),12L,"Four-state response table")
if (!all(c("n_baseline","n_strict","R","NR","response_rate","ci_lo","ci_hi","status") %in% names(rates))) stop("Response-rate schema is incomplete")

primary <- combine_replace_braun("B_primary_exact_comparisons.tsv","B_Braun_primary_exact_comparisons.tsv","B_primary_exact_comparisons.tsv")
check_grid(primary,c("cohort","contrast"),6L,"Primary exact comparison table")
primary <- recompute_bh(primary,"p","q",6L,"Primary exact comparison table")
fwrite(primary,file.path(out_dir,"B_primary_exact_comparisons.tsv"),sep="\t",na="NA")

continuous <- combine_replace_braun("B_continuous_Firth_auxiliary.tsv","B_Braun_continuous_Firth_auxiliary.tsv","B_continuous_Firth_auxiliary.tsv")
check_grid(continuous,c("cohort","axis"),6L,"Continuous auxiliary table")
continuous <- recompute_bh(continuous,"p","q",6L,"Continuous auxiliary table")
fwrite(continuous,file.path(out_dir,"B_continuous_Firth_auxiliary.tsv"),sep="\t",na="NA")

internal <- combine_replace_braun("B_C1_internal_programme_exploration.tsv","B_Braun_C1_internal_programme_exploration.tsv","B_C1_internal_programme_exploration.tsv")
check_grid(internal,c("cohort","program"),9L,"C1 internal exploration table")
internal <- recompute_bh(internal,"p","q",9L,"C1 internal exploration table")
fwrite(internal,file.path(out_dir,"B_C1_internal_programme_exploration.tsv"),sep="\t",na="NA")

crpr <- combine_replace_braun("B_CRPR_vs_PD_sensitivity.tsv","B_Braun_CRPR_vs_PD_sensitivity.tsv","B_CRPR_vs_PD_sensitivity.tsv")
check_grid(crpr,c("cohort","contrast"),6L,"CR/PR versus PD sensitivity table")
clinical <- combine_replace_braun("B_clinical_Firth_sensitivity.tsv","B_Braun_clinical_Firth_sensitivity.tsv","B_clinical_Firth_sensitivity.tsv")

# 保留 Braun 来源侧的原始参数记录，避免把尺度证据压缩到汇总统计表。
for (nm in c("B_Braun_expression_scale_and_parameters.tsv","B_Braun_duplicate_gene_handling.tsv","B_Braun_source_resolution_summary.tsv")) {
  src <- icb_input("braun",nm,optional=TRUE)
  if (file.exists(src) && file.info(src)$size > 0) file.copy(src,file.path(out_dir,nm),overwrite=FALSE)
}

cohort_n <- patients[,.(n_baseline=.N,n_response=sum(!is.na(strict_response))),by=cohort]
cat("Combined patient counts (all assigned baselines / strict response):\n")
print(cohort_n)
cat(sprintf("BH families: primary=%d slots, continuous=%d slots, C1 internal=%d slots\n",nrow(primary),nrow(continuous),nrow(internal)))
cat(sprintf("Estimated: primary=%d, continuous=%d, C1 internal=%d\n",
  sum(is.finite(primary$p)),sum(is.finite(continuous$p)),sum(is.finite(internal$p))))
cat("Versioned combined tables written to ",out_dir,"\n",sep="")

