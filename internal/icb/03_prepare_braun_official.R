#!/usr/bin/env Rscript
# Purpose: Assign pretreatment states and evaluate ICB outcomes in Riaz, Gide and Braun.
# Task: icb_three_cohort; step: s03
# Inputs (DCF fields): official_workbook,lm12_members,program_membership,gene_alias,local_clinical_crosswalk
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task icb_three_cohort --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 03_prepare_braun_official.R --config PATH [--check-config]\n       Rscript 03_prepare_braun_official.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("official_workbook","lm12_members","program_membership","gene_alias","local_clinical_crosswalk","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("official_workbook","lm12_members","program_membership","gene_alias","local_clinical_crosswalk")]
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

cfg_paths <- paths
# 普通执行才建立独立空输出目录。
dir.create(unname(cfg_paths["output_dir"]),recursive=TRUE,showWarnings=FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(data.table); library(readxl); library(GSVA)})

# Private project root replaced by explicit DCF inputs.
# Module paths are explicit DCF fields.
# Resolution output is independently configured.
input_file <- unname(cfg_paths["official_workbook"])
out <- unname(cfg_paths["output_dir"])
source_sha_expected <- "a1f8683676d416b255291dcaa007aa0fa48c736c23e6a7a566664ba528bf205b"
source_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/MediaObjects/41591_2020_839_MOESM2_ESM.xlsx"
log_event <- function(level,msg) { cat(sprintf("[%s] %s script=03_prepare_braun_official.R msg=\"%s\"\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),level,gsub('"',"'",msg))); flush.console() }

write_new <- function(x,path) {
  if(file.exists(path) && Sys.getenv("BRAUN_FORCE","false")!="true") stop("Refusing to overwrite existing output: ",path)
  fwrite(x,path,sep="\t",na="NA",quote=FALSE)
}

main <- function() {
  log_event("START","Official Braun source mapping and outcome-blind nivolumab state assignment")
  if(!file.exists(input_file)) stop("Official source workbook not found: ",input_file)
  if(!dir.exists(out)) dir.create(out,recursive=TRUE,showWarnings=FALSE)
  sha_line <- system2("sha256sum",shQuote(input_file),stdout=TRUE,stderr=TRUE)
  sha_actual <- strsplit(trimws(sha_line[1]),"[[:space:]]+")[[1]][1]
  if(!identical(sha_actual,source_sha_expected)) stop("Official workbook SHA256 mismatch: ",sha_actual)

  lm_path <- unname(cfg_paths["lm12_members"])
  pd_path <- unname(cfg_paths["program_membership"])
  alias_path <- unname(cfg_paths["gene_alias"])
  local_clin_path <- unname(cfg_paths["local_clinical_crosswalk"])
  stopifnot(file.exists(lm_path),file.exists(pd_path),file.exists(alias_path),file.exists(local_clin_path))
  lm_long <- fread(lm_path)
  stopifnot(all(c("feature_id","gene") %in% names(lm_long)))
  lm_sets <- split(toupper(trimws(lm_long$gene)),lm_long$feature_id)
  lm_sets <- lapply(lm_sets,function(z)sort(unique(z[!is.na(z)&nzchar(z)])))
  if(length(lm_sets)!=12L || any(lengths(lm_sets)==0L)) stop("Frozen LM12 membership is not the expected 12 nonempty sets")
  acid <- c("ATP6AP1","GAPDH","HIF1A","HK2","PFKP","PGK1","PKM","SLC2A1","SLC9A7")
  pd <- fread(pd_path)
  required_pd <- c("program_id","gene","include_in_definition_independent_score")
  stopifnot(all(required_pd %in% names(pd)))
  prog_ids <- c(Hypoxia="MALIGNANT_HYPOXIA",ECM="STROMAL_ECM_ORGANIZATION",Cytotoxic="TNK_CYTOTOXICITY_EFFECTOR")
  programs <- lapply(prog_ids,function(id)sort(unique(toupper(pd[program_id==id & include_in_definition_independent_score==TRUE]$gene))))
  names(programs) <- names(prog_ids)
  if(!setequal(names(programs),c("Hypoxia","ECM","Cytotoxic")) || length(programs$Cytotoxic)!=6L) stop("Frozen program membership does not match expected definitions")
  alias <- fread(alias_path)
  stopifnot(all(c("alias_symbol","canonical_symbol") %in% names(alias)))
  alias <- alias[!is.na(alias_symbol)&!is.na(canonical_symbol)&nzchar(alias_symbol)&nzchar(canonical_symbol)]
  alias[,`:=`(alias_symbol=toupper(trimws(alias_symbol)),canonical_symbol=toupper(trimws(canonical_symbol)))]
  if(anyDuplicated(alias$alias_symbol)) stop("Alias manifest has duplicate alias symbols")
  alias_map <- setNames(alias$canonical_symbol,alias$alias_symbol)

  # Read S1 once. Keep outcomes out of the assignment objects until labels are fixed.
  s1_all <- as.data.table(read_excel(input_file,sheet="S1_Clinical_and_Immune_Data",skip=1,.name_repair="minimal"))
  required_s1 <- c("SUBJID","Cohort","Arm","RNA_ID","SampleType","Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy","ORR")
  if(!all(required_s1 %in% names(s1_all))) stop("S1 is missing required clinical fields")
  s1_all[,`:=`(SUBJID=trimws(as.character(SUBJID)),Cohort=trimws(as.character(Cohort)),Arm=toupper(trimws(as.character(Arm))),RNA_ID=trimws(as.character(RNA_ID)),SampleType=trimws(as.character(SampleType)))]
  niv0 <- s1_all[Arm=="NIVOLUMAB" & Cohort %in% c("CM-009","CM-010","CM-025")]
  valid_id <- !is.na(niv0$RNA_ID) & nzchar(niv0$RNA_ID) & !toupper(niv0$RNA_ID)%in%c("NA","N/A","NONE")
  blind <- niv0[valid_id,.(patient_id=SUBJID,trial_source=Cohort,arm=Arm,sample_id=RNA_ID,sample_type=SampleType,days_collection_to_trial_start=Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy)]
  blind[,days_collection_to_trial_start:=suppressWarnings(as.numeric(as.character(days_collection_to_trial_start)))]
  if(nrow(blind)!=181L) stop("Expected 181 RNA-bearing NIVOLUMAB records, observed ",nrow(blind))
  if(anyDuplicated(blind$sample_id) || anyDuplicated(blind$patient_id)) stop("Nivolumab source has duplicate patient or RNA identifiers")
  setorder(blind,sample_id)
  source_all_ids <- s1_all[!is.na(RNA_ID)&nzchar(RNA_ID)&!toupper(RNA_ID)%in%c("NA","N/A","NONE"),unique(RNA_ID)]

  # Direct key join is by RNA_ID; no column-position inference is used.
  s4 <- read_excel(input_file,sheet="S4A_RNA_Expression",skip=1,.name_repair="minimal")
  if(!"gene_name" %in% names(s4)) stop("S4A gene_name field missing")
  s4_ids <- setdiff(names(s4),"gene_name")
  if(length(s4_ids)!=311L || !setequal(source_all_ids,s4_ids)) stop("S1/S4A RNA_ID sets are not an exact 311-ID match")
  if(!all(blind$sample_id %in% s4_ids)) stop("Some eligible NIVOLUMAB RNA_IDs are absent from S4A")
  raw_gene <- trimws(as.character(s4$gene_name))
  raw_gene_u <- toupper(raw_gene)
  frozen_genes <- unique(c(acid,unlist(lm_sets,use.names=FALSE),unlist(programs,use.names=FALSE)))
  raw_dup <- unique(raw_gene_u[!is.na(raw_gene_u)&nzchar(raw_gene_u)&(duplicated(raw_gene_u)|duplicated(raw_gene_u,fromLast=TRUE))])
  is_excel_date_artifact <- function(x) x %in% c("43525","43526") | grepl("^2019-03-0[12]",x)
  date_artifact <- !is.na(raw_gene_u) & vapply(raw_gene_u,is_excel_date_artifact,logical(1))
  duplicate_table <- rbindlist(lapply(raw_dup,function(g) {
    ii <- which(raw_gene_u==g)
    data.table(raw_gene_symbol=g,rows=length(ii),is_frozen_member=g %in% frozen_genes,handling=if(is_excel_date_artifact(g))"excluded_excel_date_serial_for_MARCH1_or_MARCH2; no frozen member" else "unresolved_duplicate")
  }),fill=TRUE)
  if(!nrow(duplicate_table)) duplicate_table <- data.table(raw_gene_symbol=character(),rows=integer(),is_frozen_member=logical(),handling=character())
  if(any(!vapply(raw_dup,is_excel_date_artifact,logical(1)))) stop("Unexpected duplicate S4A gene symbols require source resolution: ",paste(raw_dup[!vapply(raw_dup,is_excel_date_artifact,logical(1))],collapse=","))
  if(any(raw_dup %in% frozen_genes)) stop("An Excel date-converted duplicate overlaps a frozen member")
  use_row <- !is.na(raw_gene_u) & nzchar(raw_gene_u) & !date_artifact
  gene <- raw_gene_u[use_row]
  mapped <- gene %in% names(alias_map)
  gene[mapped] <- unname(alias_map[gene[mapped]])
  if(anyDuplicated(gene)) {
    collision <- unique(gene[duplicated(gene)|duplicated(gene,fromLast=TRUE)])
    stop("Alias mapping creates duplicate logged-scale rows; no averaging performed: ",paste(collision,collapse=","))
  }
  # Only the 181 eligible columns enter state construction; columns are explicitly reordered by sample_id.
  expr_df <- as.data.frame(s4[use_row,blind$sample_id,drop=FALSE],check.names=FALSE)
  expr <- as.matrix(expr_df)
  storage.mode(expr) <- "double"
  rownames(expr) <- gene
  colnames(expr) <- blind$sample_id
  incomplete <- rowSums(!is.finite(expr))>0L
  incomplete_genes <- rownames(expr)[incomplete]
  expr_complete <- expr[!incomplete,,drop=FALSE]
  if(nrow(expr_complete)<10000L) stop("Unexpectedly few complete source genes across eligible nivolumab samples: ",nrow(expr_complete))

  sets <- c(list(Acid_Score=acid),lm_sets)
  cover <- list()
  for(nm in names(sets)) {
    g <- sets[[nm]]
    seen <- intersect(g,rownames(expr_complete))
    variable_n <- if(length(seen)) sum(apply(expr_complete[seen,,drop=FALSE],1,sd)>0) else 0L
    pass <- if(nm=="Acid_Score") length(seen)==9L&&variable_n==9L else length(seen)/length(g)>=.7&&variable_n>0L
    cover[[length(cover)+1L]] <- data.table(cohort="Braun",score=nm,frozen_n=length(g),present_n=length(seen),variable_n=variable_n,coverage_pass=pass,missing=paste(setdiff(g,seen),collapse="|"),rows_removed_for_nonfinite=sum(g %in% incomplete_genes))
  }
  if(!all(vapply(tail(cover,length(sets)),function(z)isTRUE(z$coverage_pass),logical(1)))) stop("Acid/LM12 frozen membership coverage failed; see diagnostics before assignment")

  # Existing frozen assignment: Acid nine-gene z matrix, ward.D/k=3; LM12 ssGSEA, ward.D/k=2.
  gs <- c(list(Acid_Score=acid),lm_sets)
  ss <- GSVA::gsva(GSVA::ssgseaParam(exprData=expr_complete,geneSets=gs,normalize=TRUE),verbose=FALSE)
  acid_score <- as.numeric(ss["Acid_Score",]); names(acid_score)<-colnames(ss)
  lm_score <- ss[names(lm_sets),,drop=FALSE]
  inflammatory_axis <- colMeans(lm_score[setdiff(rownames(lm_score),"LM12_Antimicrobial_Peptides"),,drop=FALSE])
  az <- expr_complete[acid,,drop=FALSE]
  az_sd <- apply(az,1,sd)
  if(any(!is.finite(az_sd)|az_sd==0)) stop("Acid member failed cohort z-score variance check")
  az <- sweep(sweep(az,1,rowMeans(az),"-"),1,az_sd,"/")
  hc_a <- hclust(dist(t(az)),method="ward.D")
  ca <- cutree(hc_a,k=3L); means_a <- tapply(acid_score,ca,mean)
  if(anyDuplicated(as.numeric(means_a))) stop("Acid direction cluster means tied; cannot assign AG labels outcome-blind")
  ag_map <- setNames(paste0("AG",seq_along(means_a)),names(sort(means_a,decreasing=TRUE)))
  hc_i <- hclust(dist(t(lm_score)),method="ward.D")
  ci <- cutree(hc_i,k=2L); means_i <- tapply(inflammatory_axis,ci,mean)
  if(anyDuplicated(as.numeric(means_i))) stop("LM12 axis cluster means tied; cannot assign IG labels outcome-blind")
  ig_map <- setNames(paste0("IG",seq_along(means_i)),names(sort(means_i,decreasing=TRUE)))
  assignment <- data.table(sample_id=colnames(expr_complete),AG=unname(ag_map[as.character(ca)]),IG=unname(ig_map[as.character(ci)]),Acid_Score=acid_score[colnames(expr_complete)],LM12_inflammatory_axis=inflammatory_axis[colnames(expr_complete)])
  assignment[,C_state:=fcase(AG=="AG1"&IG=="IG1","C1",AG=="AG3"&IG=="IG2","C4",AG=="AG1"&IG=="IG2","C2",AG %in% c("AG2","AG3")&IG=="IG1","C3",AG=="AG2"&IG=="IG2","C2")]
  if(anyNA(assignment$C_state) || nrow(assignment)!=181L) stop("State mapping incomplete")
  assignment <- merge(assignment,blind,by="sample_id",all.x=TRUE,sort=FALSE)
  setorder(assignment,sample_id)

  # Fixed RNA programs are per-gene cohort z scores followed by a mean. Coverage failures remain NA.
  program_coverage <- list()
  for(nm in names(programs)) {
    g <- programs[[nm]]
    seen <- intersect(g,rownames(expr_complete))
    v <- expr_complete[seen,,drop=FALSE]
    sdg <- if(length(seen)) apply(v,1,sd) else numeric()
    variable <- names(sdg)[is.finite(sdg)&sdg>0]
    pass <- if(nm=="Cytotoxic") length(seen)==6L&&length(variable)>=5L else length(seen)/length(g)>=.8&&length(variable)>=10L
    missing_all <- setdiff(g,rownames(expr))
    incomplete_members <- intersect(g,incomplete_genes)
    program_coverage[[length(program_coverage)+1L]] <- data.table(cohort="Braun",score=nm,frozen_n=length(g),present_n=length(seen),variable_n=length(variable),coverage_pass=pass,missing=paste(unique(c(missing_all,incomplete_members)),collapse="|"),rows_removed_for_nonfinite=length(incomplete_members))
    value <- rep(NA_real_,ncol(expr_complete)); names(value)<-colnames(expr_complete)
    if(pass) {
      vv <- v[variable,,drop=FALSE]
      z <- sweep(sweep(vv,1,rowMeans(vv),"-"),1,sdg[variable],"/")
      value <- colMeans(z)
      if(is.finite(sd(value))&&sd(value)>0) value <- as.numeric(scale(value)) else value[] <- NA_real_
      names(value) <- colnames(expr_complete)
    }
    assignment[, (nm):=value[colnames(expr_complete)]]
  }
  coverage_all <- rbindlist(c(cover,program_coverage),fill=TRUE)

  # Persist outcome-blind assignments before reading/attaching ORR.
  setorder(assignment,sample_id)
  assignment[,cohort:="Braun"]
  write_new(assignment,file.path(out,"B_Braun_outcome_blind_assignments.tsv"))
  state_counts <- assignment[,.(n=.N),by=.(cohort,C_state)]
  write_new(state_counts,file.path(out,"B_Braun_state_reference_counts.tsv"))
  write_new(coverage_all,file.path(out,"B_Braun_member_coverage.tsv"))
  write_new(duplicate_table,file.path(out,"B_Braun_duplicate_gene_handling.tsv"))

  # Only now attach S1 ORR, and validate the prior local patient mapping.
  s1_outcome <- s1_all[RNA_ID %in% blind$sample_id,.(sample_id=RNA_ID,patient_id=SUBJID,trial_source=Cohort,arm=Arm,raw_response=toupper(trimws(as.character(ORR))),source_sample_type=SampleType,days_collection_to_trial_start=suppressWarnings(as.numeric(as.character(Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy))))]
  if(nrow(s1_outcome)!=181L || anyDuplicated(s1_outcome$sample_id)) stop("S1 outcome mapping is not one-to-one with 181 blind assignments")
  patients <- merge(assignment,s1_outcome,by=c("sample_id","patient_id","trial_source","arm","days_collection_to_trial_start"),all.x=TRUE,sort=FALSE)
  if(nrow(patients)!=181L || anyNA(patients$raw_response)) stop("S1 outcome join changed the assigned patient set")
  patients[,strict_response:=fcase(raw_response %in% c("CR","PR","CRPR"),"R",raw_response %in% c("SD","PD"),"NR",default=NA_character_)]
  patients[,strict_response_included:=!is.na(strict_response)]
  patients[,`:=`(cohort="Braun",regimen=trial_source,baseline_choice_rank=1L)]
  patients[,outcome_mapping_status:=fcase(strict_response=="R","CR_PR_or_CRPR_as_R",strict_response=="NR","SD_or_PD_as_NR",is.na(strict_response),paste0("excluded_strict_",raw_response))]
  local <- fread(local_clin_path)
  req_local <- c("patient_name","sample_id","Treatment","Therapy","response")
  if(!all(req_local %in% names(local))) stop("Local Braun clinical TSV lacks expected columns")
  local[,local_key:=gsub("-","_",trimws(as.character(sample_id)),fixed=TRUE)]
  patients[,local_key:=gsub("-","_",sample_id,fixed=TRUE)]
  if(anyDuplicated(local$local_key)) stop("Local Braun TSV has duplicated normalized sample keys")
  local_map <- local[match(patients$local_key,local$local_key)]
  patients[,`:=`(local_tsv_match=!is.na(local_map$sample_id),local_Treatment=local_map$Treatment,local_Therapy=local_map$Therapy,local_response=toupper(trimws(as.character(local_map$response))),local_patient_name=local_map$patient_name)]
  patients[,local_response_match:=ifelse(local_tsv_match,local_response==raw_response,NA)]
  patients[,local_mapping_status:=fcase(!local_tsv_match,"no_local_tsv_row_expected_for_NE_or_source_extra",local_Treatment!="PRE","local_timepoint_not_PRE",local_response_match==TRUE,"PRE_response_exact_match",local_response_match==FALSE,"RESPONSE_CONFLICT",default="local_mapping_unresolved")]
  if(any(patients$local_mapping_status=="RESPONSE_CONFLICT",na.rm=TRUE)) stop("Official S1 ORR conflicts with local Braun response mapping")
  # The local TSV had pre-dropped NE; retaining it only as a crosswalk must not define reference eligibility.
  local_map_summary <- patients[,.(cohort="Braun",n_source_nivolumab=.N,n_local_tsv_match=sum(local_tsv_match),n_local_response_match=sum(local_response_match,na.rm=TRUE),n_no_local_row=sum(!local_tsv_match),n_strict_response=sum(strict_response_included),n_NE_or_unmapped=sum(!strict_response_included))]
  write_new(local_map_summary,file.path(out,"B_Braun_local_source_mapping_summary.tsv"))
  write_new(patients,file.path(out,"B_Braun_patient_inclusion_and_scores.tsv"))

  timing <- patients[,.(cohort="Braun",n_baseline=.N,n_days_known=sum(is.finite(days_collection_to_trial_start)),n_days_missing=sum(!is.finite(days_collection_to_trial_start)),median_days=if(any(is.finite(days_collection_to_trial_start)))median(days_collection_to_trial_start,na.rm=TRUE) else NA_real_,q1_days=if(any(is.finite(days_collection_to_trial_start)))quantile(days_collection_to_trial_start,.25,na.rm=TRUE,names=FALSE) else NA_real_,q3_days=if(any(is.finite(days_collection_to_trial_start)))quantile(days_collection_to_trial_start,.75,na.rm=TRUE,names=FALSE) else NA_real_,n_gt365=sum(days_collection_to_trial_start>365,na.rm=TRUE),pct_gt365_known=if(sum(is.finite(days_collection_to_trial_start)))mean(days_collection_to_trial_start[is.finite(days_collection_to_trial_start)]>365) else NA_real_,n_negative_days=sum(days_collection_to_trial_start<0,na.rm=TRUE),n_FFPE=sum(source_sample_type=="FFPE",na.rm=TRUE),n_FROZEN=sum(source_sample_type=="FROZEN",na.rm=TRUE)),by=trial_source]
  write_new(timing,file.path(out,"B_Braun_baseline_timing_summary.tsv"))
  scale_info <- data.table(cohort="Braun",source_file=basename(input_file),source_url=source_url,source_sha256=sha_actual,source_scale="author-standardized log matrix: upper-quartile TPM normalized to 1000, log2 transformed, ComBat corrected across four batches, then +21 offset",analysis_scale="use S4A values directly; no additional log transform",n_nivolumab_reference=181L,n_strict_response=sum(patients$strict_response_included),n_NE=sum(!patients$strict_response_included),n_source_genes=nrow(expr),n_complete_genes_used=nrow(expr_complete),n_rows_removed_any_nonfinite=sum(incomplete),GSVA_version=as.character(packageVersion("GSVA")),readxl_version=as.character(packageVersion("readxl")),data_table_version=as.character(packageVersion("data.table")),ssGSEA_normalize=TRUE,cluster_method="Euclidean + ward.D; Acid z matrix k=3; LM12 ssGSEA k=2",sample_order="RNA_ID ascending before clustering",assignment_reference="all RNA-bearing NIVOLUMAB S1 samples including NE; everolimus excluded")
  write_new(scale_info,file.path(out,"B_Braun_expression_scale_and_parameters.tsv"))
  log_event("CHECKPOINT",sprintf("Nivolumab patients=%d; strict=%d; R=%d; NR=%d; NE=%d; state counts: %s",nrow(patients),sum(patients$strict_response_included),sum(patients$strict_response=="R",na.rm=TRUE),sum(patients$strict_response=="NR",na.rm=TRUE),sum(!patients$strict_response_included),paste(paste(state_counts$C_state,state_counts$n,sep="="),collapse=",")))
  log_event("DONE","Braun source crosswalk, outcome-blind assignments, membership coverage, and patient manifest written")
}

tryCatch(main(),error=function(e){log_event("ERROR",conditionMessage(e)); stop(e)})
