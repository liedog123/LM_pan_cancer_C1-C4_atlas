#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 14_prepare_subtype_context.R --config PATH [--check-config]\n       Rscript 14_prepare_subtype_context.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_manifest", "canonical_labels", "immune_subtypes", "molecular_subtypes", "program_membership", "mrna_rdata", "tcpa_table13", "antibody_dictionary", "rppa_matrix", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
protected <- c(input_paths[dir.exists(input_paths)],dirname(required_files),config_path)
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

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
main <- function() {
  # 普通执行才创建独立输出；源计算连续保留，成功后原子提交。
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  final_dir <- file.path(out_dir, "01_subtype_context_preparation")
  if(file.exists(final_dir) || dir.exists(final_dir)) stop("Output stage already exists.")
  out <- paste0(final_dir, ".tmp_", Sys.getpid())
  if(file.exists(out) || dir.exists(out)) stop("Temporary stage already exists.")
  if(!dir.create(out, showWarnings=FALSE)) stop("Cannot create temporary stage.")
log_event <- function(level, msg) { cat(sprintf("[%s] %s script=14_prepare_subtype_context.R msg=\"%s\"\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)); flush.console() }
log_event("START", "TCGA patient, subtype, RNA and RPPA preparation; frozen v1.0 definitions")
man <- fread(paths["patient_manifest"])
lab <- fread(paths["canonical_labels"], select=c("Sample","Cancer","C_subtype"))
stopifnot(!anyDuplicated(man$Patient_ID), !anyDuplicated(lab[, .(Sample,Cancer)]))
check <- merge(man[, .(Patient_ID,Sample,Cancer,C_subtype)], lab, by=c("Sample","Cancer"), suffixes=c("_manifest","_canonical"), all.x=TRUE)
stopifnot(nrow(check)==nrow(man), all(check$C_subtype_manifest==check$C_subtype_canonical))
imm <- fread(paths["immune_subtypes"])
mol <- fread(paths["molecular_subtypes"], select=c("sampleID","Subtype_mRNA"))
setnames(mol,"sampleID","sample")
resolve <- function(d, values, samples, patients) {
  # 优先代表肿瘤样本；仅在患者的肿瘤条目唯一一致时退到患者级。
  d <- d[grepl("^TCGA-[^-]+-[^-]+-0[1-9]", sample)]
  d[, patient := substr(sample,1,12)]
  ix <- match(samples,d$sample)
  ans <- as.character(d[[values]][ix]); how <- ifelse(!is.na(ix),"exact_sample",NA_character_)
  for (k in which(is.na(ans))) {
    z <- unique(na.omit(as.character(d[patient==patients[k]][[values]])))
    if (length(z)==1L) { ans[k] <- z; how[k] <- "unique_patient_tumour" }
    else if(length(z)>1L) how[k] <- "conflicting_patient_tumour"
    else how[k] <- "missing"
  }
  list(value=ans, route=how)
}
i <- resolve(imm,"Subtype_Immune_Model_Based",man$Sample,man$Patient_ID)
m <- resolve(mol,"Subtype_mRNA",man$Sample,man$Patient_ID)
patients <- copy(man[, .(Patient_ID,Sample,Cancer,C_subtype,AG,IG,TumorPurity,discordant_multi_sample_state)])
patients[, `:=`(immune_subtype=i$value, immune_mapping=i$route, molecular_subtype=ifelse(Cancer=="BRCA",m$value,NA_character_), molecular_mapping=ifelse(Cancer=="BRCA",m$route,NA_character_))]
patients[, immune_subtype := sub(".*\\(Immune (C[1-6])\\).*", "Immune \\1", immune_subtype)]
fwrite(patients,file.path(out,"A_patient_inclusion_and_subtype_mapping.tsv"),sep="\t",na="NA")
log_event("CHECKPOINT",sprintf("mapped %d patients across %d cancers",nrow(patients),uniqueN(patients$Cancer)))
pd <- fread(paths["program_membership"])
ids <- c(Hypoxia="MALIGNANT_HYPOXIA",ECM="STROMAL_ECM_ORGANIZATION",Cytotoxic="TNK_CYTOTOXICITY_EFFECTOR")
sets <- lapply(ids,function(id) sort(unique(pd[program_id==id & include_in_definition_independent_score==TRUE]$gene)))
stopifnot(length(sets$Cytotoxic)==6L)
fwrite(data.table(program=rep(names(sets),lengths(sets)),gene=unlist(sets,use.names=FALSE)),file.path(out,"A_frozen_program_members.tsv"),sep="\t")
e <- new.env(parent=emptyenv())
load(paths["mrna_rdata"],envir=e)
stopifnot(exists("pancancer_data_list",envir=e))
scores <- list(); cover <- list()
for (cc in sort(unique(patients$Cancer))) {
  pp <- patients[Cancer==cc]
  raw <- as.matrix(e$pancancer_data_list[[cc]]$expression_data)
  storage.mode(raw) <- "double"
  if(any(raw<0,na.rm=TRUE)) stop("negative TPM: ",cc)
  si <- match(pp$Sample,colnames(raw)); stopifnot(all(!is.na(si)))
  raw <- raw[,si,drop=FALSE]
  for (program in names(sets)) {
    genes <- sets[[program]]; present <- intersect(genes,rownames(raw))
    v <- log2(raw[present,,drop=FALSE]+1)
    s <- apply(v,1,sd,na.rm=TRUE); valid <- names(s)[is.finite(s)&s>0]
    enough <- if(program=="Cytotoxic") length(present)==6L && length(valid)>=5L else length(present)/length(genes)>=.8 && length(valid)>=10L
    cover[[length(cover)+1L]] <- data.table(Cancer=cc,program=program,frozen_n=length(genes),present_n=length(present),variable_n=length(valid),missing=paste(setdiff(genes,present),collapse="|"),zero_variance=paste(setdiff(present,valid),collapse="|"),coverage_pass=enough)
    z <- rep(NA_real_,ncol(raw))
    if(enough) {
      zz <- sweep(sweep(v[valid,,drop=FALSE],1,rowMeans(v[valid,,drop=FALSE],na.rm=TRUE),"-"),1,s[valid],"/")
      z <- colMeans(zz); z[colSums(!is.finite(zz))>0] <- NA_real_
      z <- as.numeric(scale(z))
    }
    scores[[length(scores)+1L]] <- data.table(Patient_ID=pp$Patient_ID,Cancer=cc,program=program,score=z)
  }
}
rm(e); invisible(gc())
log_event("CHECKPOINT","RNA programmes scored using log2(TPM+1) and within-cancer gene z scores")
# RPPA 五抗体按官方注释和冻结字典定向提取。
pathway <- fread(paths["tcpa_table13"])[Pathway=="EMT"]
dict <- fread(paths["antibody_dictionary"])
dd <- merge(pathway,dict,by.x="Source_Protein_Token",by.y="standardized_feature")
stopifnot(nrow(dd)==5L,all(dd$mapping_status=="OFFICIAL_SINGLE_TARGET"),!anyDuplicated(dd$raw_feature))
rppa <- fread(paths["rppa_matrix"])
rppa <- rppa[SampleID %in% dd$raw_feature]
stopifnot(nrow(rppa)==5L)
for(cc in sort(unique(patients$Cancer))) {
  pp <- patients[Cancer==cc]
  ri <- match(pp$Sample,names(rppa)); available <- !is.na(ri)
  z <- rep(NA_real_,nrow(pp))
  if(sum(available)>=2L) {
    cols <- names(rppa)[ri[available]]
    mat <- as.matrix(rppa[, ..cols])
    storage.mode(mat) <- "double"
    mat <- mat[match(dd$raw_feature,rppa$SampleID),,drop=FALSE]
    rs <- apply(mat,1,sd,na.rm=TRUE)
    if(all(is.finite(rs)&rs>0)) {
      mm <- sweep(sweep(mat,1,rowMeans(mat,na.rm=TRUE),"-"),1,rs,"/")
      avg <- as.numeric(crossprod(dd$Direction_Sign,mm))/5
      avg[colSums(!is.finite(mm))>0] <- NA_real_
      z[available] <- as.numeric(scale(avg))
    }
  }
  cover[[length(cover)+1L]] <- data.table(Cancer=cc,program="RPPA_EMT",frozen_n=5L,present_n=sum(!is.na(match(dd$raw_feature,rppa$SampleID))),variable_n=if(exists("rs")) sum(is.finite(rs)&rs>0) else NA_integer_,missing="",zero_variance="",coverage_pass=any(is.finite(z)))
  scores[[length(scores)+1L]] <- data.table(Patient_ID=pp$Patient_ID,Cancer=cc,program="RPPA_EMT",score=z)
  if(exists("rs")) rm(rs)
}
fwrite(rbindlist(scores),file.path(out,"A_patient_programme_scores.tsv"),sep="\t",na="NA")
fwrite(rbindlist(cover,fill=TRUE),file.path(out,"A_member_coverage.tsv"),sep="\t",na="NA")
log_event("DONE",sprintf("wrote patient inclusion and four programme scores; %d patient-programme rows",sum(vapply(scores,nrow,integer(1)))))

  if(file.exists(final_dir) || dir.exists(final_dir) || !file.rename(out,final_dir)) stop("Atomic output commit failed; temporary files retained: ",out)
  invisible(final_dir)
}
main()
