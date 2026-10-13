#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 01_prepare_riaz_gide.R --config PATH [--check-config]\n       Rscript 01_prepare_riaz_gide.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("lm12_members", "program_membership", "gene_alias", "prior_ipilimumab", "riaz_expression", "riaz_metadata", "gide_expression", "gide_metadata", "braun_hold_metadata", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("lm12_members", "program_membership", "gene_alias", "prior_ipilimumab", "riaz_expression", "riaz_metadata", "gide_expression", "gide_metadata", "braun_hold_metadata")]
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
# 普通运行才创建独立输出；检查模式已退出。
dir.create(unname(cfg_paths["output_dir"]), recursive = TRUE, showWarnings = FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table); library(GSVA)})
# Project paths are supplied by the checked DCF.
out <- unname(cfg_paths["output_dir"])
# Explicit cohort files replace the private input directory.
log_event <- function(level,msg) { cat(sprintf("[%s] %s script=01_prepare_riaz_gide.R msg=\"%s\"\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),level,msg)); flush.console() }
log_event("START","pretreatment expression preparation and outcome-blind within-cohort C1-C4 assignment")
lm_long <- fread(unname(cfg_paths["lm12_members"]))
stopifnot(all(c("feature_id","gene") %in% names(lm_long)))
lm_sets <- split(toupper(lm_long$gene),lm_long$feature_id)
lm_sets <- lapply(lm_sets,unique)
stopifnot(length(lm_sets)==12L)
acid <- c("ATP6AP1","GAPDH","HIF1A","HK2","PFKP","PGK1","PKM","SLC2A1","SLC9A7")
pd <- fread(unname(cfg_paths["program_membership"]))
prog_ids <- c(Hypoxia="MALIGNANT_HYPOXIA",ECM="STROMAL_ECM_ORGANIZATION",Cytotoxic="TNK_CYTOTOXICITY_EFFECTOR")
programs <- lapply(prog_ids,function(id) sort(unique(pd[program_id==id & include_in_definition_independent_score==TRUE]$gene)))
alias <- fread(unname(cfg_paths["gene_alias"]))
alias_map <- setNames(alias$canonical_symbol,alias$alias_symbol)
ipi <- fread(unname(cfg_paths["prior_ipilimumab"]))
cover <- list(); cases <- list(); scale_info <- list(); duplicate_info <- list()
spec <- data.table(cohort=c("Riaz","Gide"),prefix=c("Melanoma-GSE91061","Melanoma-PRJEB23709"))
for(k in seq_len(nrow(spec))) {
  cohort <- spec$cohort[k]; prefix <- spec$prefix[k]
  meta <- fread(unname(cfg_paths[paste0(tolower(cohort),"_metadata")]))
  pre <- meta[Treatment=="PRE" & !is.na(patient_name) & !is.na(sample_id)]
  setorder(pre,patient_name,sample_id)
  pre[, baseline_choice_rank:=seq_len(.N),by=patient_name]
  pre <- pre[baseline_choice_rank==1L]
  # 响应列直到结局盲赋类完成才连接。
  blind <- pre[,.(sample_id,patient_id=patient_name,regimen=Therapy,baseline_choice_rank)]
  setorder(blind,sample_id)
  raw <- as.data.frame(readRDS(unname(cfg_paths[paste0(tolower(cohort),"_expression")])),check.names=FALSE)
  stopifnot(all(blind$sample_id %in% names(raw)))
  mat <- as.matrix(raw[,blind$sample_id,drop=FALSE]); storage.mode(mat) <- "double"
  gene <- toupper(trimws(as.character(raw$GENE_SYMBOL)))
  gene[gene %in% names(alias_map)] <- alias_map[gene[gene %in% names(alias_map)]]
  keep <- !is.na(gene) & nzchar(gene)
  mat <- mat[keep,,drop=FALSE]; rownames(mat) <- gene[keep]
  if(any(!is.finite(mat)) || any(mat<0)) stop(cohort,": expected finite nonnegative FPKM")
  dup <- sum(duplicated(rownames(mat)))
  duplicate_info[[k]] <- data.table(cohort=cohort,duplicate_gene_rows=dup,aggregation="equal_weight_mean_on_linear_scale")
  if(dup) {
    n <- rowsum(matrix(1,nrow(mat),ncol(mat)),group=rownames(mat),reorder=FALSE)
    mat <- rowsum(mat,group=rownames(mat),reorder=FALSE)/n
  }
  fpkm_sum <- colSums(mat)
  if(any(fpkm_sum<=0)) stop(cohort,": zero FPKM column sum")
  # tigeR 来源矩阵为 FPKM；在明确线性尺度上逐样本换算 TPM。
  tpm <- sweep(mat,2,fpkm_sum,"/")*1e6
  expr <- log2(tpm+1)
  scale_info[[k]] <- data.table(cohort=cohort,source_scale="tigeR FPKM",analysis_scale="TPM=FPKM/sum(FPKM)*1e6; log2(TPM+1)",n_patients=ncol(expr),fpkm_sum_min=min(fpkm_sum),fpkm_sum_max=max(fpkm_sum),min_raw=min(mat),max_raw=max(mat),GSVA_version=as.character(packageVersion("GSVA")),ssGSEA_normalize=TRUE)
  sets <- c(list(Acid_Score=acid),lm_sets)
  present <- lapply(sets,intersect,y=rownames(expr))
  for(nm in names(sets)) {
    v <- expr[present[[nm]],,drop=FALSE]
    variable <- sum(apply(v,1,sd)>0)
    enough <- if(nm=="Acid_Score") length(present[[nm]])==9L && variable==9L else length(present[[nm]])/length(sets[[nm]])>=.7 && variable>0
    cover[[length(cover)+1L]] <- data.table(cohort=cohort,score=nm,frozen_n=length(sets[[nm]]),present_n=length(present[[nm]]),variable_n=variable,coverage_pass=enough,missing=paste(setdiff(sets[[nm]],present[[nm]]),collapse="|"))
  }
  if(!all(vapply(tail(cover,length(sets)),function(z)z$coverage_pass,logical(1)))) stop(cohort,": state gene-set coverage failed")
  ss <- GSVA::gsva(GSVA::ssgseaParam(exprData=expr,geneSets=present,normalize=TRUE),verbose=FALSE)
  lm <- ss[names(lm_sets),,drop=FALSE]
  axis <- colMeans(lm[setdiff(rownames(lm),"LM12_Antimicrobial_Peptides"),,drop=FALSE])
  acid_score <- as.numeric(ss["Acid_Score",]); names(acid_score) <- colnames(ss)
  az <- expr[acid,,drop=FALSE]
  az <- sweep(sweep(az,1,rowMeans(az),"-"),1,apply(az,1,sd),"/")
  hc_a <- hclust(dist(t(az)),method="ward.D")
  ca <- cutree(hc_a,k=3); means_a <- tapply(acid_score,ca,mean)
  if(anyDuplicated(as.numeric(means_a))) stop(cohort,": tied AG direction")
  ag_map <- setNames(paste0("AG",1:3),names(sort(means_a,decreasing=TRUE)))
  hc_i <- hclust(dist(t(lm)),method="ward.D")
  ci <- cutree(hc_i,k=2); means_i <- tapply(axis,ci,mean)
  if(anyDuplicated(as.numeric(means_i))) stop(cohort,": tied IG direction")
  ig_map <- setNames(paste0("IG",1:2),names(sort(means_i,decreasing=TRUE)))
  assign <- data.table(sample_id=colnames(expr),AG=unname(ag_map[as.character(ca)]),IG=unname(ig_map[as.character(ci)]),Acid_Score=acid_score[colnames(expr)],LM12_inflammatory_axis=axis[colnames(expr)])
  assign[, C_state:=fcase(AG=="AG1"&IG=="IG1","C1",AG=="AG3"&IG=="IG2","C4",AG=="AG1"&IG=="IG2","C2",AG %in% c("AG2","AG3")&IG=="IG1","C3",AG=="AG2"&IG=="IG2","C2")]
  if(anyNA(assign$C_state)) stop(cohort,": state map incomplete")
  for(nm in names(programs)) {
    gs <- programs[[nm]]; seen <- intersect(gs,rownames(expr)); v <- expr[seen,,drop=FALSE]
    sdg <- apply(v,1,sd); variable <- names(sdg)[is.finite(sdg)&sdg>0]
    pass <- if(nm=="Cytotoxic") length(seen)==6L&&length(variable)>=5L else length(seen)/length(gs)>=.8&&length(variable)>=10L
    cover[[length(cover)+1L]] <- data.table(cohort=cohort,score=nm,frozen_n=length(gs),present_n=length(seen),variable_n=length(variable),coverage_pass=pass,missing=paste(setdiff(gs,seen),collapse="|"))
    val <- rep(NA_real_,ncol(expr))
    if(pass) {
      v <- v[variable,,drop=FALSE]
      val <- colMeans(sweep(sweep(v,1,rowMeans(v),"-"),1,sdg[variable],"/"))
      val <- as.numeric(scale(val))
    }
    assign[, (nm):=val]
  }
  assign <- merge(assign,blind,by="sample_id",sort=FALSE)
  setorder(assign,sample_id)
  # 最后才连接严格 CR/PR vs SD/PD 结局。
  clin <- pre[,.(sample_id,raw_response=response)]
  assign <- merge(assign,clin,by="sample_id",sort=FALSE)
  assign[, strict_response:=fcase(raw_response %in% c("CR","PR"),"R",raw_response %in% c("SD","PD"),"NR",default=NA_character_)]
  assign[, prior_ipilimumab:=if(cohort=="Riaz") ipi$prior_anti_ctla4_status[match(patient_id,ipi$patient_id)] else NA_character_]
  assign[, cohort:=cohort]
  cases[[k]] <- assign
  log_event("CHECKPOINT",sprintf("%s: %d pretreatment patients assigned; C1=%d C2=%d C3=%d C4=%d",cohort,nrow(assign),sum(assign$C_state=="C1"),sum(assign$C_state=="C2"),sum(assign$C_state=="C3"),sum(assign$C_state=="C4")))
  rm(raw,mat,tpm,expr); invisible(gc())
}
fwrite(rbindlist(cases,fill=TRUE),file.path(out,"B_baseline_patient_inclusion_and_scores.tsv"),sep="\t",na="NA")
fwrite(rbindlist(cover,fill=TRUE),file.path(out,"B_member_coverage.tsv"),sep="\t",na="NA")
fwrite(rbindlist(scale_info),file.path(out,"B_expression_scale_and_parameters.tsv"),sep="\t")
fwrite(rbindlist(duplicate_info),file.path(out,"B_duplicate_gene_handling.tsv"),sep="\t")
br <- fread(unname(cfg_paths["braun_hold_metadata"]))[Treatment=="PRE"]
br[, `:=`(cohort="Braun",source_status="HOLD_UNRESOLVED_EXPRESSION_SCALE_AND_TRIAL_SOURCE")]
fwrite(br[,.(cohort,patient_id=patient_name,sample_id,Treatment,Therapy,response,source_status)],file.path(out,"B_Braun_source_hold_patient_list.tsv"),sep="\t",na="NA")
log_event("DONE","Riaz/Gide baseline states and historical Braun HOLD inclusion list written")
