#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 07_analyze_liu_exploratory.R --config PATH [--check-config]\n       Rscript 07_analyze_liu_exploratory.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("clinical_workbook","outcome_blind_assignments","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("clinical_workbook","outcome_blind_assignments")]
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
# 普通运行才创建独立空输出；检查模式已退出。
dir.create(unname(cfg_paths["output_dir"]),recursive=TRUE,showWarnings=FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table); library(readxl)})
# Private project root replaced by explicit DCF inputs.
# Module paths are supplied by the checked DCF.
source_clin <- unname(cfg_paths["clinical_workbook"])
res_dir <- unname(cfg_paths["output_dir"])
log_event <- function(level,msg) {
  cat(sprintf("[%s] %s script=07_analyze_liu_exploratory.R msg=\"%s\"\n",
              format(Sys.time(),"%Y-%m-%d %H:%M:%S"),level,msg)); flush.console()
}
log_event("START","link best RECIST response only after outcome-blind assignments; estimate fixed comparisons")
blind_all <- fread(unname(cfg_paths["outcome_blind_assignments"]))
primary <- blind_all[sample_set=="primary_Pre_PD1_context3_n103"]
broad <- blind_all[sample_set=="sensitivity_all_pre_PD1_context1to3_n119"]
stopifnot(nrow(primary)==103L,nrow(broad)==119L)

# 仅在两套状态赋类完成后读取 Table S1 的 BR 列并连接患者结局。
read_col <- function(rng) readxl::read_xlsx(source_clin,sheet="Supplemental Table 1",range=rng,col_names=FALSE)
patient_id <- as.character(read_col("A4:A161")[[1]])
br <- as.character(read_col("L4:L161")[[1]])
clin <- data.table(patient_id=patient_id,raw_BR=br)
clin <- clin[!is.na(patient_id)&nzchar(patient_id)]
stopifnot(!anyDuplicated(clin$patient_id))
join_br <- function(x) {
  z <- merge(x,clin,by="patient_id",all.x=TRUE,sort=FALSE)
  if(nrow(z)!=nrow(x)) stop("Response mapping changed row count")
  if(anyNA(z$raw_BR)) stop("At least one assigned patient lacks a best-response mapping")
  if(any(!z$raw_BR %in% c("CR","PR","SD","PD","MR"))) stop("Unexpected BR class in Liu Table S1")
  z[,strict_response:=fcase(raw_BR %in% c("CR","PR"),"R",raw_BR %in% c("SD","PD"),"NR",default=NA_character_)]
  z[,CRPR_vs_PD:=fcase(raw_BR %in% c("CR","PR"),"R",raw_BR=="PD","NR",default=NA_character_)]
  z
}
primary <- join_br(primary)
broad <- join_br(broad)
fwrite(primary,file.path(res_dir,"Liu_primary103_patient_response_mapping.tsv"),sep="\t",na="NA")
fwrite(broad,file.path(res_dir,"Liu_broad119_patient_response_mapping.tsv"),sep="\t",na="NA")

states <- c("C1","C2","C3","C4")
response_rates <- function(d,set_name) {
  rbindlist(lapply(states,function(st) {
    z <- d[C_state==st]
    e <- z[!is.na(strict_response)]
    nr <- sum(e$strict_response=="R"); nnr <- sum(e$strict_response=="NR")
    ci <- if(nr+nnr>0L) binom.test(nr,nr+nnr)$conf.int else c(NA_real_,NA_real_)
    data.table(analysis_set=set_name,C_state=st,n_baseline=nrow(z),n_CRPR_SD_PD=nr+nnr,
               CRPR=nr,SD_PD=nnr,MR=sum(z$raw_BR=="MR"),unmapped_response=sum(is.na(z$raw_BR)),
               response_rate=if(nr+nnr>0L)nr/(nr+nnr) else NA_real_,ci_lo=ci[1],ci_hi=ci[2])
  }))
}
rate_main <- response_rates(primary,"Primary Pre-PD1 context 3 (n=103)")
rate_broad <- response_rates(broad,"Timing sensitivity contexts 1-3 (n=119)")
fwrite(rate_main,file.path(res_dir,"Liu_primary103_four_state_response_rates.tsv"),sep="\t",na="NA")
fwrite(rate_broad,file.path(res_dir,"Liu_broad119_four_state_response_rates.tsv"),sep="\t",na="NA")

fisher_row <- function(d,ref,response_col="strict_response") {
  y <- d[get(response_col) %in% c("R","NR") & C_state %in% c("C1",ref)]
  a <- y[C_state=="C1"]; b <- y[C_state==ref]
  base <- list(n_C1=nrow(a),n_ref=nrow(b),R_C1=sum(a[[response_col]]=="R"),NR_C1=sum(a[[response_col]]=="NR"),
               R_ref=sum(b[[response_col]]=="R"),NR_ref=sum(b[[response_col]]=="NR"),
               status="estimated",OR=NA_real_,OR_lo=NA_real_,OR_hi=NA_real_,p=NA_real_)
  if(min(nrow(a),nrow(b))<10L) {base$status <- "under_10_endpoint_patients_per_state"; return(base)}
  if(sum(y[[response_col]]=="R")==0L || sum(y[[response_col]]=="NR")==0L) {base$status <- "single_endpoint_in_pooled_pair"; return(base)}
  tab <- matrix(c(base$R_C1,base$NR_C1,base$R_ref,base$NR_ref),nrow=2,byrow=TRUE)
  ft <- fisher.test(tab,alternative="two.sided")
  base$OR <- unname(ft$estimate); base$OR_lo <- ft$conf.int[1]; base$OR_hi <- ft$conf.int[2]; base$p <- ft$p.value
  base
}
contrasts <- c("C4","C3")
main_cmp <- rbindlist(lapply(contrasts,function(ref) {
  z <- fisher_row(primary,ref)
  as.data.table(c(list(analysis_set="primary_Pre_PD1_context3_n103",contrast=paste0("C1/",ref)),z))
}))
main_cmp[,q:=NA_real_]
ii <- which(is.finite(main_cmp$p)); if(length(ii)) main_cmp$q[ii] <- p.adjust(main_cmp$p[ii],method="BH",n=2L)
fwrite(main_cmp,file.path(res_dir,"Liu_primary103_exact_comparisons.tsv"),sep="\t",na="NA")

# 119人按 context 1-3 重新独立赋类，只作方向/效应敏感性，不作第二套显著性胜负。
broad_cmp <- rbindlist(lapply(contrasts,function(ref) {
  z <- fisher_row(broad,ref)
  as.data.table(c(list(analysis_set="timing_sensitivity_contexts1to3_n119",contrast=paste0("C1/",ref)),z))
}))
broad_cmp[,q:=NA_real_]
fwrite(broad_cmp,file.path(res_dir,"Liu_broad119_timewindow_sensitivity.tsv"),sep="\t",na="NA")

# 固定终点敏感性：CR/PR 对 PD；SD、MR 不纳入本次终点，仍保留在主状态参考中。
pd_sens <- rbindlist(lapply(contrasts,function(ref) {
  z <- fisher_row(primary,ref,response_col="CRPR_vs_PD")
  as.data.table(c(list(analysis_set="primary_Pre_PD1_context3_n103",contrast=paste0("C1/",ref)),z))
}))
fwrite(pd_sens,file.path(res_dir,"Liu_primary103_CRPR_vs_PD_sensitivity.tsv"),sep="\t",na="NA")

# Jeffreys 惩罚似然 Firth logistic，置信区间用 profile penalized likelihood。
firth_fit <- function(X,y,index) {
  X <- as.matrix(X); storage.mode(X) <- "double"; y <- as.numeric(y)
  if(qr(X)$rank<ncol(X)) stop("Firth design is rank deficient")
  ll <- function(b) {
    eta <- drop(X%*%b); pr <- plogis(eta); w <- pmax(pr*(1-pr),1e-15)
    inf <- crossprod(X*sqrt(w)); ld <- determinant(inf,logarithm=TRUE)
    if(ld$sign<=0 || !is.finite(ld$modulus)) return(-1e50)
    sum(y*eta-(pmax(eta,0)+log1p(exp(-abs(eta)))))+.5*as.numeric(ld$modulus)
  }
  full <- optim(rep(0,ncol(X)),function(b)-ll(b),method="BFGS",control=list(maxit=1000,reltol=1e-10))
  if(full$convergence!=0) stop("Firth full optimization failed")
  bh <- full$par; top <- ll(bh)
  profile_ll <- function(v) {
    rem <- setdiff(seq_along(bh),index)
    fun <- function(z) { b<-bh; b[index]<-v; b[rem]<-z; -ll(b) }
    op <- optim(bh[rem],fun,method="BFGS",control=list(maxit=1000,reltol=1e-9))
    if(op$convergence!=0) stop("Firth profile optimization failed")
    -op$value
  }
  target <- top-qchisq(.95,1)/2
  bound <- function(side) {
    step <- .5; b <- bh[index]+side*step
    while(profile_ll(b)>target && step<32) {step<-step*2;b<-bh[index]+side*step}
    if(profile_ll(b)>target) return(if(side<0)-Inf else Inf)
    if(side<0) uniroot(function(z)profile_ll(z)-target,c(b,bh[index]),tol=1e-5)$root
    else uniroot(function(z)profile_ll(z)-target,c(bh[index],b),tol=1e-5)$root
  }
  p <- pchisq(2*(top-profile_ll(0)),df=1,lower.tail=FALSE)
  c(OR=exp(bh[index]),lo=exp(bound(-1)),hi=exp(bound(1)),p=p)
}

clinical <- rbindlist(lapply(contrasts,function(ref) {
  x <- primary[C_state %in% c("C1",ref) & !is.na(strict_response)]
  covariate <- "priorCTLA4"
  z <- x[!is.na(get(covariate))]
  complete <- if(nrow(x)>0L)nrow(z)/nrow(x) else NA_real_
  nlev <- uniqueN(z[[covariate]])
  npar <- 1L+as.integer(nlev>1L)
  minor <- min(sum(z$strict_response=="R"),sum(z$strict_response=="NR"))
  status <- if(is.na(complete)||complete<.8) "covariate_completeness_under_80pct" else if(minor<5L*npar) "minor_outcome_under_5_per_parameter" else if(min(nrow(z[C_state=="C1"]),nrow(z[C_state==ref]))<10L) "under_10_per_state" else "estimated"
  row <- data.table(contrast=paste0("C1/",ref),covariate=covariate,n_main=nrow(x),n_complete=nrow(z),covariate_completeness=complete,
                    n_priorCTLA4_0=sum(z$priorCTLA4==0),n_priorCTLA4_1=sum(z$priorCTLA4==1),status=status,
                    unadjusted_same_case_OR=NA_real_,unadjusted_lo=NA_real_,unadjusted_hi=NA_real_,firth_OR=NA_real_,firth_lo=NA_real_,firth_hi=NA_real_,firth_p=NA_real_)
  if(status=="estimated") {
    u <- fisher_row(z,ref)
    row[,`:=`(unadjusted_same_case_OR=u$OR,unadjusted_lo=u$OR_lo,unadjusted_hi=u$OR_hi)]
    Z <- data.frame(state=as.integer(z$C_state=="C1"),cov=factor(z[[covariate]]))
    X <- if(nlev>1L) model.matrix(~state+cov,Z) else model.matrix(~state,Z)
    fit <- tryCatch(firth_fit(X,as.integer(z$strict_response=="R"),which(colnames(X)=="state")),error=function(e)e)
    if(inherits(fit,"error")) row[,status:=paste0("firth_fit_failed: ",conditionMessage(fit))] else row[,`:=`(firth_OR=fit["OR"],firth_lo=fit["lo"],firth_hi=fit["hi"],firth_p=fit["p"])]
  }
  row
}))
fwrite(clinical,file.path(res_dir,"Liu_primary103_priorCTLA4_Firth_sensitivity.tsv"),sep="\t",na="NA")

# 两个连续轴是固定辅助分析；按患者级总体响应建模，不寻优切点。
continuous <- rbindlist(lapply(c("Acid_Score","LM12_inflammatory_axis"),function(axis) {
  z_all <- as.numeric(scale(primary[[axis]]))
  d <- copy(primary); d[,z:=z_all]
  x <- d[!is.na(strict_response)&is.finite(z)]
  nR <- sum(x$strict_response=="R"); nNR <- sum(x$strict_response=="NR")
  row <- data.table(axis=axis,n_R=nR,n_NR=nNR,status="under_10_per_outcome",OR_per_SD=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,q=NA_real_)
  if(min(nR,nNR)>=10L) {
    fit <- tryCatch(firth_fit(model.matrix(~z,x),as.integer(x$strict_response=="R"),2L),error=function(e)e)
    if(inherits(fit,"error")) row[,status:=paste0("firth_fit_failed: ",conditionMessage(fit))] else row[,`:=`(status="estimated",OR_per_SD=fit["OR"],lo=fit["lo"],hi=fit["hi"],p=fit["p"])]
  }
  row
}))
ii <- which(is.finite(continuous$p)); if(length(ii)) continuous$q[ii] <- p.adjust(continuous$p[ii],method="BH",n=2L)
fwrite(continuous,file.path(res_dir,"Liu_primary103_continuous_Firth_auxiliary.tsv"),sep="\t",na="NA")

# C1 内部三程序固定采用 Wilcoxon + 秩二列效应量，2000次分层 bootstrap。
set.seed(20260923)
rank_biserial <- function(a,b) {
  U <- unname(wilcox.test(a,b,exact=FALSE,correct=FALSE)$statistic)
  2*U/(length(a)*length(b))-1
}
internal <- rbindlist(lapply(c("Hypoxia","ECM","Cytotoxic"),function(program) {
  x <- primary[C_state=="C1" & !is.na(strict_response) & is.finite(get(program))]
  a <- x[strict_response=="R"][[program]]; b <- x[strict_response=="NR"][[program]]
  row <- data.table(program=program,n_R=length(a),n_NR=length(b),status="under_10_per_response_group_or_missing_score",
                    rank_biserial=NA_real_,ci_lo=NA_real_,ci_hi=NA_real_,p=NA_real_,q=NA_real_)
  if(min(length(a),length(b))>=10L) {
    boots <- replicate(2000L,rank_biserial(sample(a,length(a),replace=TRUE),sample(b,length(b),replace=TRUE)))
    ci <- quantile(boots,c(.025,.975),names=FALSE)
    row[,`:=`(status="estimated",rank_biserial=rank_biserial(a,b),ci_lo=ci[1],ci_hi=ci[2],p=wilcox.test(a,b,alternative="two.sided",exact=FALSE)$p.value)]
  }
  row
}))
ii <- which(is.finite(internal$p)); if(length(ii)) internal$q[ii] <- p.adjust(internal$p[ii],method="BH",n=3L)
fwrite(internal,file.path(res_dir,"Liu_primary103_C1_internal_programme_exploration.tsv"),sep="\t",na="NA")

log_event("DONE",sprintf("Liu primary n=103; exact slots %d/2; continuous %d/2; C1 programmes %d/3; broad timing sensitivity n=119",sum(main_cmp$status=="estimated"),sum(continuous$status=="estimated"),sum(internal$status=="estimated")))
