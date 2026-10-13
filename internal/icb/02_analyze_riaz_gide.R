#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 02_analyze_riaz_gide.R --config PATH [--check-config]\n       Rscript 02_analyze_riaz_gide.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_scores", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient_scores")]
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
suppressPackageStartupMessages(library(data.table))
# Project paths are supplied by the checked DCF.
out <- unname(cfg_paths["output_dir"])
log_event <- function(level,msg) { cat(sprintf("[%s] %s script=02_analyze_riaz_gide.R msg=\"%s\"\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),level,msg)); flush.console() }
log_event("START","Riaz/Gide analyses with historical unresolved Braun slots; exact tests, Firth sensitivity and programme exploration")
d <- fread(unname(cfg_paths["patient_scores"]))
cohorts <- c("Riaz","Gide","Braun")
states <- c("C1","C2","C3","C4")
safe_ci <- function(r,n) if(n>0L) binom.test(r,n)$conf.int else c(NA_real_,NA_real_)
rate <- rbindlist(lapply(cohorts,function(cc) rbindlist(lapply(states,function(st) {
  x <- d[cohort==cc & C_state==st]; x <- x[!is.na(strict_response)]
  n <- nrow(x); r <- sum(x$strict_response=="R")
  ci <- safe_ci(r,n)
  data.table(cohort=cc,C_state=st,n_baseline=if(cc=="Braun") NA_integer_ else sum(d$cohort==cc & d$C_state==st),n_strict=if(cc=="Braun")NA_integer_ else n,R=if(cc=="Braun")NA_integer_ else r,NR=if(cc=="Braun")NA_integer_ else n-r,response_rate=if(n>0)r/n else NA_real_,ci_lo=ci[1],ci_hi=ci[2],status=if(cc=="Braun")"source_unresolved" else "reported")
}))))
fwrite(rate,file.path(out,"B_four_state_response_rates.tsv"),sep="\t",na="NA")
fisher_row <- function(x,ref) {
  a <- x[C_state=="C1"]; b <- x[C_state==ref]
  if(min(nrow(a),nrow(b))<10L) return(list(status="under_10_per_state",or=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_))
  ar <- sum(a$strict_response=="R"); br <- sum(b$strict_response=="R")
  if(ar+br==0L || ar+br==nrow(a)+nrow(b)) return(list(status="single_outcome",or=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_))
  tab <- matrix(c(ar,nrow(a)-ar,br,nrow(b)-br),nrow=2,byrow=TRUE)
  ft <- fisher.test(tab,alternative="two.sided")
  list(status="estimated",or=unname(ft$estimate),lo=ft$conf.int[1],hi=ft$conf.int[2],p=ft$p.value)
}
primary <- list(); pdonly <- list()
for(cc in cohorts) for(ref in c("C4","C3")) {
  x <- d[cohort==cc & C_state %in% c("C1",ref) & !is.na(strict_response)]
  z <- if(cc=="Braun") list(status="source_unresolved",or=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_) else fisher_row(x,ref)
  primary[[length(primary)+1L]] <- data.table(cohort=cc,contrast=paste0("C1/",ref),n_C1=sum(x$C_state=="C1"),n_ref=sum(x$C_state==ref),R_C1=sum(x$C_state=="C1" & x$strict_response=="R"),R_ref=sum(x$C_state==ref & x$strict_response=="R"),status=z$status,OR=z$or,OR_lo=z$lo,OR_hi=z$hi,p=z$p,q=NA_real_)
  if(cc=="Braun") primary[[length(primary)]][,c("n_C1","n_ref","R_C1","R_ref"):=list(NA_integer_,NA_integer_,NA_integer_,NA_integer_)]
  y <- x[raw_response!="SD"]
  h <- if(cc=="Braun") list(status="source_unresolved",or=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_) else fisher_row(y,ref)
  pdonly[[length(pdonly)+1L]] <- data.table(cohort=cc,contrast=paste0("C1/",ref),n_C1=nrow(y[C_state=="C1"]),n_ref=nrow(y[C_state==ref]),status=h$status,OR=h$or,OR_lo=h$lo,OR_hi=h$hi,p=h$p)
  if(cc=="Braun") pdonly[[length(pdonly)]][,c("n_C1","n_ref"):=list(NA_integer_,NA_integer_)]
}
main <- rbindlist(primary); ii <- which(is.finite(main$p)); main$q[ii] <- p.adjust(main$p[ii],method="BH",n=6L)
fwrite(main,file.path(out,"B_primary_exact_comparisons.tsv"),sep="\t",na="NA")
fwrite(rbindlist(pdonly),file.path(out,"B_CRPR_vs_PD_sensitivity.tsv"),sep="\t",na="NA")
# Jeffreys 惩罚似然实现 Firth logistic；CI 为 profile penalized-likelihood CI。
firth <- function(X,y,index) {
  X <- as.matrix(X); storage.mode(X) <- "double"; y <- as.numeric(y)
  if(qr(X)$rank<ncol(X)) stop("Firth design is rank deficient")
  ll <- function(b) {
    eta <- drop(X%*%b); pr <- plogis(eta); w <- pmax(pr*(1-pr),1e-15)
    inf <- crossprod(X*sqrt(w))
    ld <- determinant(inf,logarithm=TRUE)
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
    while(profile_ll(b)>target && step<32) { step<-step*2; b<-bh[index]+side*step }
    if(profile_ll(b)>target) return(if(side<0)-Inf else Inf)
    if(side<0) uniroot(function(z)profile_ll(z)-target,c(b,bh[index]),tol=1e-5)$root
    else uniroot(function(z)profile_ll(z)-target,c(bh[index],b),tol=1e-5)$root
  }
  p <- pchisq(2*(top-profile_ll(0)),df=1,lower.tail=FALSE)
  c(OR=exp(bh[index]),lo=exp(bound(-1)),hi=exp(bound(1)),p=p)
}
# 2x2 无协变量情形须与 Jeffreys 半计数闭式 OR 一致。
self_x <- model.matrix(~state,data.frame(state=c(rep(1,15),rep(0,19))))
self_y <- c(rep(1,10),rep(0,5),rep(1,7),rep(0,12))
self_or <- firth(self_x,self_y,2L)["OR"]
stopifnot(abs(log(self_or/((10.5*12.5)/(5.5*7.5))))<1e-4)
clinical <- list()
for(cc in cohorts) for(ref in c("C4","C3")) {
  x <- d[cohort==cc & C_state %in% c("C1",ref) & !is.na(strict_response)]
  cov <- if(cc=="Riaz") "prior_ipilimumab" else if(cc=="Gide") "regimen" else "trial_source"
  row <- data.table(cohort=cc,contrast=paste0("C1/",ref),covariate=cov,n_main=nrow(x),n_complete=NA_integer_,covariate_completeness=NA_real_,status="source_unresolved",unadjusted_same_case_OR=NA_real_,unadjusted_lo=NA_real_,unadjusted_hi=NA_real_,firth_OR=NA_real_,firth_lo=NA_real_,firth_hi=NA_real_)
  if(cc=="Braun") row[,n_main:=NA_integer_]
  if(cc!="Braun") {
    z <- x[!is.na(get(cov))]; row[,`:=`(n_complete=nrow(z),covariate_completeness=nrow(z)/nrow(x))]
    nlev <- uniqueN(z[[cov]])
    npar <- 1L+as.integer(nlev>1L)
    minor <- min(sum(z$strict_response=="R"),sum(z$strict_response=="NR"))
    status_code <- if(nrow(z)/nrow(x)<.8) "covariate_completeness_under_80pct" else if(minor<npar*5) "minor_outcome_under_5_per_parameter" else if(min(table(z$C_state))<10L) "under_10_per_state" else "estimated"
    row[,status:=status_code]
    if(status_code=="estimated") {
      f <- fisher_row(z,ref)
      row[,`:=`(unadjusted_same_case_OR=f$or,unadjusted_lo=f$lo,unadjusted_hi=f$hi)]
      Z <- data.frame(state=as.integer(z$C_state=="C1"),cov=factor(z[[cov]]))
      X <- if(nlev>1L) model.matrix(~state+cov,Z) else model.matrix(~state,Z)
      ff <- firth(X,as.integer(z$strict_response=="R"),which(colnames(X)=="state"))
      row[,`:=`(firth_OR=ff["OR"],firth_lo=ff["lo"],firth_hi=ff["hi"])]
      if(nlev==1L) row[,covariate:=paste0(cov,"_one_level_omitted")]
    }
  }
  clinical[[length(clinical)+1L]] <- row
}
fwrite(rbindlist(clinical),file.path(out,"B_clinical_Firth_sensitivity.tsv"),sep="\t",na="NA")
continuous <- list()
for(cc in cohorts) for(axis in c("Acid_Score","LM12_inflammatory_axis")) {
  x <- d[cohort==cc]
  row <- data.table(cohort=cc,axis=axis,n_R=sum(x$strict_response=="R"),n_NR=sum(x$strict_response=="NR"),status="source_unresolved",OR_per_SD=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,q=NA_real_)
  if(cc=="Braun") row[,c("n_R","n_NR"):=list(NA_integer_,NA_integer_)]
  if(cc!="Braun") {
    x[,z:=as.numeric(scale(get(axis)))]
    x <- x[!is.na(strict_response)&is.finite(z)]
    if(min(sum(x$strict_response=="R"),sum(x$strict_response=="NR"))>=10L) {
      X <- model.matrix(~z,x)
      ff <- firth(X,as.integer(x$strict_response=="R"),2L)
      row[,`:=`(status="estimated",OR_per_SD=ff["OR"],lo=ff["lo"],hi=ff["hi"],p=ff["p"])]
    } else row[,status:="under_10_per_outcome"]
  }
  continuous[[length(continuous)+1L]] <- row
}
cont <- rbindlist(continuous); ii <- which(is.finite(cont$p)); cont$q[ii] <- p.adjust(cont$p[ii],method="BH",n=6L)
fwrite(cont,file.path(out,"B_continuous_Firth_auxiliary.tsv"),sep="\t",na="NA")
set.seed(20260923)
inside <- list()
rb <- function(a,b) {
  U <- unname(wilcox.test(a,b,exact=FALSE)$statistic)
  2*U/(length(a)*length(b))-1
}
for(cc in cohorts) for(program in c("Hypoxia","ECM","Cytotoxic")) {
  x <- d[cohort==cc & C_state=="C1" & !is.na(strict_response)]
  a <- x[strict_response=="R"][[program]]; b <- x[strict_response=="NR"][[program]]
  row <- data.table(cohort=cc,program=program,n_R=length(a),n_NR=length(b),status="source_unresolved",rank_biserial=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,q=NA_real_)
  if(cc=="Braun") row[,c("n_R","n_NR"):=list(NA_integer_,NA_integer_)]
  if(cc!="Braun") {
    if(min(length(a),length(b))>=10L && all(is.finite(c(a,b)))) {
      boots <- replicate(2000L,rb(sample(a,length(a),replace=TRUE),sample(b,length(b),replace=TRUE)))
      ci <- quantile(boots,c(.025,.975),names=FALSE)
      row[,`:=`(status="estimated",rank_biserial=rb(a,b),lo=ci[1],hi=ci[2],p=wilcox.test(a,b,alternative="two.sided",exact=FALSE)$p.value)]
    } else row[,status:="under_10_per_outcome_or_missing_programme"]
  }
  inside[[length(inside)+1L]] <- row
}
internal <- rbindlist(inside); ii <- which(is.finite(internal$p)); internal$q[ii] <- p.adjust(internal$p[ii],method="BH",n=9L)
fwrite(internal,file.path(out,"B_C1_internal_programme_exploration.tsv"),sep="\t",na="NA")
log_event("DONE",sprintf("primary estimated %d/6; continuous estimated %d/6; C1 internal estimated %d/9",sum(main$status=="estimated"),sum(cont$status=="estimated"),sum(internal$status=="estimated")))
