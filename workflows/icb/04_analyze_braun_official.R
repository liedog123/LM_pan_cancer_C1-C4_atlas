#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 04_analyze_braun_official.R --config PATH [--check-config]\n       Rscript 04_analyze_braun_official.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_scores","member_coverage","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient_scores","member_coverage")]
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
suppressPackageStartupMessages(library(data.table))
# Private project root replaced by explicit DCF inputs.
# Module paths are explicit DCF fields.
out <- unname(cfg_paths["output_dir"])
log_event <- function(level,msg) { cat(sprintf("[%s] %s script=04_analyze_braun_official.R msg=\"%s\"\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),level,gsub('"',"'",msg))); flush.console() }
write_new <- function(x,path) {
  if(file.exists(path) && Sys.getenv("BRAUN_FORCE","false")!="true") stop("Refusing to overwrite existing output: ",path)
  fwrite(x,path,sep="\t",na="NA",quote=FALSE)
}
safe_ci <- function(r,n) if(n>0L) binom.test(r,n)$conf.int else c(NA_real_,NA_real_)

fisher_row <- function(x,ref) {
  a <- x[C_state=="C1"]; b <- x[C_state==ref]
  n1 <- nrow(a); n0 <- nrow(b); r1 <- sum(a$strict_response=="R"); r0 <- sum(b$strict_response=="R")
  if(min(n1,n0)<10L) return(list(status="under_10_per_state",or=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_))
  if(r1+r0==0L || (n1-r1)+(n0-r0)==0L) return(list(status="single_outcome_across_contrast",or=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_))
  ft <- fisher.test(matrix(c(r1,n1-r1,r0,n0-r0),nrow=2,byrow=TRUE),alternative="two.sided")
  list(status="estimated",or=unname(ft$estimate),lo=unname(ft$conf.int[1]),hi=unname(ft$conf.int[2]),p=unname(ft$p.value))
}

# Jeffreys penalized-likelihood logistic regression; confidence limits are profile likelihood.
firth <- function(X,y,index) {
  X <- as.matrix(X); storage.mode(X) <- "double"; y <- as.numeric(y)
  if(qr(X)$rank<ncol(X)) stop("rank_deficient_design")
  ll <- function(b) {
    eta <- drop(X%*%b); pr <- plogis(eta); w <- pmax(pr*(1-pr),1e-15)
    inf <- crossprod(X*sqrt(w)); ld <- determinant(inf,logarithm=TRUE)
    if(ld$sign<=0 || !is.finite(ld$modulus)) return(-1e50)
    sum(y*eta-(pmax(eta,0)+log1p(exp(-abs(eta)))))+.5*as.numeric(ld$modulus)
  }
  full <- optim(rep(0,ncol(X)),function(b)-ll(b),method="BFGS",control=list(maxit=1000,reltol=1e-10))
  if(full$convergence!=0) stop("firth_full_optimization_failed")
  bh <- full$par; top <- ll(bh)
  profile_ll <- function(v) {
    rem <- setdiff(seq_along(bh),index)
    fun <- function(z) { b<-bh; b[index]<-v; if(length(rem)) b[rem]<-z; -ll(b) }
    if(!length(rem)) return(-fun(numeric()))
    op <- optim(bh[rem],fun,method="BFGS",control=list(maxit=1000,reltol=1e-9))
    if(op$convergence!=0) stop("firth_profile_optimization_failed")
    -op$value
  }
  target <- top-qchisq(.95,1)/2
  bound <- function(side) {
    step <- .5; b <- bh[index]+side*step
    while(profile_ll(b)>target && step<32) { step<-step*2; b<-bh[index]+side*step }
    if(profile_ll(b)>target) return(if(side<0)-Inf else Inf)
    if(side<0) uniroot(function(z)profile_ll(z)-target,c(b,bh[index]),tol=1e-5)$root else uniroot(function(z)profile_ll(z)-target,c(bh[index],b),tol=1e-5)$root
  }
  p <- pchisq(2*(top-profile_ll(0)),df=1,lower.tail=FALSE)
  c(OR=exp(bh[index]),lo=exp(bound(-1)),hi=exp(bound(1)),p=p)
}

main <- function() {
  log_event("START","Fixed Braun state-response statistics and bounded auxiliary analyses")
  f <- unname(cfg_paths["patient_scores"])
  if(!file.exists(f)) stop("Patient assignment manifest is missing")
  d <- fread(f)
  stopifnot(nrow(d)==181L,all(d$cohort=="Braun"),sum(d$strict_response=="R",na.rm=TRUE)==39L,sum(d$strict_response=="NR",na.rm=TRUE)==133L,sum(is.na(d$strict_response))==9L)
  cohorts <- "Braun"; states <- c("C1","C2","C3","C4")

  rates <- rbindlist(lapply(states,function(st) {
    x <- d[C_state==st]; y <- x[!is.na(strict_response)]; n <- nrow(y); r <- sum(y$strict_response=="R"); ci <- safe_ci(r,n)
    data.table(cohort="Braun",C_state=st,n_baseline=nrow(x),n_strict=n,R=r,NR=n-r,response_rate=if(n)r/n else NA_real_,ci_lo=ci[1],ci_hi=ci[2],status=if(n)"estimated" else "no_strict_outcomes")
  }))
  write_new(rates,file.path(out,"B_Braun_four_state_response_rates.tsv"))

  primary <- rbindlist(lapply(c("C4","C3"),function(ref) {
    x <- d[C_state %in% c("C1",ref) & !is.na(strict_response)]
    z <- fisher_row(x,ref)
    a <- x[C_state=="C1"]; b <- x[C_state==ref]
    data.table(cohort="Braun",contrast=paste0("C1/",ref),n_C1=nrow(a),n_ref=nrow(b),R_C1=sum(a$strict_response=="R"),R_ref=sum(b$strict_response=="R"),OR=z$or,OR_lo=z$lo,OR_hi=z$hi,p=z$p,q=NA_real_,status=z$status)
  }))
  write_new(primary,file.path(out,"B_Braun_primary_exact_comparisons.tsv"))

  # Sensitivity restricted to objectively responding CR/PR/CRPR versus PD; SD and NE omitted.
  pdrows <- rbindlist(lapply(c("C4","C3"),function(ref) {
    x <- d[C_state %in% c("C1",ref) & raw_response %in% c("CR","PR","CRPR","PD")]
    x[,strict_response_pd:=fcase(raw_response %in% c("CR","PR","CRPR"),"R",raw_response=="PD","NR")]
    x[,strict_response:=strict_response_pd]
    z <- fisher_row(x,ref)
    a <- x[C_state=="C1"]; b <- x[C_state==ref]
    data.table(cohort="Braun",contrast=paste0("C1/",ref),n_C1=nrow(a),n_ref=nrow(b),R_C1=sum(a$strict_response=="R"),R_ref=sum(b$strict_response=="R"),OR=z$or,OR_lo=z$lo,OR_hi=z$hi,p=z$p,status=z$status)
  }))
  write_new(pdrows,file.path(out,"B_Braun_CRPR_vs_PD_sensitivity.tsv"))

  clinical <- rbindlist(lapply(c("C4","C3"),function(ref) {
    cov <- "trial_source"
    x <- d[C_state %in% c("C1",ref) & !is.na(strict_response)]
    z <- x[!is.na(get(cov))]
    nlev <- uniqueN(z[[cov]])
    npar <- 1L+max(0L,nlev-1L)
    minor <- if(nrow(z)) min(sum(z$strict_response=="R"),sum(z$strict_response=="NR")) else 0L
    status <- if(!nrow(x))"no_strict_cases" else if(nrow(z)/nrow(x)<.8)"covariate_completeness_under_80pct" else if(min(nrow(z[C_state=="C1"]),nrow(z[C_state==ref]))<10L)"under_10_per_state" else if(minor<5L*npar)"minor_outcome_under_5_per_parameter" else "estimated"
    same <- list(or=NA_real_,lo=NA_real_,hi=NA_real_,status="not_estimable")
    if(nrow(z) && min(nrow(z[C_state=="C1"]),nrow(z[C_state==ref]))>=10L) same <- fisher_row(z,ref)
    row <- data.table(cohort="Braun",contrast=paste0("C1/",ref),covariate=cov,n_main=nrow(x),n_complete=nrow(z),covariate_completeness=if(nrow(x))nrow(z)/nrow(x) else NA_real_,n_nonintercept_parameters=npar,minor_outcome_n=minor,status=status,unadjusted_same_case_OR=same$or,unadjusted_lo=same$lo,unadjusted_hi=same$hi,firth_OR=NA_real_,firth_lo=NA_real_,firth_hi=NA_real_)
    if(status=="estimated") {
      z[,state:=as.integer(C_state=="C1")]
      z[,trial_source:=droplevels(factor(trial_source))]
      X <- if(nlevels(z$trial_source)>1L) model.matrix(~state+trial_source,z) else model.matrix(~state,z)
      idx <- which(colnames(X)=="state")
      fit <- tryCatch(firth(X,as.integer(z$strict_response=="R"),idx),error=function(e)e)
      if(inherits(fit,"error")) row[,status:=paste0("firth_failed:",conditionMessage(fit))] else row[,`:=`(firth_OR=fit["OR"],firth_lo=fit["lo"],firth_hi=fit["hi"],status=if(nlevels(z$trial_source)==1L)"estimated_trial_source_one_level_omitted" else "estimated")]
    }
    row
  }))
  write_new(clinical,file.path(out,"B_Braun_clinical_Firth_sensitivity.tsv"))

  continuous <- rbindlist(lapply(c("Acid_Score","LM12_inflammatory_axis"),function(axis) {
    allx <- copy(d)
    if(!all(is.finite(allx[[axis]]))) {
      return(data.table(cohort="Braun",axis=axis,n_R=sum(allx$strict_response=="R",na.rm=TRUE),n_NR=sum(allx$strict_response=="NR",na.rm=TRUE),status="missing_axis_values",OR_per_SD=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,q=NA_real_))
    }
    allx[,z:=as.numeric(scale(get(axis)))]
    x <- allx[!is.na(strict_response)&is.finite(z)]
    nr <- sum(x$strict_response=="R"); nnr <- sum(x$strict_response=="NR")
    row <- data.table(cohort="Braun",axis=axis,n_R=nr,n_NR=nnr,status="under_10_per_outcome",OR_per_SD=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,q=NA_real_)
    if(min(nr,nnr)>=10L) {
      X <- model.matrix(~z,x)
      fit <- tryCatch(firth(X,as.integer(x$strict_response=="R"),which(colnames(X)=="z")),error=function(e)e)
      if(inherits(fit,"error")) row[,status:=paste0("firth_failed:",conditionMessage(fit))] else row[,`:=`(status="estimated",OR_per_SD=fit["OR"],lo=fit["lo"],hi=fit["hi"],p=fit["p"])]
    }
    row
  }))
  write_new(continuous,file.path(out,"B_Braun_continuous_Firth_auxiliary.tsv"))

  set.seed(20260923)
  rank_biserial <- function(r, nr) {
    u <- unname(wilcox.test(r,nr,exact=FALSE)$statistic)
    2*u/(length(r)*length(nr))-1
  }
  internal <- rbindlist(lapply(c("Hypoxia","ECM","Cytotoxic"),function(program) {
    x <- d[C_state=="C1" & !is.na(strict_response)]
    n_r <- sum(x$strict_response=="R"); n_nr <- sum(x$strict_response=="NR")
    avail <- is.finite(x[[program]])
    r <- x[strict_response=="R" & avail][[program]]
    nr <- x[strict_response=="NR" & avail][[program]]
    row <- data.table(cohort="Braun",program=program,n_R=n_r,n_NR=n_nr,n_R_scored=length(r),n_NR_scored=length(nr),status="under_10_per_outcome",rank_biserial=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,q=NA_real_)
    cov <- fread(unname(cfg_paths["member_coverage"]))[score==program]
    if(!nrow(cov)||!isTRUE(cov$coverage_pass[1])) row[,status:=paste0("coverage_failed:",if(nrow(cov))cov$missing[1] else "coverage_record_missing")]
    else if(min(n_r,n_nr)<10L) row[,status:="under_10_per_outcome"]
    else if(min(length(r),length(nr))<10L) row[,status:="under_10_scored_per_outcome"]
    else {
      fit <- tryCatch({
        boots <- replicate(2000L,rank_biserial(sample(r,length(r),replace=TRUE),sample(nr,length(nr),replace=TRUE)))
        ci <- quantile(boots,c(.025,.975),names=FALSE)
        list(effect=rank_biserial(r,nr),lo=ci[1],hi=ci[2],p=wilcox.test(r,nr,alternative="two.sided",exact=FALSE)$p.value)
      },error=function(e)e)
      if(inherits(fit,"error")) row[,status:=paste0("test_failed:",conditionMessage(fit))] else row[,`:=`(status="estimated",rank_biserial=fit$effect,lo=fit$lo,hi=fit$hi,p=fit$p)]
    }
    row
  }))
  write_new(internal,file.path(out,"B_Braun_C1_internal_programme_exploration.tsv"))

  log_event("CHECKPOINT",sprintf("primary estimated=%d/2; continuous estimated=%d/2; C1 exploration estimated=%d/3",sum(primary$status=="estimated"),sum(continuous$status=="estimated"),sum(internal$status=="estimated")))
  log_event("DONE","Braun four-state, exact primary, sensitivity, continuous, and C1 exploration tables written")
}
tryCatch(main(),error=function(e){log_event("ERROR",conditionMessage(e));stop(e)})
