#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 15_analyze_subtype_context.R --config PATH [--check-config]\n       Rscript 15_analyze_subtype_context.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_mapping", "programme_scores", "output_dir")
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

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table); library(sandwich)})
main <- function() {
  # 普通执行才创建独立输出；源计算连续保留，成功后原子提交。
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  final_dir <- file.path(out_dir, "02_subtype_context_analysis")
  if(file.exists(final_dir) || dir.exists(final_dir)) stop("Output stage already exists.")
  out <- paste0(final_dir, ".tmp_", Sys.getpid())
  if(file.exists(out) || dir.exists(out)) stop("Temporary stage already exists.")
  if(!dir.create(out, showWarnings=FALSE)) stop("Cannot create temporary stage.")
log_event <- function(level,msg) { cat(sprintf("[%s] %s script=15_analyze_subtype_context.R msg=\"%s\"\n",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),level,msg)); flush.console() }
log_event("START","A: conditional state differences with HC3 and fixed BH families")
p <- fread(paths["patient_mapping"])
s <- fread(paths["programme_scores"])
d <- merge(s,p,by=c("Patient_ID","Cancer"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(d)==nrow(s),!anyNA(d$C_subtype))
fwrite(d,file.path(out,"A_analysis_patient_long.tsv"),sep="\t",na="NA")
make_map <- function(subcol,background) {
  x <- p[!is.na(get(subcol)), .N,by=.(Cancer,C_subtype,known_subtype=get(subcol))]
  x[, `:=`(background=background,state_n=sum(N),fraction=N/sum(N)),by=.(Cancer,C_subtype)]
  x
}
fwrite(rbindlist(list(make_map("immune_subtype","immune"),make_map("molecular_subtype","BRCA_mRNA"))),file.path(out,"A_state_known_subtype_counts.tsv"),sep="\t",na="NA")
fit_hc3 <- function(x,adjust=TRUE,purity=FALSE) {
  x <- copy(x); x[, state:=as.integer(C_subtype=="C1")]
  if(purity) x[, purity:=as.numeric(scale(TumorPurity))]
  f <- if(adjust && uniqueN(x$known_subtype)>1L) score ~ state + known_subtype else score ~ state
  if(purity) f <- update(f,. ~ . + purity)
  mod <- lm(f,data=x)
  se <- unname(sqrt(diag(sandwich::vcovHC(mod,type="HC3")))["state"])
  est <- unname(coef(mod)["state"]); df <- df.residual(mod)
  c(est=est,se=se,lo=est-qt(.975,df)*se,hi=est+qt(.975,df)*se,p=2*pt(-abs(est/se),df),df=df)
}
eligible <- function(x) {
  tb <- x[, .N,by=.(known_subtype,C_subtype)]
  wide <- dcast(tb,known_subtype~C_subtype,value.var="N",fill=0)
  if(!all(c("C1",unique(x$C_subtype[x$C_subtype!="C1"])) %in% names(wide))) return(list(x=x[0],reason="one_state_absent",strata=character()))
  ref <- setdiff(unique(x$C_subtype),"C1")
  keep <- wide[get("C1")>=5 & get(ref)>=5]$known_subtype
  xx <- x[known_subtype %in% keep]
  n1 <- sum(xx$C_subtype=="C1"); n0 <- sum(xx$C_subtype==ref)
  df <- nrow(xx)-2L-max(0L,length(keep)-1L)
  why <- if(!length(keep)) "no_shared_stratum_with_5_per_state" else if(min(n1,n0)<10L) "under_10_per_state" else if(df<20L) "residual_df_under_20" else "eligible"
  list(x=xx,reason=why,strata=keep)
}
rows <- list(); desc <- list(); sens <- list()
for (bg in c("immune","BRCA_mRNA")) {
  ccset <- if(bg=="immune") sort(unique(p$Cancer)) else "BRCA"
  subcol <- if(bg=="immune") "immune_subtype" else "molecular_subtype"
  for(cc in ccset) for(ref in c("C4","C3")) for(prog in c("Hypoxia","ECM","Cytotoxic","RPPA_EMT")) {
    x <- d[Cancer==cc & program==prog & C_subtype %in% c("C1",ref) & is.finite(score) & !is.na(get(subcol))]
    x[, known_subtype:=get(subcol)]
    e <- eligible(x); xx <- e$x
    row <- data.table(background=bg,Cancer=cc,contrast=paste0("C1/",ref),program=prog,status=e$reason,
                      n_C1=sum(xx$C_subtype=="C1"),n_ref=sum(xx$C_subtype==ref),
                      n_raw_C1=sum(x$C_subtype=="C1"),n_raw_ref=sum(x$C_subtype==ref),
                      n_strata=length(e$strata),strata=paste(e$strata,collapse="|"),
                      adjusted_effect=NA_real_,adjusted_lo=NA_real_,adjusted_hi=NA_real_,p=NA_real_,q=NA_real_,
                      unadjusted_effect=NA_real_,unadjusted_lo=NA_real_,unadjusted_hi=NA_real_,
                      heterogeneous_direction=NA)
    if(nrow(xx)) {
      st <- xx[, .(n=.N,mean_score=mean(score),median_score=median(score)),by=.(known_subtype,C_subtype)]
      st[, `:=`(background=bg,Cancer=cc,contrast=paste0("C1/",ref),program=prog)]
      desc[[length(desc)+1L]] <- st
    }
    if(e$reason=="eligible") {
      a <- fit_hc3(xx,TRUE); u <- fit_hc3(xx,FALSE)
      row[, `:=`(adjusted_effect=a["est"],adjusted_lo=a["lo"],adjusted_hi=a["hi"],p=a["p"],
                 unadjusted_effect=u["est"],unadjusted_lo=u["lo"],unadjusted_hi=u["hi"])]
      stmean <- dcast(xx[,.(v=mean(score)),by=.(known_subtype,C_subtype)],known_subtype~C_subtype,value.var="v")
      dif <- stmean[["C1"]]-stmean[[ref]]
      row[, heterogeneous_direction:=any(dif>0,na.rm=TRUE)&&any(dif<0,na.rm=TRUE)]
      # 同一纯度完整病例上并列主模型和加纯度模型。
      pur <- xx[is.finite(TumorPurity)]
      ep <- eligible(pur)
      if(ep$reason=="eligible" && sd(ep$x$TumorPurity)>0) {
        base <- fit_hc3(ep$x,TRUE); adj <- fit_hc3(ep$x,TRUE,TRUE)
        sens[[length(sens)+1L]] <- data.table(background=bg,Cancer=cc,contrast=paste0("C1/",ref),program=prog,sensitivity="purity",n=nrow(ep$x),
          base_effect=base["est"],base_lo=base["lo"],base_hi=base["hi"],sensitivity_effect=adj["est"],sensitivity_lo=adj["lo"],sensitivity_hi=adj["hi"],status="eligible")
      } else sens[[length(sens)+1L]] <- data.table(background=bg,Cancer=cc,contrast=paste0("C1/",ref),program=prog,sensitivity="purity",n=nrow(ep$x),status=paste0("not_estimable_",ep$reason))
      if(any(xx$discordant_multi_sample_state,na.rm=TRUE)) {
        ed <- eligible(xx[discordant_multi_sample_state==FALSE])
        if(ed$reason=="eligible") {
          base <- fit_hc3(ed$x,TRUE); sens[[length(sens)+1L]] <- data.table(background=bg,Cancer=cc,contrast=paste0("C1/",ref),program=prog,sensitivity="exclude_discordant",n=nrow(ed$x),sensitivity_effect=base["est"],sensitivity_lo=base["lo"],sensitivity_hi=base["hi"],status="eligible")
        } else sens[[length(sens)+1L]] <- data.table(background=bg,Cancer=cc,contrast=paste0("C1/",ref),program=prog,sensitivity="exclude_discordant",n=nrow(ed$x),status=paste0("not_estimable_",ed$reason))
      }
    }
    rows[[length(rows)+1L]] <- row
  }
}
res <- rbindlist(rows,fill=TRUE)
for(bg in c("immune","BRCA_mRNA")) {
  ii <- which(res$background==bg & is.finite(res$p))
  res$q[ii] <- p.adjust(res$p[ii],method="BH")
}
fwrite(res,file.path(out,"A_conditional_state_contrasts.tsv"),sep="\t",na="NA")
fwrite(rbindlist(desc,fill=TRUE),file.path(out,"A_within_subtype_descriptive.tsv"),sep="\t",na="NA")
fwrite(rbindlist(sens,fill=TRUE),file.path(out,"A_purity_and_discordance_sensitivity.tsv"),sep="\t",na="NA")
log_event("DONE",sprintf("%d contrast slots; %d eligible; %d q<0.05",nrow(res),sum(res$status=="eligible"),sum(res$q<.05,na.rm=TRUE)))

  if(file.exists(final_dir) || dir.exists(final_dir) || !file.rename(out,final_dir)) stop("Atomic output commit failed; temporary files retained: ",out)
  invisible(final_dir)
}
main()
