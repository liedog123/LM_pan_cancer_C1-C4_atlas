# Purpose: subtype context analysis stages.
# Configuration: config/molecular.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Assess programme differences after immune-subtype and BRCA PAM50 adjustment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/14_prepare_subtype_context.R
stage_transcriptome_proteome_14_prepare_subtype_context <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
}

# Purpose: Assess programme differences after immune-subtype and BRCA PAM50 adjustment.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/transcriptome_proteome/15_analyze_subtype_context.R
stage_transcriptome_proteome_15_analyze_subtype_context <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
}
