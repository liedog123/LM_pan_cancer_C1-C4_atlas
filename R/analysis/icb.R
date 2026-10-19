# Purpose: icb analysis stages.
# Configuration: config/icb.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Assign pretreatment states and evaluate ICB outcomes in Riaz, Gide and Braun.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/01_prepare_riaz_gide.R
stage_icb_01_prepare_riaz_gide <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("lm12_members", "program_membership", "gene_alias", "prior_ipilimumab", "riaz_expression", "riaz_metadata", "gide_expression", "gide_metadata", "braun_hold_metadata")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
cover <- list()
cases <- list()
scale_info <- list()
duplicate_info <- list()
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
}

# Purpose: Assign pretreatment states and evaluate ICB outcomes in Riaz, Gide and Braun.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/02_analyze_riaz_gide.R
stage_icb_02_analyze_riaz_gide <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient_scores")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`safe_ci` <- bind_stage_helper(shared_safe_ci_030, environment())
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
primary <- list()
pdonly <- list()
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
main <- rbindlist(primary)
ii <- which(is.finite(main$p))
main$q[ii] <- p.adjust(main$p[ii],method="BH",n=6L)
fwrite(main,file.path(out,"B_primary_exact_comparisons.tsv"),sep="\t",na="NA")
fwrite(rbindlist(pdonly),file.path(out,"B_CRPR_vs_PD_sensitivity.tsv"),sep="\t",na="NA")
`firth` <- bind_stage_helper(shared_firth_002, environment())
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
cont <- rbindlist(continuous)
ii <- which(is.finite(cont$p))
cont$q[ii] <- p.adjust(cont$p[ii],method="BH",n=6L)
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
internal <- rbindlist(inside)
ii <- which(is.finite(internal$p))
internal$q[ii] <- p.adjust(internal$p[ii],method="BH",n=9L)
fwrite(internal,file.path(out,"B_C1_internal_programme_exploration.tsv"),sep="\t",na="NA")
log_event("DONE",sprintf("primary estimated %d/6; continuous estimated %d/6; C1 internal estimated %d/9",sum(main$status=="estimated"),sum(cont$status=="estimated"),sum(internal$status=="estimated")))
}

# Purpose: Assign pretreatment states and evaluate ICB outcomes in Riaz, Gide and Braun.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/03_prepare_braun_official.R
stage_icb_03_prepare_braun_official <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("official_workbook","lm12_members","program_membership","gene_alias","local_clinical_crosswalk")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`write_new` <- bind_stage_helper(shared_write_new_051, environment())

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
}

# Purpose: Assign pretreatment states and evaluate ICB outcomes in Riaz, Gide and Braun.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/04_analyze_braun_official.R
stage_icb_04_analyze_braun_official <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient_scores","member_coverage")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`write_new` <- bind_stage_helper(shared_write_new_051, environment())
`safe_ci` <- bind_stage_helper(shared_safe_ci_030, environment())

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
}

# Purpose: Assign pretreatment states and evaluate ICB outcomes in Riaz, Gide and Braun.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/05_integrate_three_cohort.R
stage_icb_05_integrate_three_cohort <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
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
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
}

# Purpose: Evaluate the separate post hoc Liu ICB and timing-sensitivity comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/06_prepare_liu_exploratory.R
stage_icb_06_prepare_liu_exploratory <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("source_rda","official_s2_tpm","clinical_workbook","gene_alias","lm12_members","program_membership")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

cfg_paths <- paths
# 普通运行才创建独立空输出；检查模式已退出。
dir.create(unname(cfg_paths["output_dir"]),recursive=TRUE,showWarnings=FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(readxl)
  library(SummarizedExperiment)
})

# Private project root replaced by explicit DCF inputs.
# Module paths are supplied by the checked DCF.
# Prepared source files are independently configured.
out_dir <- unname(cfg_paths["output_dir"])
source_rda <- unname(cfg_paths["source_rda"])
source_tpm <- unname(cfg_paths["official_s2_tpm"])
source_clin <- unname(cfg_paths["clinical_workbook"])

log_event <- function(level, msg) {
  cat(sprintf("[%s] %s script=06_prepare_liu_exploratory.R msg=\"%s\"\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
  flush.console()
}
log_event("START", "tigeR full-gene FPKM conversion and outcome-blind assignment; primary 103, timing sensitivity 119")
required <- c(source_rda, source_tpm, source_clin,
              unname(cfg_paths["gene_alias"]),
              unname(cfg_paths["lm12_members"]),
              unname(cfg_paths["program_membership"]))
if (any(!file.exists(required))) stop("Missing required input: ", paste(required[!file.exists(required)], collapse = " | "))
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L)
  stop("Output directory is non-empty; refusing to overwrite: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
`read_col` <- bind_stage_helper(shared_read_col_034, environment())
patient_id <- trimws(as.character(read_col("A4:A161")[[1]]))
months_pre <- suppressWarnings(as.numeric(read_col("K4:K161")[[1]]))
tx <- trimws(as.character(read_col("R4:R161")[[1]]))
prior_ctla4 <- suppressWarnings(as.integer(read_col("AJ4:AJ161")[[1]]))
ctx_days <- read_col("AP4:AQ161")
biopsy_context <- suppressWarnings(as.integer(ctx_days[[1]]))
days_to_pd1 <- suppressWarnings(as.numeric(ctx_days[[2]]))
site <- trimws(as.character(read_col("J4:J161")[[1]]))
meta_all <- data.table(patient_id = patient_id, monthsBiopsyPreTx = months_pre, Tx = tx,
                       priorCTLA4 = prior_ctla4, biopsyContext = biopsy_context,
                       daysBiopsyToPD1 = days_to_pd1, biopsy_site = site)
meta_all <- meta_all[!is.na(patient_id) & nzchar(patient_id)]
if (anyDuplicated(meta_all$patient_id)) stop("Table S1 patient IDs are not unique")
meta_all[, timepoint_label := fcase(biopsyContext == 1L, "Pre-Ipi", biopsyContext == 2L, "On-Ipi",
                                    biopsyContext == 3L, "Pre-PD1", biopsyContext == 4L, "On-PD1",
                                    default = "Unknown")]

# RDA 中只抽取样本身份两列；结局列不参与本脚本的状态构建。
e <- new.env(parent = globalenv())
loaded <- load(source_rda, envir = e)
if (!"MEL_phs000452" %in% loaded) stop("Author RDA did not contain MEL_phs000452")
obj <- e$MEL_phs000452
fpkm0 <- as.matrix(SummarizedExperiment::assay(obj))
storage.mode(fpkm0) <- "double"
cd <- as.data.frame(SummarizedExperiment::colData(obj)[, c("sample_id", "patient_name"), drop = FALSE])
if (!all(c("sample_id", "patient_name") %in% names(cd))) stop("tigeR colData is missing sample_id/patient_name")
source_map <- data.table(sample_id = as.character(cd$sample_id), patient_name = as.character(cd$patient_name))
if (!setequal(colnames(fpkm0), source_map$sample_id)) stop("FPKM assay columns do not match tigeR sample_id")
source_map[, patient_id := sub("_.*$", "", patient_name)]
source_map <- source_map[grepl("^Patient[0-9]+$", patient_id)]
if (anyDuplicated(source_map$sample_id) || anyDuplicated(source_map$patient_id))
  stop("tigeR Patient### keys are not one-to-one with samples")
if (nrow(fpkm0) != 56269L || nrow(source_map) != 121L) stop("Unexpected tigeR object dimensions or Patient### count")
if (any(!is.finite(fpkm0)) || any(fpkm0 < 0)) stop("tigeR FPKM contains non-finite or negative values")
if (!setequal(rownames(fpkm0), rownames(obj))) stop("Assay gene rownames mismatch")

# 官方 S2 只用于 Patient ID 交集和表达量尺度对照，不作为状态评分矩阵。
official <- fread(source_tpm, check.names = FALSE, showProgress = FALSE)
if (nrow(official) != 121L || ncol(official) != 20849L) stop("Unexpected official S2 dimensions")
official_ids <- as.character(official[[1]])
official_gene <- toupper(trimws(names(official)[-1]))
if (anyDuplicated(official_ids) || anyDuplicated(official_gene)) stop("Duplicate official S2 sample or gene IDs")
official_tpm <- as.matrix(official[, -1, with = FALSE])
storage.mode(official_tpm) <- "double"
rownames(official_tpm) <- official_ids
colnames(official_tpm) <- official_gene
if (any(!is.finite(official_tpm)) || any(official_tpm < 0)) stop("Official S2 TPM contains invalid values")

# 完整 tigeR FPKM 按每样本全部保留基因换算 TPM；不按 S2 过滤后的 20,848 行再次缩放。
gene_raw <- toupper(trimws(rownames(fpkm0)))
if (anyNA(gene_raw) || any(!nzchar(gene_raw)) || anyDuplicated(gene_raw)) stop("tigeR gene symbols are missing or duplicated")
rownames(fpkm0) <- gene_raw
fpkm_sum_raw <- colSums(fpkm0)
if (any(fpkm_sum_raw <= 0)) stop("tigeR has a zero-sum FPKM sample")
tpm_for_source_check <- sweep(fpkm0, 2, fpkm_sum_raw, "/") * 1e6
common_gene <- intersect(rownames(tpm_for_source_check), colnames(official_tpm))
common_patient <- sort(intersect(source_map$patient_id, official_ids))
if (length(common_patient) != 120L || length(common_gene) < 10000L)
  stop(sprintf("Unexpected source concordance dimensions: patients=%d genes=%d", length(common_patient), length(common_gene)))
cor_row <- function(local_id, official_id, comparison) {
  a <- log2(tpm_for_source_check[common_gene, source_map$sample_id[match(local_id, source_map$patient_id)]] + 1)
  b <- log2(official_tpm[official_id, common_gene] + 1)
  data.table(comparison = comparison, local_patient_id = local_id, official_S2_patient_id = official_id,
             shared_gene_n = length(common_gene), pearson_log2TPM = cor(a, b, method = "pearson"),
             spearman_log2TPM = cor(a, b, method = "spearman"))
}
matched_cor <- rbindlist(lapply(common_patient, function(id) cor_row(id, id, "exact_patient_key")))
shift_ids <- c(common_patient[-1L], common_patient[1L])
shift_cor <- rbindlist(lapply(seq_along(common_patient), function(i)
  cor_row(common_patient[i], shift_ids[i], "fixed_one_position_shift_control")))
patient192_test <- if ("Patient192" %in% source_map$patient_id && "Patient20" %in% official_ids)
  cor_row("Patient192", "Patient20", "unmapped_candidate_Patient192_vs_Patient20") else
  data.table(comparison = "unmapped_candidate_Patient192_vs_Patient20", local_patient_id = "Patient192",
             official_S2_patient_id = "Patient20", shared_gene_n = length(common_gene),
             pearson_log2TPM = NA_real_, spearman_log2TPM = NA_real_)
concordance <- rbindlist(list(matched_cor, shift_cor, patient192_test), fill = TRUE)
source_crosswalk <- data.table(official_S2_patient_id = sort(official_ids))
source_crosswalk[, exact_tigeR_sample_match := source_map$sample_id[match(official_S2_patient_id, source_map$patient_id)]]
source_crosswalk[, exact_tigeR_patient_key_match := !is.na(exact_tigeR_sample_match)]
source_crosswalk[, mapping_status := fifelse(exact_tigeR_patient_key_match, "exact Patient### prefix match", "no exact tigeR Patient### key")]
source_extra <- source_map[!patient_id %in% official_ids,
                           .(official_S2_patient_id = NA_character_, exact_tigeR_sample_match = sample_id,
                             exact_tigeR_patient_key_match = FALSE, mapping_status = paste0("tigeR-only exact key: ", patient_id))]
source_crosswalk <- rbindlist(list(source_crosswalk, source_extra), fill = TRUE)
if (sum(source_crosswalk$exact_tigeR_patient_key_match) != 120L) stop("Expected 120 exact tigeR/S2 Patient### matches")
if (any(source_crosswalk$official_S2_patient_id == "Patient20" & source_crosswalk$exact_tigeR_patient_key_match)) stop("Patient20 unexpectedly mapped")

# 临床主集合与宽时间窗先按 Table S1 定义，再要求 tigeR 有精确患者键。
meta_s2 <- meta_all[patient_id %in% official_ids]
if (any(meta_s2$biopsyContext == 3L & meta_s2$daysBiopsyToPD1 > 0, na.rm = TRUE))
  stop("Positive day offset conflicts with Table S1 Pre-PD1 context")
if (any(meta_s2$biopsyContext == 4L & meta_s2$daysBiopsyToPD1 <= 0, na.rm = TRUE))
  stop("On-PD1 context conflicts with nonpositive day offset")
sample_map <- source_map[, .(patient_id, sample_id, patient_name)]
manifest <- merge(meta_all, sample_map, by = "patient_id", all.x = TRUE, sort = FALSE)
manifest[, has_exact_tigeR_patient_key := !is.na(sample_id)]
manifest[, has_official_S2_TPM_patient := patient_id %in% official_ids]
manifest[, official_primary_timepoint := has_official_S2_TPM_patient & !is.na(biopsyContext) & biopsyContext == 3L & !is.na(daysBiopsyToPD1) & daysBiopsyToPD1 <= 0]
manifest[, official_broad_timepoint := has_official_S2_TPM_patient & !is.na(biopsyContext) & biopsyContext %in% c(1L, 2L, 3L) & !is.na(daysBiopsyToPD1) & daysBiopsyToPD1 <= 0]
manifest[, include_primary103 := official_primary_timepoint & has_exact_tigeR_patient_key]
manifest[, include_broad119 := official_broad_timepoint & has_exact_tigeR_patient_key]
manifest[, primary_reason := fcase(!has_official_S2_TPM_patient, "excluded_not_in_official_S2_expression_cohort",
                                    !official_primary_timepoint, "not_Pre-PD1_context3_day_le_0",
                                    !has_exact_tigeR_patient_key, "excluded_no_exact_tigeR_patient_key",
                                    default = "included_primary_Pre-PD1_context3")]
manifest[, broad_reason := fcase(!has_official_S2_TPM_patient, "excluded_not_in_official_S2_expression_cohort",
                                  !official_broad_timepoint, "not_context1to3_day_le_0",
                                  !has_exact_tigeR_patient_key, "excluded_no_exact_tigeR_patient_key",
                                  default = "included_broad_pre-PD1_context1to3")]
manifest[, day0_prePD1_context := biopsyContext == 3L & daysBiopsyToPD1 == 0]
primary_official <- manifest[official_primary_timepoint == TRUE]
broad_official <- manifest[official_broad_timepoint == TRUE]
primary_meta <- manifest[include_primary103 == TRUE]
broad_meta <- manifest[include_broad119 == TRUE]
if (nrow(primary_official) != 104L || nrow(broad_official) != 120L || nrow(primary_meta) != 103L || nrow(broad_meta) != 119L)
  stop(sprintf("Unexpected eligible/mapped set sizes: official primary=%d broad=%d; mapped primary=%d broad=%d",
               nrow(primary_official), nrow(broad_official), nrow(primary_meta), nrow(broad_meta)))
if (!setequal(primary_official$patient_id[!primary_official$has_exact_tigeR_patient_key], "Patient20"))
  stop("The sole unmapped primary patient is not Patient20")
if (sum(primary_meta$daysBiopsyToPD1 == 0, na.rm = TRUE) != 2L) stop("Unexpected mapped day-zero primary count")

alias <- fread(unname(cfg_paths["gene_alias"]))
alias_map <- setNames(toupper(alias$canonical_symbol), toupper(alias$alias_symbol))
normalize_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- ""
  hit <- x %in% names(alias_map)
  x[hit] <- unname(alias_map[x[hit]])
  x
}
gene <- normalize_symbol(rownames(fpkm0))
keep_gene <- nzchar(gene)
fpkm <- fpkm0[keep_gene, , drop = FALSE]
rownames(fpkm) <- gene[keep_gene]
duplicate_gene_rows <- sum(duplicated(rownames(fpkm)))
if (duplicate_gene_rows > 0L) {
  n_gene_rows <- rowsum(matrix(1, nrow(fpkm), ncol(fpkm)), group = rownames(fpkm), reorder = FALSE)
  fpkm <- rowsum(fpkm, group = rownames(fpkm), reorder = FALSE) / n_gene_rows
}
fpkm_sum <- colSums(fpkm)
if (any(fpkm_sum <= 0)) stop("Zero-sum FPKM sample after gene-symbol normalization")
tpm <- sweep(fpkm, 2, fpkm_sum, "/") * 1e6
logexpr <- log2(tpm + 1)

lm_long <- fread(unname(cfg_paths["lm12_members"]))
lm_sets <- lapply(split(normalize_symbol(lm_long$gene), lm_long$feature_id), unique)
lm_sets <- lm_sets[sort(names(lm_sets))]
if (length(lm_sets) != 12L || !"LM12_Antimicrobial_Peptides" %in% names(lm_sets)) stop("Unexpected frozen LM12 signature catalog")
acid <- c("ATP6AP1", "GAPDH", "HIF1A", "HK2", "PFKP", "PGK1", "PKM", "SLC2A1", "SLC9A7")
pd <- fread(unname(cfg_paths["program_membership"]))
prog_ids <- c(Hypoxia = "MALIGNANT_HYPOXIA", ECM = "STROMAL_ECM_ORGANIZATION", Cytotoxic = "TNK_CYTOTOXICITY_EFFECTOR")
programs <- lapply(prog_ids, function(id) sort(unique(normalize_symbol(pd[program_id == id & include_in_definition_independent_score == TRUE]$gene))))
if (length(programs$Cytotoxic) != 6L) stop("Frozen cytotoxic program does not contain the expected six members")

assign_one_set <- function(set_name, md) {
  setorder(md, patient_id)
  ids <- md$patient_id
  sample_ids <- md$sample_id
  expr <- logexpr[, sample_ids, drop = FALSE]
  colnames(expr) <- ids
  sets <- c(list(Acid_Score = acid), lm_sets)
  coverage <- list()
  present_sets <- list()
  for (nm in names(sets)) {
    gs <- unique(sets[[nm]])
    seen <- intersect(gs, rownames(expr))
    vv <- expr[seen, , drop = FALSE]
    sdv <- if (nrow(vv) > 0L) apply(vv, 1, sd) else numeric()
    variable <- names(sdv)[is.finite(sdv) & sdv > 0]
    pass <- if (nm == "Acid_Score") length(seen) == 9L && length(variable) == 9L else length(seen) / length(gs) >= 0.70 && length(variable) > 0L
    coverage[[length(coverage) + 1L]] <- data.table(sample_set = set_name, score = nm, frozen_n = length(gs),
      present_n = length(seen), variable_n = length(variable), coverage_pass = pass,
      missing = paste(setdiff(gs, seen), collapse = "|"), zero_variance = paste(setdiff(seen, variable), collapse = "|"))
    present_sets[[nm]] <- seen
  }
  cov <- rbindlist(coverage)
  if (!all(cov$coverage_pass)) stop(set_name, ": fixed Acid/LM12 gene coverage failed")
  ss <- GSVA::gsva(GSVA::ssgseaParam(exprData = expr, geneSets = present_sets, normalize = TRUE), verbose = FALSE)
  lm <- ss[names(lm_sets), , drop = FALSE]
  axis <- colMeans(lm[setdiff(rownames(lm), "LM12_Antimicrobial_Peptides"), , drop = FALSE])
  acid_score <- as.numeric(ss["Acid_Score", ]); names(acid_score) <- colnames(ss)

  # 使用已冻结的九基因 z-score + ward.D/k=3 和 LM12 ssGSEA + ward.D/k=2。
  az <- expr[acid, , drop = FALSE]
  az <- sweep(sweep(az, 1, rowMeans(az), "-"), 1, apply(az, 1, sd), "/")
  ca <- cutree(hclust(dist(t(az)), method = "ward.D"), k = 3)
  means_a <- tapply(acid_score, ca, mean)
  if (anyDuplicated(as.numeric(means_a))) stop(set_name, ": Acid cluster means tie; AG direction unresolved")
  ag_map <- setNames(paste0("AG", 1:3), names(sort(means_a, decreasing = TRUE)))
  ci <- cutree(hclust(dist(t(lm)), method = "ward.D"), k = 2)
  means_i <- tapply(axis, ci, mean)
  if (anyDuplicated(as.numeric(means_i))) stop(set_name, ": LM12 cluster means tie; IG direction unresolved")
  ig_map <- setNames(paste0("IG", 1:2), names(sort(means_i, decreasing = TRUE)))
  ans <- data.table(patient_id = ids, sample_id = sample_ids, AG = unname(ag_map[as.character(ca)]),
    IG = unname(ig_map[as.character(ci)]), Acid_Score = acid_score[ids], LM12_inflammatory_axis = axis[ids])
  ans[, C_state := fcase(AG == "AG1" & IG == "IG1", "C1", AG == "AG3" & IG == "IG2", "C4",
                         AG == "AG1" & IG == "IG2", "C2", AG %in% c("AG2", "AG3") & IG == "IG1", "C3",
                         AG == "AG2" & IG == "IG2", "C2")]
  if (anyNA(ans$C_state)) stop(set_name, ": C1-C4 map incomplete")
  for (nm in names(programs)) {
    gs <- unique(programs[[nm]])
    seen <- intersect(gs, rownames(expr))
    vv <- expr[seen, , drop = FALSE]
    sdv <- if (nrow(vv) > 0L) apply(vv, 1, sd) else numeric()
    variable <- names(sdv)[is.finite(sdv) & sdv > 0]
    pass <- if (nm == "Cytotoxic") length(seen) == 6L && length(variable) >= 5L else length(seen) / length(gs) >= 0.80 && length(variable) >= 10L
    cov <- rbind(cov, data.table(sample_set = set_name, score = nm, frozen_n = length(gs), present_n = length(seen),
      variable_n = length(variable), coverage_pass = pass, missing = paste(setdiff(gs, seen), collapse = "|"),
      zero_variance = paste(setdiff(seen, variable), collapse = "|")))
    val <- rep(NA_real_, length(ids)); names(val) <- ids
    if (pass) {
      use <- vv[variable, , drop = FALSE]
      val <- colMeans(sweep(sweep(use, 1, rowMeans(use), "-"), 1, sdv[variable], "/"))
      if (!is.finite(sd(val)) || sd(val) == 0) pass <- FALSE
      else {
        val <- as.numeric(scale(val))
        names(val) <- ids
      }
    }
    if (!pass) val[] <- NA_real_
    ans[, (nm) := val[ids]]
  }
  scored <- merge(ans, md, by = c("patient_id", "sample_id"), all.x = TRUE, sort = FALSE)
  setorder(scored, patient_id)
  scored[, sample_set := set_name]
  list(data = scored, coverage = cov)
}

res103 <- assign_one_set("primary_Pre_PD1_context3_n103", primary_meta)
res119 <- assign_one_set("sensitivity_all_pre_PD1_context1to3_n119", broad_meta)
state_stability <- merge(res103$data[, .(patient_id, primary_state = C_state)],
                         res119$data[, .(patient_id, broad_state = C_state)], by = "patient_id")
state_table <- as.data.table(as.table(table(factor(state_stability$primary_state, levels = c("C1", "C2", "C3", "C4")),
                                             factor(state_stability$broad_state, levels = c("C1", "C2", "C3", "C4")))))
setnames(state_table, c("primary_state", "broad_state", "n"))
stability_summary <- data.table(n_shared = nrow(state_stability), n_same_state = sum(state_stability$primary_state == state_stability$broad_state),
  exact_agreement = mean(state_stability$primary_state == state_stability$broad_state),
  n_primary_day0 = sum(primary_meta$daysBiopsyToPD1 == 0, na.rm = TRUE),
  n_primary_day_negative = sum(primary_meta$daysBiopsyToPD1 < 0, na.rm = TRUE),
  n_broad_Pre_Ipi = sum(broad_meta$biopsyContext == 1L), n_broad_On_Ipi = sum(broad_meta$biopsyContext == 2L),
  n_broad_Pre_PD1 = sum(broad_meta$biopsyContext == 3L))

audit <- rbindlist(list(
  data.table(metric = c("tigeR_source_object", "source_expression_scale", "analysis_expression_scale",
                        "official_S2_use", "TPM_scaling_universe", "official_S2_exact_patient_matches",
                        "official_S2_patients_without_exact_tigeR_match", "tigeR_only_patient_key",
                        "shared_exact_gene_symbols_for_concordance", "matched_log2TPM_Pearson_median",
                        "matched_log2TPM_Pearson_P10", "matched_log2TPM_Spearman_median",
                        "shift_control_log2TPM_Pearson_median", "Patient192_vs_Patient20_log2TPM_Pearson",
                        "Table_S1_primary_context3_day_le0", "primary_103_actual", "Table_S1_broad_context1to3_day_le0",
                        "broad_119_actual", "primary_day0_official", "primary_day0_actual", "duplicate_gene_rows_after_alias",
                        "GSVA_version", "ssGSEA_normalize"),
             value = c("YuLab-SMU/tigeR.data MEL_phs000452.rda", "FPKM", "FPKM/sum(FPKM over full 56,269-gene object)*1e6 then log2(TPM+1)",
                       "patient ID and expression-scale concordance only; not used for state scores", "all tigeR gene rows after frozen symbol alias aggregation by mean",
                       sum(source_crosswalk$exact_tigeR_patient_key_match), "Patient20", "Patient192", length(common_gene),
                       median(matched_cor$pearson_log2TPM), unname(quantile(matched_cor$pearson_log2TPM, .10)),
                       median(matched_cor$spearman_log2TPM), median(shift_cor$pearson_log2TPM), patient192_test$pearson_log2TPM[1],
                       nrow(primary_official), nrow(primary_meta), nrow(broad_official), nrow(broad_meta),
                       sum(primary_official$daysBiopsyToPD1 == 0, na.rm = TRUE), sum(primary_meta$daysBiopsyToPD1 == 0, na.rm = TRUE),
                       duplicate_gene_rows, as.character(packageVersion("GSVA")), "TRUE"))), fill = TRUE)
tx_summary <- primary_meta[, .N, by = Tx][order(Tx)]
tx_summary[, analysis_set := "primary_103"]
prior_summary <- primary_meta[, .N, by = priorCTLA4][order(priorCTLA4)]
prior_summary[, analysis_set := "primary_103"]

fwrite(source_crosswalk, file.path(out_dir, "Liu_tigeR_S2_exact_patient_crosswalk.tsv"), sep = "\t", na = "NA")
fwrite(concordance, file.path(out_dir, "Liu_tigeR_S2_expression_concordance.tsv"), sep = "\t", na = "NA")
fwrite(manifest, file.path(out_dir, "Liu_103_119_patient_inclusion_manifest.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(list(res103$data, res119$data), fill = TRUE), file.path(out_dir, "Liu_103_119_outcome_blind_state_assignments.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(list(res103$coverage, res119$coverage)), file.path(out_dir, "Liu_103_119_member_coverage.tsv"), sep = "\t", na = "NA")
fwrite(state_table, file.path(out_dir, "Liu_103_vs_119_assignment_crosstab.tsv"), sep = "\t", na = "NA")
fwrite(stability_summary, file.path(out_dir, "Liu_103_vs_119_assignment_stability.tsv"), sep = "\t", na = "NA")
fwrite(audit, file.path(out_dir, "Liu_tigeR_input_and_source_audit.tsv"), sep = "\t", na = "NA")
fwrite(tx_summary, file.path(out_dir, "Liu_103_Tx_levels.tsv"), sep = "\t", na = "NA")
fwrite(prior_summary, file.path(out_dir, "Liu_103_priorCTLA4_levels.tsv"), sep = "\t", na = "NA")

log_event("CHECKPOINT", sprintf("source mapping: 120/121 exact Patient### IDs; Patient20 has no exact expression key; Patient192 remains unmapped"))
log_event("DONE", sprintf("outcome-blind assignment complete: primary=%d/%d; broad=%d/%d; shared state agreement=%.3f; output=%s",
                           nrow(res103$data), nrow(primary_official), nrow(res119$data), nrow(broad_official),
                           stability_summary$exact_agreement, out_dir))
}

# Purpose: Evaluate the separate post hoc Liu ICB and timing-sensitivity comparisons.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/07_analyze_liu_exploratory.R
stage_icb_07_analyze_liu_exploratory <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("clinical_workbook","outcome_blind_assignments")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`read_col` <- bind_stage_helper(shared_read_col_034, environment())
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
ii <- which(is.finite(main_cmp$p))
if(length(ii)) main_cmp$q[ii] <- p.adjust(main_cmp$p[ii],method="BH",n=2L)
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
`firth_fit` <- bind_stage_helper(shared_firth_002, environment())

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
ii <- which(is.finite(continuous$p))
if(length(ii)) continuous$q[ii] <- p.adjust(continuous$p[ii],method="BH",n=2L)
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
ii <- which(is.finite(internal$p))
if(length(ii)) internal$q[ii] <- p.adjust(internal$p[ii],method="BH",n=3L)
fwrite(internal,file.path(res_dir,"Liu_primary103_C1_internal_programme_exploration.tsv"),sep="\t",na="NA")

log_event("DONE",sprintf("Liu primary n=103; exact slots %d/2; continuous %d/2; C1 programmes %d/3; broad timing sensitivity n=119",sum(main_cmp$status=="estimated"),sum(continuous$status=="estimated"),sum(internal$status=="estimated")))
}

# Purpose: Evaluate the exploratory within-C1 Gide programme-response associations.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/icb/08_run_gide_c1_sensitivity.R
stage_icb_08_run_gide_c1_sensitivity <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient_scores","frozen_program_members","lm12_members")]
required_files <- unname(input_paths)

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
input_dirs <- character(0L)
if(any(!dir.exists(input_dirs)) || any(file.access(input_dirs,5L)!=0L)) stop("Input directories must be readable and searchable.")
package_root <- context$root
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

cfg_paths <- paths
# 普通运行才创建独立空输出；检查模式已退出。
dir.create(unname(cfg_paths["output_dir"]),recursive=TRUE,showWarnings=FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# Private project root replaced by explicit DCF inputs.
# Module paths are supplied by the checked DCF.
src <- unname(cfg_paths["patient_scores"])
out <- unname(cfg_paths["output_dir"])
res <- file.path(out, "results")
fig <- file.path(out, "figures")
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

log_event <- function(level, msg) {
  cat(sprintf("[%s] %s script=08_run_gide_c1_sensitivity.R msg=\"%s\"\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
  flush.console()
}

log_event("START", paste0(
  "post-hoc sensitivity in Gide C1; input=", src,
  "; target=regimen distribution and Cytotoxic-LM12 axis relation",
  "; Firth profile 95% CI; Spearman patient bootstrap B=5000 seed=20260923",
  "; outputs=", out))

d <- fread(src)
x <- d[cohort == "Gide" & C_state == "C1" & !is.na(strict_response)]
if (nrow(x) != 25L || uniqueN(x$patient_id) != 25L || uniqueN(x$sample_id) != 25L ||
    sum(x$strict_response == "R") != 13L || sum(x$strict_response == "NR") != 12L) {
  stop("Gide C1 differs from the fixed 25-patient input summary; review inclusion before analysis")
}
if (anyNA(x$regimen) || any(!x$strict_response %in% c("R", "NR")) ||
    any(!is.finite(x$Cytotoxic)) || any(!is.finite(x$LM12_inflammatory_axis))) {
  stop("Gide C1 has missing or invalid values in required outcome, regimen, or score fields")
}

# 沿原始治疗标签归并为单药/联合；保留原始 regimen 供核查。
regimen_lower <- tolower(trimws(x$regimen))
is_combo <- grepl("ctla|ipilimumab|combination|combined", regimen_lower)
is_mono <- grepl("pd.?1|nivolumab|pembrolizumab", regimen_lower) & !is_combo
if (any(!is_combo & !is_mono)) {
  stop(paste0("Unresolved Gide C1 regimen label(s): ",
              paste(unique(x$regimen[!is_combo & !is_mono]), collapse = " | ")))
}
x[, regimen_group := ifelse(is_combo, "Combination", "Monotherapy")]
x[, regimen_group := factor(regimen_group, levels = c("Monotherapy", "Combination"))]
x[, response_label := factor(strict_response, levels = c("NR", "R"),
                              labels = c("SD/PD", "CR/PR"))]

# 列出每种方案内响应率，并同时显示每种方案占响应者/非响应者的比例。
regimen_summary <- x[, .(
  n = .N,
  n_R = sum(strict_response == "R"),
  n_NR = sum(strict_response == "NR"),
  regimen_source = paste(sort(unique(regimen)), collapse = " | ")
), by = .(regimen_group)]
regimen_summary[, response_rate := n_R / n]
regimen_summary[, `:=`(
  response_ci_low = mapply(function(r, n) binom.test(r, n)$conf.int[1], n_R, n),
  response_ci_high = mapply(function(r, n) binom.test(r, n)$conf.int[2], n_R, n),
  share_among_R = n_R / sum(x$strict_response == "R"),
  share_among_NR = n_NR / sum(x$strict_response == "NR")
)]
setorder(regimen_summary, regimen_group)
fwrite(regimen_summary, file.path(res, "Gide_C1_regimen_response_distribution.tsv"),
       sep = "\t", na = "NA")

# 只作小样本描述：展示各治疗/响应格子的冻结 Cytotoxic 分数分布。
cell_grid <- CJ(regimen_group = levels(x$regimen_group),
                strict_response = c("R", "NR"), unique = TRUE)
cyt_summary <- x[, .(
  n = .N,
  mean = mean(Cytotoxic),
  median = median(Cytotoxic),
  min = min(Cytotoxic),
  max = max(Cytotoxic),
  q1 = unname(quantile(Cytotoxic, .25)),
  q3 = unname(quantile(Cytotoxic, .75))
), by = .(regimen_group = as.character(regimen_group), strict_response)]
cyt_summary <- merge(cell_grid, cyt_summary,
                     by = c("regimen_group", "strict_response"), all.x = TRUE, sort = FALSE)
cyt_summary[is.na(n), n := 0L]
setorder(cyt_summary, regimen_group, strict_response)
fwrite(cyt_summary, file.path(res, "Gide_C1_cytotoxic_by_regimen_response.tsv"),
       sep = "\t", na = "NA")

# 在同一完整病例中比较原始与 regimen 调整的 Firth logistic OR/区间，不作第二套显著性筛选。
cc <- x[complete.cases(x[, .(strict_response, regimen_group, Cytotoxic)])]
cc[, response_binary := as.integer(strict_response == "R")]
cc[, Cytotoxic_cohort_SD := Cytotoxic]

firth_profile_ci <- function(X, y, index, profile = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)
  if (qr(X)$rank < ncol(X)) stop("Firth design matrix is rank deficient")
  penalized_loglik <- function(beta) {
    eta <- drop(X %*% beta)
    prob <- plogis(eta)
    weight <- pmax(prob * (1 - prob), 1e-15)
    info <- crossprod(X * sqrt(weight))
    det_info <- determinant(info, logarithm = TRUE)
    if (det_info$sign <= 0 || !is.finite(det_info$modulus)) return(-1e50)
    loglik <- sum(y * eta - ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta))))
    loglik + 0.5 * as.numeric(det_info$modulus)
  }
  full <- optim(rep(0, ncol(X)), function(beta) -penalized_loglik(beta),
                method = "BFGS", control = list(maxit = 1500, reltol = 1e-10))
  if (full$convergence != 0L) stop("Firth full model optimization failed")
  beta_hat <- full$par
  if (!profile) return(c(OR = exp(beta_hat[index]), logOR = beta_hat[index]))
  ll_top <- penalized_loglik(beta_hat)
  remaining <- setdiff(seq_along(beta_hat), index)
  profile_ll <- function(value) {
    beta <- beta_hat
    beta[index] <- value
    if (!length(remaining)) return(penalized_loglik(beta))
    opt <- optim(beta_hat[remaining], function(rest) {
      candidate <- beta
      candidate[remaining] <- rest
      -penalized_loglik(candidate)
    }, method = "BFGS", control = list(maxit = 1500, reltol = 1e-9))
    if (opt$convergence != 0L) stop("Firth profile optimization failed")
    -opt$value
  }
  target <- ll_top - qchisq(.95, df = 1) / 2
  find_bound <- function(side) {
    step <- .5
    bound <- beta_hat[index] + side * step
    while (profile_ll(bound) > target && step < 32) {
      step <- step * 2
      bound <- beta_hat[index] + side * step
    }
    if (profile_ll(bound) > target) return(if (side < 0) -Inf else Inf)
    interval <- if (side < 0) c(bound, beta_hat[index]) else c(beta_hat[index], bound)
    uniroot(function(value) profile_ll(value) - target, interval, tol = 1e-5)$root
  }
  c(OR = exp(beta_hat[index]),
    CI_low = exp(find_bound(-1)),
    CI_high = exp(find_bound(1)))
}

fit_one <- function(dat, adjusted) {
  if (adjusted && uniqueN(dat$regimen_group) > 1L) {
    design_data <- data.frame(Cytotoxic_cohort_SD = dat$Cytotoxic_cohort_SD,
                              regimen_group = droplevels(dat$regimen_group))
    X <- model.matrix(~ Cytotoxic_cohort_SD + regimen_group, design_data)
    term <- which(colnames(X) == "Cytotoxic_cohort_SD")
    adjustment <- "regimen group"
  } else {
    design_data <- data.frame(Cytotoxic_cohort_SD = dat$Cytotoxic_cohort_SD)
    X <- model.matrix(~ Cytotoxic_cohort_SD, design_data)
    term <- which(colnames(X) == "Cytotoxic_cohort_SD")
    adjustment <- if (adjusted) "regimen has one level; omitted" else "none"
  }
  est <- firth_profile_ci(X, dat$response_binary, term)
  c1_sd <- sd(dat$Cytotoxic_cohort_SD)
  data.table(model = if (adjusted) "Regimen-adjusted" else "Unadjusted",
             adjustment = adjustment, n_complete = nrow(dat),
             n_R = sum(dat$response_binary == 1), n_NR = sum(dat$response_binary == 0),
             score_scale = "one SD in all Gide baseline patients",
             OR_per_SD = est["OR"], CI_low = est["CI_low"], CI_high = est["CI_high"],
             C1_25_patient_SD_in_cohort_SD_units = c1_sd,
             OR_per_C1_25_patient_SD = exp(log(est["OR"]) * c1_sd),
             CI_low_per_C1_25_patient_SD = est["CI_low"] ^ c1_sd,
             CI_high_per_C1_25_patient_SD = est["CI_high"] ^ c1_sd,
             interval = "95% profile penalized-likelihood CI")
}
firth_rows <- rbindlist(list(fit_one(cc, FALSE), fit_one(cc, TRUE)))
fwrite(firth_rows, file.path(res, "Gide_C1_cytotoxic_Firth_regimen_sensitivity.tsv"),
       sep = "\t", na = "NA")

# 散点中可见高 Cytotoxic 响应者；逐一剔除患者检查 regimen 调整点估计是否由单人支配。
loo_rows <- rbindlist(lapply(seq_len(nrow(cc)), function(i) {
  keep <- cc[-i]
  design_data <- data.frame(Cytotoxic_cohort_SD = keep$Cytotoxic_cohort_SD,
                            regimen_group = droplevels(keep$regimen_group))
  X <- model.matrix(~ Cytotoxic_cohort_SD + regimen_group, design_data)
  term <- which(colnames(X) == "Cytotoxic_cohort_SD")
  est <- firth_profile_ci(X, keep$response_binary, term, profile = FALSE)
  data.table(omitted_patient_id = cc$patient_id[i],
             omitted_response = cc$strict_response[i],
             omitted_regimen = as.character(cc$regimen_group[i]),
             omitted_Cytotoxic = cc$Cytotoxic[i],
             n_remaining = nrow(keep),
             OR_per_cohort_SD = est["OR"],
             OR_per_full_C1_SD = exp(est["logOR"] * sd(cc$Cytotoxic)))
}))
fwrite(loo_rows, file.path(res, "Gide_C1_cytotoxic_regimen_adjusted_leave_one_out.tsv"),
       sep = "\t", na = "NA")
lo_idx <- which.min(loo_rows$OR_per_cohort_SD)
hi_idx <- which.max(loo_rows$OR_per_cohort_SD)
loo_summary <- data.table(
  n_leave_one_out_fits = nrow(loo_rows),
  full_data_OR_per_cohort_SD = firth_rows[model == "Regimen-adjusted", OR_per_SD],
  leave_one_out_OR_min = loo_rows$OR_per_cohort_SD[lo_idx],
  omitted_patient_at_min = loo_rows$omitted_patient_id[lo_idx],
  leave_one_out_OR_max = loo_rows$OR_per_cohort_SD[hi_idx],
  omitted_patient_at_max = loo_rows$omitted_patient_id[hi_idx],
  full_C1_SD_in_cohort_SD_units = sd(cc$Cytotoxic),
  interpretation = "point-estimate influence check only; no p values or intervals"
)
fwrite(loo_summary, file.path(res, "Gide_C1_cytotoxic_leave_one_out_summary.tsv"),
       sep = "\t", na = "NA")

# 以患者为单位 bootstrap Spearman 区间；该轴是既有 11 个非 AMP LM12 ssGSEA 分数均值。
set.seed(20260923)
n_boot <- 5000L
rho <- suppressWarnings(cor(x$Cytotoxic, x$LM12_inflammatory_axis, method = "spearman"))
boot_rho <- replicate(n_boot, {
  idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
  suppressWarnings(cor(x$Cytotoxic[idx], x$LM12_inflammatory_axis[idx], method = "spearman"))
})
rho_ci <- quantile(boot_rho[is.finite(boot_rho)], c(.025, .975), names = FALSE)
corr_result <- data.table(
  cohort = "Gide", C_state = "C1", n = nrow(x),
  x_score = "Cytotoxic effector programme score",
  y_score = "LM12 inflammatory axis (mean of 11 non-AMP ssGSEA scores)",
  spearman_rho = rho, CI_low = rho_ci[1], CI_high = rho_ci[2],
  interval = paste0("95% percentile patient bootstrap; B=", n_boot),
  seed = 20260923L, finite_bootstrap = sum(is.finite(boot_rho)),
  p_value = NA_real_, multiplicity_adjustment = "not performed; descriptive association only"
)
fwrite(corr_result, file.path(res, "Gide_C1_cytotoxic_LM12_axis_correlation.tsv"),
       sep = "\t", na = "NA")

# 检查 Cytotoxic 六基因是否直接参与 LM12 轴成员定义，以限定相关性的解释。
program_members <- fread(unname(cfg_paths["frozen_program_members"]))
if (!all(c("program", "gene") %in% names(program_members))) {
  stop("Unexpected frozen programme member table schema")
}
cytogenes <- unique(toupper(program_members[program == "Cytotoxic", gene]))
lm_long <- fread(unname(cfg_paths["lm12_members"]))
axis_genes <- unique(toupper(lm_long[feature_id != "LM12_Antimicrobial_Peptides", gene]))
overlap <- data.table(gene = sort(cytogenes),
                      member_of_LM12_inflammatory_axis = sort(cytogenes) %in% axis_genes,
                      axis_definition = "11 non-AMP LM12 gene sets")
fwrite(overlap, file.path(res, "Gide_C1_cytotoxic_LM12_member_overlap.tsv"),
       sep = "\t", na = "NA")

# 保留本次实际纳入的 25 名患者及用于本探索的冻结分数/临床字段。
patient_subset <- x[, .(patient_id, sample_id, strict_response, regimen_source = regimen,
                         regimen_group, Cytotoxic, LM12_inflammatory_axis, Acid_Score)]
setorder(patient_subset, patient_id)
fwrite(patient_subset, file.path(res, "Gide_C1_patients_used.tsv"), sep = "\t", na = "NA")

# 单张诊断图：散点同时编码疗效与方案，旁图呈现各方案下的 Cytotoxic 分布。
response_cols <- c("SD/PD" = "#5279a6", "CR/PR" = "#b85a50")
regimen_shapes <- c("Monotherapy" = 16, "Combination" = 17)
theme_explore <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9),
        legend.position = "bottom", legend.title = element_text(face = "bold"))
p_axis <- ggplot(x, aes(LM12_inflammatory_axis, Cytotoxic,
                        colour = response_label, shape = regimen_group)) +
  geom_point(size = 2.8, alpha = .9) +
  scale_colour_manual(values = response_cols, drop = FALSE) +
  scale_shape_manual(values = regimen_shapes, drop = FALSE) +
  labs(title = "Gide C1: cytotoxicity and inflammation scores",
       subtitle = sprintf("n=25; Spearman rho=%.2f (95%% bootstrap CI %.2f to %.2f)",
                          rho, rho_ci[1], rho_ci[2]),
       x = "LM12 inflammatory axis (mean of 11 ssGSEA scores)",
       y = "Cytotoxic programme score (cohort SD)",
       colour = "Best response", shape = "Pretreatment regimen") +
  theme_explore

p_regimen <- ggplot(x, aes(regimen_group, Cytotoxic,
                           colour = response_label)) +
  geom_boxplot(aes(group = regimen_group), width = .45,
               outlier.shape = NA, alpha = .15, colour = "#555555") +
  geom_point(position = position_jitter(width = .09, height = 0, seed = 20260923),
             size = 2.5, alpha = .9) +
  scale_colour_manual(values = response_cols, drop = FALSE) +
  scale_x_discrete(labels = c("Monotherapy" = paste0("Anti-PD-1\nmonotherapy\nn=",
                                                       regimen_summary[regimen_group == "Monotherapy", n]),
                              "Combination" = paste0("Anti-PD-1 +\nanti-CTLA-4\nn=",
                                                      regimen_summary[regimen_group == "Combination", n])),
                   drop = FALSE) +
  labs(title = "Cytotoxic scores by regimen",
       subtitle = "Each point is one pretreatment C1 patient",
       x = NULL, y = "Cytotoxic programme score (cohort SD)",
       colour = "Best response") + theme_explore

draw_diagnostic <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.15, 1), "null"))))
  print(p_axis, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_regimen, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  popViewport()
}
png(file.path(fig, "Gide_C1_regimen_axis_diagnostic.png"),
    width = 2600, height = 1350, res = 180)
draw_diagnostic()
dev.off()
pdf(file.path(fig, "Gide_C1_regimen_axis_diagnostic.pdf"), width = 13, height = 6.75)
draw_diagnostic()
dev.off()

expected <- c(
  "Gide_C1_regimen_response_distribution.tsv",
  "Gide_C1_cytotoxic_by_regimen_response.tsv",
  "Gide_C1_cytotoxic_Firth_regimen_sensitivity.tsv",
  "Gide_C1_cytotoxic_regimen_adjusted_leave_one_out.tsv",
  "Gide_C1_cytotoxic_leave_one_out_summary.tsv",
  "Gide_C1_cytotoxic_LM12_axis_correlation.tsv",
  "Gide_C1_cytotoxic_LM12_member_overlap.tsv",
  "Gide_C1_patients_used.tsv")
if (any(!file.exists(file.path(res, expected))) ||
    any(file.info(file.path(res, expected))$size <= 0) ||
    any(file.info(file.path(fig, c("Gide_C1_regimen_axis_diagnostic.png",
                                  "Gide_C1_regimen_axis_diagnostic.pdf")))$size <= 0)) {
  stop("One or more exploratory outputs are missing or empty")
}
log_event("CHECKPOINT", sprintf(
  "Gide C1 n=%d; R=%d NR=%d; regimen levels=%s; cytotoxic-LM12 rho=%.3f",
  nrow(x), sum(x$strict_response == "R"), sum(x$strict_response == "NR"),
  paste(levels(droplevels(x$regimen_group)), collapse = "/"), rho))
log_event("DONE", paste0("outputs written under ", out,
                         "; exploratory sensitivity only; fixed C1 three-program results unchanged"))
}
