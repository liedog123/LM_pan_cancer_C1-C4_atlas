#!/usr/bin/env Rscript
# Purpose: Check code syntax, task wiring and metadata guards without running scientific analyses.
# Inputs: current code checkout; optional output report path outside this checkout.
# Outputs: TSV check report at the explicitly supplied new path, or standard output.
# Usage: Rscript checks/run_checks.R /path/to/new/report.tsv
options(stringsAsFactors=FALSE, warn=1)
args <- commandArgs(TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(self), ".."), mustWork=TRUE)
if (length(args)>1L) stop("Expected at most one report path.")
if(length(args)) {
  report <- normalizePath(dirname(args[1L]),mustWork=TRUE)
  report <- file.path(report,basename(args[1L]))
  if(file.exists(report) || startsWith(report,paste0(root,"/"))) stop("Report must be a new file outside the checkout.")
}
idx <- read.delim(file.path(root,"docs/task_index.tsv"),check.names=FALSE)
si <- read.delim(file.path(root,"docs/step_index.tsv"),check.names=FALSE)
rs <- file.path(root,"internal",unique(si$script))
answers <- list()
record <- function(name,pass,detail="") {
 answers[[length(answers)+1L]] <<- data.frame(check=name,pass=isTRUE(pass),detail=detail)
}
rscript <- file.path(R.home("bin"),"Rscript")
call <- function(script, arg=character()) {
 out <- suppressWarnings(system2(rscript,c(shQuote(script),vapply(arg,shQuote,character(1))),stdout=TRUE,stderr=TRUE))
 code <- attr(out,"status"); if(is.null(code))code<-0L
 list(status=code,out=out)
}
for(p in c(file.path(root,"run.R"),rs,normalizePath(self))) {
 result <- tryCatch({parse(p);TRUE},error=function(e)FALSE)
 record(paste0("parse_",basename(p)),result)
}
record("all_implementations_indexed",setequal(normalizePath(rs),normalizePath(list.files(file.path(root,"internal"),pattern="[.]R$",recursive=TRUE,full.names=TRUE))))
record("tasks_match_steps",setequal(idx$task_id,si$task_id))
record("scope_register_covers_tasks",setequal(idx$task_id,read.delim(file.path(root,"docs/manuscript_scope.tsv"))$task_id))
record("step_counts",all(vapply(idx$task_id,function(id)sum(si$task_id==id),integer(1))==idx$steps))
for(p in rs) {
 h <- readLines(p,n=8L,warn=FALSE)
 record(paste0("header_",basename(p)),all(vapply(c("# Purpose:","# Task:","# Inputs (DCF fields):","# Outputs:","# Usage:"),function(x)any(startsWith(h,x)),logical(1))))
 z <- call(p,"--help"); record(paste0("help_",basename(p)),z$status==0L)
}
z<-call(file.path(root,"run.R"),"--help");record("dispatcher_help",z$status==0L)
z<-call(file.path(root,"run.R"),"--list");record("dispatcher_list",z$status==0L&&length(z$out)==nrow(idx)+1L)
for(id in idx$task_id) {
 z<-call(file.path(root,"run.R"),c("--task",id,"--plan"))
 record(paste0("plan_",id),z$status==0L)
}
for(arg in list(c("--task","bassez_cellchat","--plan"),c("--task","drug_crossplatform","--steps","s02","--plan"),c("--task","ag_survival","--steps","missing","--plan"))) {
 z<-call(file.path(root,"run.R"),arg);record(paste0("reject_",paste(arg,collapse="_")),z$status!=0L)
}
# 伪文件只验证元信息早退，不包含研究行，也不进入统计或绘图代码。
fixture<-tempfile("c1c4_metadata_");dir.create(fixture)
inp<-file.path(fixture,"inputs");cfgdir<-file.path(fixture,"configs");pred<-file.path(fixture,"prediction")
dir.create(inp);dir.create(cfgdir);dir.create(pred)
fake<-function(n){p<-file.path(inp,n);writeLines("metadata fixture; not a study object",p);p}
make_cfg<-function(vals,name) {
 p<-file.path(cfgdir,paste0(name,".dcf"));write.dcf(as.data.frame(as.list(vals),check.names=FALSE),file=p,width=100000L);p
}
prism<-c(gdsc_response=fake("GDSC2.xlsx"),prism_curve=fake("secondary-screen-dose-response-top8-filtered.csv"),prism_cell_metadata=fake("secondary-screen-cell-line-info.csv"),output_dir=file.path(fixture,"prism_output"))
axis<-c(drug_complete=fake("complete.tsv"),shared_patient=fake("shared.tsv"),prediction_audit=fake("audit.tsv"),prediction_dir=pred,output_dir=file.path(fixture,"axis_output"))
surv<-c(projects_data=fake("projects.rdata"),clinical_groups=fake("clinical.rdata"),output_dir=file.path(fixture,"survival_output"))
numeric_fields <- c("patient","mrna","program_genes","program_clusters","rppa","rppa_dict","tcpa","tmb","mutation","cnv","estimate","cibersort","mcpcounter","drug")
numeric <- setNames(vapply(numeric_fields,function(x)fake(paste0(x,".tsv")),character(1)),numeric_fields)
numeric <- c(numeric,output_dir=file.path(fixture,"numeric_output"))
for(item in list(list(id="prism",script="drug_response/03_run_prism_concordance.R",vals=prism),list(id="axis",script="drug_response/04_run_shared_axis_adjustment.R",vals=axis),list(id="survival",script="ag_survival/00_prepare_ag_survival_inputs.R",vals=surv),list(id="numeric",script="drug_response/05_build_programme_numeric_inputs.R",vals=numeric))) {
 config<-make_cfg(item$vals,item$id)
 z<-call(file.path(root,"internal",item$script),c("--config",config,"--check-config"))
 record(paste0(item$id,"_metadata_early_exit"),z$status==0L&&!file.exists(item$vals[["output_dir"]]))
 wrong<-c(item$vals,unexpected="obsolete_dependency")
 config<-make_cfg(wrong,paste0(item$id,"_wrong_fields"))
 z<-call(file.path(root,"internal",item$script),c("--config",config,"--check-config"))
 record(paste0(item$id,"_reject_unknown_field"),z$status!=0L)
 vals<-item$vals;vals[["output_dir"]]<-inp;config<-make_cfg(vals,paste0(item$id,"_overlap"))
 z<-call(file.path(root,"internal",item$script),c("--config",config,"--check-config"))
 record(paste0(item$id,"_reject_output_overlap"),z$status!=0L)
}
# 调度器必须接受新的精简字段，同时保护非空输出。
p<-make_cfg(prism,"dispatcher_prism")
header<-paste0("Task: drug_crossplatform\nLayer: supporting\nOutput-Root: ",file.path(fixture,"dispatcher_logs"),"\n\nStep: s03\nScript: drug_response/03_run_prism_concordance.R\n")
writeLines(paste0(header,paste(readLines(p),collapse="\n")),p)
z<-call(file.path(root,"run.R"),c("--task","drug_crossplatform","--config",p,"--check-config"))
record("dispatcher_prism_metadata_only",z$status==0L&&!dir.exists(prism[["output_dir"]]))
dir.create(prism[["output_dir"]]);sentinel<-file.path(prism[["output_dir"]],"keep.txt");writeLines("keep",sentinel)
z<-call(file.path(root,"run.R"),c("--task","drug_crossplatform","--config",p,"--check-config"))
record("dispatcher_protect_nonempty_output",z$status!=0L&&identical(readLines(sentinel),"keep"))
unlink(fixture,recursive=TRUE)
ans<-do.call(rbind,answers)
if(length(args))write.table(ans,report,sep="\t",quote=TRUE,row.names=FALSE) else write.table(ans,stdout(),sep="\t",quote=TRUE,row.names=FALSE)
cat(sprintf("Checks: %d; passed: %d; failed: %d\n",nrow(ans),sum(ans$pass),sum(!ans$pass)))
if(any(!ans$pass))quit(status=1L)
