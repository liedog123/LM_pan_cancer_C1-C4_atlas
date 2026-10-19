#!/usr/bin/env Rscript
# Purpose: Validate syntax, contracts, metadata guards and dispatcher isolation.
# Inputs: code checkout only. Metadata files and synthetic stages contain no study data.
# Outputs: stdout or a new TSV report outside the checkout. No scientific analyses run.
# 使用伪文件和合成阶段检查接口，不读取研究内容或调用研究模型。
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(TRUE)
self<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
root<-normalizePath(file.path(dirname(self),".."),mustWork=TRUE)
if(length(args)>1L)stop("Expected at most one report path.")
for(n in c("config","io","statistics","plotting","runtime"))sys.source(file.path(root,"R/common",paste0(n,".R")),envir=.GlobalEnv)
if(length(args)){
 report<-resolve_path(args[1],getwd())
 if(file.exists(report)||is_within(report,root)||!dir.exists(dirname(report)))die("Report must be a new file outside the checkout.")
}
answers<-list();record<-function(name,ok,detail="")answers[[length(answers)+1L]]<<-data.frame(check=name,pass=isTRUE(ok),detail=detail)
rscript<-file.path(R.home("bin"),"Rscript")
call<-function(script,arg=character()){
 out<-suppressWarnings(system2(rscript,vapply(c(script,arg),shQuote,character(1)),stdout=TRUE,stderr=TRUE))
 status<-attr(out,"status");if(is.null(status))status<-0L
 list(status=status,out=out)
}
for(p in list.files(root,pattern="[.]R$",full.names=TRUE,recursive=TRUE))record(paste0("parse_",sub(paste0(root,"/"),"",p,fixed=TRUE)),tryCatch({parse(p);TRUE},error=function(e)FALSE))
reg<-read.delim(file.path(root,"config/stages.tsv"),colClasses="character",check.names=FALSE)
record("18_analysis_modules",length(list.files(file.path(root,"R/analysis"),pattern="[.]R$"))==18L)
record("5_common_files",length(list.files(file.path(root,"R/common"),pattern="[.]R$"))==5L)
record("8_configuration_templates",length(list.files(file.path(root,"config"),pattern="[.]dcf$"))==8L)
record("87_unique_stage_records",nrow(reg)==87L&&!anyDuplicated(reg$stage))
record("86_distinct_implementations",length(unique(reg$handler))==86L)
record("35_task_aliases",length(unique(reg$task_id))==35L)
for(a in unique(reg$analysis)){
 env<-new.env(parent=globalenv());load_stage_code(root,a,env)
 handlers<-unique(reg$handler[reg$analysis==a])
 record(paste0("handlers_",a),all(vapply(handlers,function(h)exists(h,envir=env,inherits=FALSE)&&is.function(get(h,envir=env)),logical(1))))
 z<-call(file.path(root,"run.R"),c("--module",a,"--plan"))
 record(paste0("analysis_plan_",a),z$status==0L&&any(grepl(handlers[1],z$out,fixed=TRUE)|grepl(a,z$out,fixed=TRUE)))
}
for(d in unique(reg$module)){
 r<-read_records(file.path(root,"config",paste0(d,".dcf")))
 expected<-reg[reg$module==d,,drop=FALSE]
 record(paste0("stage_order_",d),identical(vapply(r[-1],function(x)x[["Stage"]],character(1)),expected$stage))
 for(j in seq_len(nrow(expected))){
  fields<-strsplit(expected$fields[j],",",fixed=TRUE)[[1]]
  record(paste0("schema_",expected$stage[j]),setequal(names(r[[j+1L]]),c("Stage",fields)))
 }
 z<-call(file.path(root,"run.R"),c("--module",d,"--plan"));record(paste0("domain_plan_",d),z$status==0L)
}
for(t in unique(reg$task_id)){
 z<-call(file.path(root,"run.R"),c("--task",t,"--plan"));record(paste0("task_alias_",t),z$status==0L&&all(vapply(reg$stage[reg$task_id==t],function(k)any(grepl(k,z$out,fixed=TRUE)),logical(1))))
}
fixture<-tempfile("c1c4_contract_checks_");dir.create(fixture)
inp<-file.path(fixture,"inputs");cfgdir<-file.path(fixture,"configs");dir.create(inp);dir.create(cfgdir)
fake<-function(n){p<-file.path(inp,n);writeLines("metadata fixture; not a study object",p);p}
make_cfg<-function(vals,name){p<-file.path(cfgdir,paste0(name,".dcf"));write.dcf(as.data.frame(as.list(vals),check.names=FALSE),p,width=100000L);p}
pred<-file.path(inp,"prediction");dir.create(pred)
prism<-c(gdsc_response=fake("GDSC2.xlsx"),prism_curve=fake("secondary-screen-dose-response-top8-filtered.csv"),prism_cell_metadata=fake("secondary-screen-cell-line-info.csv"),output_dir=file.path(fixture,"prism_output"))
axis<-c(drug_complete=fake("complete.tsv"),shared_patient=fake("shared.tsv"),prediction_audit=fake("audit.tsv"),prediction_dir=pred,output_dir=file.path(fixture,"axis_output"))
surv<-c(projects_data=fake("projects.rdata"),clinical_groups=fake("clinical.rdata"),output_dir=file.path(fixture,"survival_output"))
ff<-c("patient","mrna","program_genes","program_clusters","rppa","rppa_dict","tcpa","tmb","mutation","cnv","estimate","cibersort","mcpcounter","drug")
numeric<-c(setNames(vapply(ff,function(n)fake(paste0(n,".tsv")),character(1)),ff),output_dir=file.path(fixture,"numeric_output"))
items<-list(list(script="drug_response/03_run_prism_concordance.R",vals=prism),list(script="drug_response/04_run_shared_axis_adjustment.R",vals=axis),list(script="ag_survival/00_prepare_ag_survival_inputs.R",vals=surv),list(script="drug_response/05_build_programme_numeric_inputs.R",vals=numeric))
for(item in items){
 row<-reg[match(item$script,reg$legacy_script),,drop=FALSE]
 invoke<-function(vals,label){p<-make_cfg(vals,paste0(row$task_id,"_",label));call(file.path(root,"run.R"),c("--worker",row$stage,row$handler,p,"check"))}
 z<-invoke(item$vals,"valid");record(paste0("metadata_only_",row$stage),z$status==0L&&!file.exists(item$vals[["output_dir"]]),paste(z$out,collapse=" "))
 z<-invoke(c(item$vals,unexpected="old_field"),"unknown");record(paste0("unknown_field_",row$stage),z$status!=0L)
 vals<-item$vals;vals[["output_dir"]]<-file.path(root,"forbidden_output");z<-invoke(vals,"package");record(paste0("package_overlap_",row$stage),z$status!=0L)
 vals<-item$vals;vals[["output_dir"]]<-inp;z<-invoke(vals,"overlap");record(paste0("input_overlap_",row$stage),z$status!=0L)
 vals<-item$vals;vals[[names(vals)[1]]]<-file.path(inp,"missing");z<-invoke(vals,"missing");record(paste0("missing_input_",row$stage),z$status!=0L)
}
# Check module dispatch to the same metadata-only worker with all other template paths unresolved.
row<-reg[match(items[[1]]$script,reg$legacy_script),,drop=FALSE]
r<-read_records(file.path(root,"config",paste0(row$module,".dcf")))
r[[1]][["Output-Root"]]<-file.path(fixture,"logs")
k<-which(vapply(r[-1],function(z)z[["Stage"]],character(1))==row$stage)+1L
r[[k]]<-c(Stage=row$stage,prism)
p<-file.path(cfgdir,"dispatch.dcf");writeLines(unlist(lapply(r,function(v)c(paste0(names(v),": ",v),""))),p)
z<-call(file.path(root,"run.R"),c("--task",row$task_id,"--config",p,"--check-config"))
record("dispatcher_metadata_only",z$status==0L&&!dir.exists(prism[["output_dir"]]))
dir.create(prism[["output_dir"]]);sentinel<-file.path(prism[["output_dir"]],"keep.txt");writeLines("keep",sentinel)
z<-call(file.path(root,"run.R"),c("--task",row$task_id,"--config",p,"--check-config"))
record("protect_nonempty_output",z$status!=0L&&identical(readLines(sentinel),"keep"))
# Duplicate config keys must not silently overwrite earlier values.
pdup<-file.path(cfgdir,"duplicate.dcf");writeLines(c("Module: drug","Module: clinical","Output-Root: x"),pdup)
record("reject_duplicate_fields",inherits(try(read_records(pdup),silent=TRUE),"try-error"))
record("path_root_and_parent",is_within("/tmp","/")&&identical(resolve_path("../x",inp),file.path(fixture,"x")))
dangling<-file.path(fixture,"dangling")
file.symlink(file.path(fixture,"absent_target"),dangling)
record("reject_dangling_symlink",inherits(try(resolve_path(dangling,fixture),silent=TRUE),"try-error"))
# Independent child processes and failure-stop behavior are tested with synthetic functions only.
mini<-file.path(fixture,"mini");dir.create(mini)
for(d in c("R/common","R/analysis","config"))dir.create(file.path(mini,d),recursive=TRUE,showWarnings=FALSE)
file.copy(file.path(root,"run.R"),mini)
file.copy(list.files(file.path(root,"R/common"),full.names=TRUE),file.path(mini,"R/common"))
syn<-data.frame(module="probe",analysis="probe",stage=c("probe.s01","probe.s02"),task_id="probe",layer="synthetic",step_id=c("s01","s02"),handler=c("first","second"),fields="output_dir",legacy_script="synthetic")
write.table(syn,file.path(mini,"config/stages.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c('first <- function(cfg,context,check_only=FALSE) { if(check_only)return(invisible(TRUE)); options(c1c4_probe=TRUE); dir.create(cfg[["output_dir"]]); writeLines(as.character(Sys.getpid()),file.path(cfg[["output_dir"]],"pid")) }','second <- function(cfg,context,check_only=FALSE) { if(check_only)return(invisible(TRUE)); stopifnot(is.null(getOption("c1c4_probe"))); dir.create(cfg[["output_dir"]]); writeLines(as.character(Sys.getpid()),file.path(cfg[["output_dir"]],"pid")) }'),file.path(mini,"R/analysis/probe.R"))
make_probe<-function(label){outs<-file.path(fixture,paste0(label,c("_one","_two")));p<-file.path(cfgdir,paste0(label,".dcf"));writeLines(c("Module: probe",paste0("Output-Root: ",file.path(fixture,paste0(label,"_logs"))),"","Stage: probe.s01",paste0("output_dir: ",outs[1]),"","Stage: probe.s02",paste0("output_dir: ",outs[2])),p);list(config=p,outs=outs)}
probe<-make_probe("success");z<-call(file.path(mini,"run.R"),c("--module","probe","--config",probe$config))
record("independent_process_state",z$status==0L&&all(file.exists(file.path(probe$outs,"pid")))&&length(unique(vapply(file.path(probe$outs,"pid"),readLines,character(1))))==2L)
writeLines(c('first <- function(cfg,context,check_only=FALSE) { if(check_only)return(invisible(TRUE)); stop("synthetic failure") }','second <- function(cfg,context,check_only=FALSE) { dir.create(cfg[["output_dir"]]) }'),file.path(mini,"R/analysis/probe.R"))
probe<-make_probe("failure");z<-call(file.path(mini,"run.R"),c("--module","probe","--config",probe$config))
record("stop_after_failed_stage",z$status!=0L&&!any(dir.exists(probe$outs))&&any(grepl("No later stage was run",z$out,fixed=TRUE)))
unlink(fixture,recursive=TRUE)
ans<-do.call(rbind,answers)
if(length(args))write.table(ans,report,sep="\t",quote=TRUE,row.names=FALSE) else write.table(ans,stdout(),sep="\t",quote=TRUE,row.names=FALSE)
cat(sprintf("Checks: %d; passed: %d; failed: %d\n",nrow(ans),sum(ans$pass),sum(!ans$pass)))
if(any(!ans$pass)) {print(ans[!ans$pass,,drop=FALSE]);quit(status=1L)}
