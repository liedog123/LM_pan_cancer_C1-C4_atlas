#!/usr/bin/env Rscript
# Purpose: Dispatch manuscript-mapped analysis tasks in isolated R processes.
# Inputs: task DCFs and the study inputs configured for each selected step.
# Outputs: independent step directories and execution logs under Output-Root.
# Usage: Rscript run.R --list; see docs/task_usage.md.
# 任务调度保持原进程隔离，不修改统计范围、随机种子或输入文件。
options(stringsAsFactors = FALSE, warn = 1)
die <- function(...) stop(..., call. = FALSE)
usage <- function() cat(
"Usage: Rscript run.R --list\n",
"       Rscript run.R --task ID [--config PATH] [--steps STEP1,STEP2] [--plan | --check-config]\n",
"       Rscript run.R --help\n\n",
"--plan validates task structure and prints paths without reading study contents or creating outputs.\n",
"--check-config invokes only the retained step metadata checks; downstream outputs must already exist.\n",
"Ordinary invocation runs selected steps in configuration order. Existing matching intermediates may be supplied.\n",
"Input definitions and acquisition boundaries are documented in docs/reproduction_routes.md.\n", sep="")
argv <- commandArgs(trailingOnly=TRUE)
if(identical(argv,"--help")) {usage();quit(status=0L)}
full <- commandArgs(trailingOnly=FALSE)
self <- sub("^--file=","",full[grepl("^--file=",full)])
if(length(self)!=1L) die("Invoke this entry with Rscript.")
root <- dirname(normalizePath(self,mustWork=TRUE))
inside <- function(p,d) identical(p,d) || startsWith(p,paste0(d,"/"))
resolve <- function(p,base) {
 p <- path.expand(p)
 if(!grepl("^/",p)) p <- file.path(base,p)
 if(file.exists(p)||dir.exists(p)) return(normalizePath(p,mustWork=TRUE))
 parent <- dirname(p)
 if(identical(parent,p)) die("Cannot resolve path: ",p)
 pp <- resolve(parent,base);leaf <- basename(p)
 if(leaf==".")return(pp)
 if(leaf=="..")return(dirname(pp))
 file.path(pp,leaf)
}
read_records <- function(p) {
 if(!file_test("-f",p)||file.access(p,4L)!=0L||is.na(file.info(p)$size)||file.info(p)$size<=0)die("Readable nonempty task config required.")
 ln <- readLines(p,warn=FALSE)
 ln <- ln[!grepl("^\\s*#",ln)]
 con <- textConnection(ln);on.exit(close(con))
 x <- read.dcf(con)
 lapply(seq_len(nrow(x)),function(i) {
  r <- setNames(as.character(x[i,!is.na(x[i,]),drop=TRUE]),colnames(x)[!is.na(x[i,])])
  r[] <- gsub("[\\r\\n]+[ \\t]*", " ", r, perl=TRUE)
  r
 })
}
index <- read.delim(file.path(root,"docs","task_index.tsv"),check.names=FALSE,quote="",na.strings="",colClasses="character")
steps_index <- read.delim(file.path(root,"docs","step_index.tsv"),check.names=FALSE,quote="",na.strings="",colClasses="character")
if(identical(argv,"--list")) {write.table(index[,c("task_id","layer","title")],stdout(),sep="\t",quote=FALSE,row.names=FALSE);quit(status=0L)}
opt <- list(task=NULL,config=NULL,steps=NULL,mode="run")
i <- 1L;seen <- character()
while(i<=length(argv)) {
 a<-argv[i]
 if(a %in% c("--plan","--check-config")) {
  if(opt$mode!="run")die("Choose one mode.")
  opt$mode<-sub("^--","",a);i<-i+1L;next
 }
 if(!a %in% c("--task","--config","--steps")||i==length(argv)||a %in% seen)die("Invalid or repeated argument: ",a)
 k<-sub("^--","",a);opt[[k]]<-argv[i+1L];seen<-c(seen,a);i<-i+2L
}
if(is.null(opt$task)||!opt$task %in% index$task_id)die("Unknown task. Use --list.")
row <- index[index$task_id==opt$task,,drop=FALSE]
config <- if(is.null(opt$config)) file.path(root,row$config) else opt$config
config <- normalizePath(config,mustWork=TRUE);base <- dirname(config)
records <- read_records(config)
if(length(records)<2L)die("Task config needs a header and at least one step.")
head <- records[[1L]];st <- records[-1L]
required_head <- c("Task","Layer","Output-Root")
if(!all(required_head %in% names(head))||!identical(head["Task"][[1L]],opt$task)||!identical(head["Layer"][[1L]],row$layer))die("Task header differs from task index.")
if(any(!names(head) %in% c(required_head,names(head)[startsWith(names(head),"Input-")])))die("Unknown header field.")
keys <- vapply(st,function(x)if("Step" %in% names(x))x[["Step"]] else "",character(1))
expected <- steps_index[steps_index$task_id==opt$task,,drop=FALSE]
if(anyDuplicated(keys)||!identical(keys,expected$step_id))die("Step identifiers/order differ from retained task map.")
names(st)<-keys
for(j in seq_along(st)) {
 r<-st[[j]];schema<-strsplit(expected$fields[j],",",fixed=TRUE)[[1L]]
 if(!setequal(names(r),c("Step","Script",schema))||length(r)!=length(schema)+2L)die("Invalid fields in step ",keys[j])
 if(r[["Script"]]!=expected$script[j])die("Implementation differs from retained task map.")
 p<-normalizePath(file.path(root,"internal",r[["Script"]]),mustWork=TRUE)
 if(!inside(p,normalizePath(file.path(root,"internal"),mustWork=TRUE))||!file_test("-f",p))die("Implementation escapes internal directory.")
}
# 引用仅替换路径；不修改任何研究参数，也不自动替换冻结标签。
expand <- function(v,stack=character()) {
 if(!startsWith(v,"@"))return(v)
 if(v %in% stack)die("Cyclic config reference: ",v)
 m<-regmatches(v,regexec("^@([^/]+)/?(.*)$",v))[[1L]]
 if(length(m)!=3L)die("Malformed reference: ",v)
 name<-m[2L];tail<-m[3L]
 if(name=="input") {
  k<-paste0("Input-",tail)
  if(!k %in% names(head))die("Unknown shared input: ",tail)
  return(expand(head[[k]],c(stack,v)))
 }
 if(grepl("(^|/)\\.\\.(/|$)",tail)||startsWith(tail,"/"))die("Reference suffix escapes its producer.")
 if(name=="root") prefix<-head[["Output-Root"]] else {
  if(!name %in% keys)die("Unknown producer step: ",name)
  prefix<-expand(st[[name]][["output_dir"]],c(stack,v))
 }
 if(nzchar(tail))file.path(prefix,tail) else prefix
}
fields <- lapply(st,function(r) r[!names(r) %in% c("Step","Script")])
resolved <- lapply(fields,function(r) {
 vals<-vapply(r,expand,character(1))
 path_keys<-setdiff(names(vals),"edge_id")
 vals[path_keys]<-vapply(vals[path_keys],resolve,character(1),base=base)
 vals
})
for(j in seq_along(st)) {
 vals<-st[[j]][!names(st[[j]]) %in% c("Step","Script","output_dir")]
 refs<-vals[grepl("^@s",vals)]
 parents<-unique(sub("^@([^/]+).*$","\\1",refs))
 if(any(match(parents,keys)>=j)||anyNA(match(parents,keys)))die("Producer must precede consumer: ",keys[j])
}
selected <- if(is.null(opt$steps))keys else strsplit(opt$steps,",",fixed=TRUE)[[1L]]
if(!length(selected)||anyDuplicated(selected)||any(!selected %in% keys))die("Unknown or repeated --steps selection.")
selected <- keys[keys %in% selected]
plan <- data.frame(step=selected,script=vapply(st[selected],function(r)r[["Script"]],character(1)),
 output_dir=vapply(resolved[selected],function(r)r[["output_dir"]],character(1)),
 unresolved_fields=vapply(resolved[selected],function(r)paste(names(r)[grepl("SET_ME",r,fixed=TRUE)],collapse=","),character(1)))
write.table(plan,stdout(),sep="\t",quote=FALSE,row.names=FALSE)
if(opt$mode=="plan")quit(status=0L)
if(any(nzchar(plan$unresolved_fields)))die("Replace SET_ME placeholders in the selected steps.")
out <- vapply(resolved[selected],function(r)r[["output_dir"]],character(1))
for(p in out) {
 if(p=="/"||p==normalizePath(path.expand("~"))||inside(p,root)||inside(root,p))die("Output must be independent of the package.")
 if(file.exists(p)&&!dir.exists(p))die("Output is a file: ",p)
 if(dir.exists(p)&&length(list.files(p,all.files=TRUE,no..=TRUE)))die("Output is nonempty: ",p)
}
if(length(out)>1L)for(a in seq_along(out))for(b in seq_along(out))if(a!=b&&(inside(out[a],out[b])||inside(out[b],out[a])))die("Step outputs must be disjoint.")
# 将路径写成绝对路径后交给旧接口，避免临时配置位置改变路径含义。
rscript<-file.path(R.home("bin"),"Rscript")
run_child<-function(k,check) {
 tmp<-tempfile("c1c4_step_",fileext=".dcf");on.exit(unlink(tmp))
 write.dcf(as.data.frame(as.list(resolved[[k]]),check.names=FALSE),file=tmp,width=100000L)
 script<-file.path(root,"internal",st[[k]][["Script"]])
 args<-c(shQuote(script),"--config",shQuote(tmp),if(check)"--check-config")
 status<-system2(rscript,args)
 if(!identical(as.integer(status),0L))die("Step failed: ",k," (exit ",status,"). No later step was run.")
}
if(opt$mode=="check-config") {
 for(k in selected)run_child(k,TRUE)
 cat("TASK_METADATA_CHECK_PASS: no analysis outputs created.\n");quit(status=0L)
}
log_root<-resolve(file.path(head[["Output-Root"]],".task_logs",opt$task),base)
if(inside(log_root,root)||inside(root,log_root))die("Log directory overlaps package.")
dir.create(log_root,recursive=TRUE,showWarnings=FALSE)
log_file<-file.path(log_root,paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid(),".tsv"))
if(file.exists(log_file))die("Run log already exists.")
writeLines("step\tevent\ttime",log_file)
for(k in selected) {
 dir.create(dirname(resolved[[k]][["output_dir"]]),recursive=TRUE,showWarnings=FALSE)
 cat(paste(k,"STARTED",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),sep="\t"),"\n",file=log_file,append=TRUE)
 run_child(k,TRUE)
 run_child(k,FALSE)
 cat(paste(k,"COMPLETED",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),sep="\t"),"\n",file=log_file,append=TRUE)
}
cat("TASK_COMPLETED: ",opt$task,"\n",sep="")
