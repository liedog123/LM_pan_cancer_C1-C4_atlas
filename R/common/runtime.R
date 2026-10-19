# Purpose: Shared runtime helpers; exact matching implementations only.
# Helpers are rebound to the calling stage scope to retain original lexical lookup.
# 仅合并实现完全相同的函数；作用域保留在各独立阶段。

shared_log_event_008 <- function(level, msg, ...) {
    extras <- list(...)
    extra_text <- if (length(extras) == 0L)
        ""
    else paste(vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)), collapse = " ")
    cat(sprintf("[%s] %s script=%s msg=%s%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg), if (nzchar(extra_text))
        " "
    else "", extra_text))
    flush.console()
}

shared_log_event_009 <- function(level, msg, ...) {
    extras <- list(...)
    extra_text <- if (length(extras) == 0L)
        ""
    else paste(vapply(names(extras), function(nm) sprintf("%s=%s", nm, shQuote(as.character(extras[[nm]]))), character(1)), collapse = " ")
    cat(sprintf("[%s] %s script=%s msg=%s%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, shQuote(msg), ifelse(nzchar(extra_text), " ", ""), extra_text))
    flush.console()
}

shared_log_event_010 <- function(level, msg, fields = list()) {
    clean <- function(x) gsub("\"", "'", as.character(x), fixed = TRUE)
    extra <- if (!length(fields))
        ""
    else paste0(" ", paste(sprintf("%s=\"%s\"", names(fields), vapply(fields, clean, character(1))), collapse = " "))
    cat(sprintf("[%s] %s script=\"%s\" msg=\"%s\"%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, clean(msg), extra))
    flush.console()
}

shared_log_event_011 <- function(level, msg, fields = list()) {
    extra <- if (length(fields))
        paste(vapply(names(fields), function(nm) {
            sprintf("%s=\"%s\"", nm, gsub("[\r\n\t]", " ", as.character(fields[[nm]])))
        }, character(1)), collapse = " ")
    else ""
    cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, gsub("[\r\n]", " ", msg), if (nzchar(extra))
        " "
    else "", extra))
    flush.console()
}

shared_log_event_012 <- function(level, msg, fields = list()) {
    extra <- if (length(fields)) {
        paste(vapply(names(fields), function(nm) {
            sprintf("%s=\"%s\"", nm, gsub("[\r\n\t]", " ", as.character(fields[[nm]])))
        }, character(1)), collapse = " ")
    }
    else ""
    cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, gsub("[\r\n]", " ", msg), if (nzchar(extra))
        " "
    else "", extra))
    flush.console()
}

shared_log_event_013 <- function(level, msg, fields = list()) {
    extra <- if (length(fields)) {
        paste(vapply(names(fields), function(nm) {
            sprintf("%s=\"%s\"", nm, gsub("[\r\n\t]", " ", as.character(fields[[nm]])))
        }, character(1)), collapse = " ")
    }
    else {
        ""
    }
    cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, gsub("[\r\n]", " ", msg), if (nzchar(extra))
        " "
    else "", extra))
    flush.console()
}

shared_log_event_014 <- function(level, msg, fields = list()) {
    extra <- if (length(fields)) {
        paste(vapply(names(fields), function(nm) {
            value <- gsub("[\r\n\t]", " ", as.character(fields[[nm]]))
            sprintf("%s=\"%s\"", nm, value)
        }, character(1)), collapse = " ")
    }
    else {
        ""
    }
    cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, gsub("[\r\n]", " ", msg), if (nzchar(extra))
        " "
    else "", extra))
    flush.console()
}

shared_log_event_015 <- function(level, msg, fields = list()) {
    extras <- if (length(fields))
        paste(vapply(names(fields), function(k) {
            sprintf("%s=\"%s\"", k, gsub("\"", "'", as.character(fields[[k]]), fixed = TRUE))
        }, character(1)), collapse = " ")
    else ""
    cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, gsub("\"", "'", msg, fixed = TRUE), if (nzchar(extras))
        " "
    else "", extras))
    flush.console()
}

shared_log_msg_016 <- function(level, step, ...) {
    detail <- paste(..., collapse = " ")
    cat(sprintf("%s\t%s\t%s\t%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, step, detail))
    flush.console()
}

shared_fail_018 <- function(msg) {
    log_event("ERROR", msg)
    stop(msg, call. = FALSE)
}

shared_save_session_info_019 <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(capture.output(sessionInfo()), path)
}

# Purpose: Dispatch a selected stage in a fresh process; fail before later stages on error.
# 仅调度器负责日志和进程；加载模块本身不加载研究输入或执行计算。
load_stage_code <- function(root,analysis,scope) {
 path <- file.path(root,"R/analysis",paste0(analysis,".R"))
 if (!file_test("-f",path)) die("Missing analysis module: ",analysis)
 sys.source(path,envir=scope,keep.source=FALSE)
}
worker_main <- function(args,root) {
 if (length(args)!=4L || !args[4] %in% c("check","run")) die("Invalid worker invocation.")
 reg <- read.delim(file.path(root,"config/stages.tsv"),check.names=FALSE,colClasses="character")
 row <- reg[reg$stage==args[1],,drop=FALSE]
 if (nrow(row)!=1L || row$handler!=args[2]) die("Unknown worker stage or handler.")
 config <- normalizePath(args[3],mustWork=TRUE)
 cfg <- read_stage_config(config,strsplit(row$fields,",",fixed=TRUE)[[1]])
 output <- resolve_path(cfg[["output_dir"]],dirname(config))
 if (output=="/" || output==normalizePath(path.expand("~")) || is_within(output,root) || is_within(root,output)) die("Output must be independent of the package.")
 scope <- new.env(parent=globalenv())
 load_stage_code(root,row$analysis,scope)
 get(row$handler,envir=scope,inherits=FALSE)(cfg,list(root=root,config=config),args[4]=="check")
 invisible(NULL)
}
usage <- function() cat(
 "Usage: Rscript run.R --list\n",
 "       Rscript run.R --module MODULE [--stage TASK_OR_STAGE,...] [--config PATH] [--plan | --check-config]\n",
 "       Rscript run.R --task LEGACY_TASK [--steps STEP,...] [--config PATH] [--plan | --check-config]\n",
 "MODULE accepts one of the 8 configuration domains or 18 analysis names.\n",
 "Copy the corresponding template outside the checkout and replace SET_ME before executing.\n",
 "--plan reads configuration only; --check-config checks filesystem metadata only.\n",
 "See docs/reproduction.md for input boundaries and prepared-input handoffs.\n",sep="")
cli_main <- function(argv,root) {
 if (identical(argv,"--help")) {usage();return(invisible(NULL))}
 if (length(argv) && argv[1]=="--worker") return(worker_main(argv[-1],root))
 reg <- read.delim(file.path(root,"config/stages.tsv"),check.names=FALSE,colClasses="character")
 if (identical(argv,"--list")) {
  write.table(unique(reg[,c("module","analysis","task_id","layer")]),stdout(),sep="\t",quote=FALSE,row.names=FALSE)
  return(invisible(NULL))
 }
 opt <- list(module=NULL,task=NULL,stage=NULL,steps=NULL,config=NULL,mode="run")
 i <- 1L; seen <- character()
 while(i<=length(argv)) {
  a <- argv[i]
  if (a %in% c("--plan","--check-config")) {
   if (opt$mode!="run") die("Choose one mode.")
   opt$mode <- sub("^--","",a);i<-i+1L;next
  }
  if (!a %in% c("--module","--task","--stage","--steps","--config") || i==length(argv) || a %in% seen) die("Invalid or repeated argument: ",a)
  opt[[sub("^--","",a)]] <- argv[i+1L];seen<-c(seen,a);i<-i+2L
 }
 if (is.null(opt$module)==is.null(opt$task)) die("Choose exactly one --module or --task.")
 if (!is.null(opt$task)) {
  if (!opt$task %in% reg$task_id || !is.null(opt$stage)) die("Unknown task or incompatible --stage.")
  route <- reg[reg$task_id==opt$task,,drop=FALSE]
 } else {
  if (!is.null(opt$steps)) die("--steps requires --task.")
  if (opt$module %in% reg$module) route <- reg[reg$module==opt$module,,drop=FALSE] else if(opt$module %in% reg$analysis) route <- reg[reg$analysis==opt$module,,drop=FALSE] else die("Unknown module.")
 }
 domain <- unique(route$module);stopifnot(length(domain)==1L)
 expected <- reg[reg$module==domain,,drop=FALSE]
 config <- if(is.null(opt$config))file.path(root,"config",paste0(domain,".dcf")) else opt$config
 config <- normalizePath(config,mustWork=TRUE);base <- dirname(config)
 records <- read_records(config)
 if (length(records)<2L) die("Module config needs a header and stages.")
 head <- records[[1]];st <- records[-1]
 if (!all(c("Module","Output-Root") %in% names(head)) || head[["Module"]]!=domain || any(!names(head) %in% c("Module","Output-Root",names(head)[startsWith(names(head),"Input-")]))) die("Module header differs from the contract.")
 keys <- vapply(st,function(r)if("Stage" %in% names(r))r[["Stage"]] else "",character(1))
 if (anyDuplicated(keys) || !identical(keys,expected$stage)) die("Stage identifiers/order differ from config/stages.tsv.")
 names(st)<-keys
 for(j in seq_along(st)) {
  schema <- strsplit(expected$fields[j],",",fixed=TRUE)[[1]]
  if (!setequal(names(st[[j]]),c("Stage",schema)) || length(st[[j]])!=length(schema)+1L) die("Invalid fields in stage: ",keys[j])
 }
 expand <- function(v,stack=character()) {
  if(!startsWith(v,"@"))return(v)
  if(v %in% stack)die("Cyclic config reference: ",v)
  m <- regmatches(v,regexec("^@([^/]+)/?(.*)$",v))[[1]]
  if(length(m)!=3L)die("Malformed reference: ",v)
  name<-m[2];tail<-m[3]
  if(name=="input") {
   k<-paste0("Input-",tail)
   if(!k %in% names(head))die("Unknown shared input: ",tail)
   return(expand(head[[k]],c(stack,v)))
  }
  if(grepl("(^|/)\\.\\.(/|$)",tail) || startsWith(tail,"/"))die("Reference suffix escapes its producer.")
  if(name=="root")prefix<-expand(head[["Output-Root"]],c(stack,v)) else {
   if(!name %in% keys)die("Unknown producer: ",name)
   prefix<-expand(st[[name]][["output_dir"]],c(stack,v))
  }
  if(nzchar(tail))file.path(prefix,tail) else prefix
 }
 resolved <- lapply(st,function(r) {
  vals <- vapply(r[names(r)!="Stage"],expand,character(1))
  path_keys<-setdiff(names(vals),"edge_id")
  vals[path_keys]<-vapply(vals[path_keys],resolve_path,character(1),base_dir=base)
  vals
 })
 for(j in seq_along(st)) {
  vals<-st[[j]][!names(st[[j]]) %in% c("Stage","output_dir")]
  refs<-vals[startsWith(vals,"@") & !grepl("^@(input|root)/",vals)]
  parents<-unique(sub("^@([^/]+).*$","\\1",refs))
  if(anyNA(match(parents,keys)) || any(match(parents,keys)>=j))die("Producer must precede consumer: ",keys[j])
 }
 selected <- route$stage
 if(!is.null(opt$stage)) {
  tokens<-strsplit(opt$stage,",",fixed=TRUE)[[1]]
  if(!length(tokens)||anyDuplicated(tokens)||any(!tokens %in% c(route$stage,route$task_id)))die("Unknown or repeated --stage selection.")
  selected<-route$stage[route$stage %in% tokens | route$task_id %in% tokens]
 }
 if(!is.null(opt$steps)) {
  tokens<-strsplit(opt$steps,",",fixed=TRUE)[[1]]
  if(!length(tokens)||anyDuplicated(tokens)||any(!tokens %in% route$step_id))die("Unknown or repeated --steps selection.")
  selected<-route$stage[route$step_id %in% tokens]
 }
 if(!length(selected))die("No stages selected.")
 plan <- data.frame(stage=selected,analysis=reg$analysis[match(selected,reg$stage)],
  output_dir=vapply(resolved[selected],function(r)r[["output_dir"]],character(1)),
  unresolved_fields=vapply(resolved[selected],function(r)paste(names(r)[grepl("SET_ME",r,fixed=TRUE)],collapse=","),character(1)))
 write.table(plan,stdout(),sep="\t",quote=FALSE,row.names=FALSE)
 if(opt$mode=="plan")return(invisible(NULL))
 if(any(nzchar(plan$unresolved_fields)))die("Replace SET_ME placeholders in the selected stages.")
 outputs<-vapply(resolved[selected],function(r)r[["output_dir"]],character(1))
 for(p in outputs) {
  if(p=="/" || p==normalizePath(path.expand("~")) || is_within(p,root)||is_within(root,p))die("Output must be independent of the package.")
  if(file.exists(p)&&!dir.exists(p))die("Output is a file: ",p)
  if(dir.exists(p)&&length(list.files(p,all.files=TRUE,no..=TRUE)))die("Output is nonempty: ",p)
 }
 if(length(outputs)>1L)for(a in seq_along(outputs))for(b in seq_along(outputs))if(a!=b&&(is_within(outputs[a],outputs[b])||is_within(outputs[b],outputs[a])))die("Stage outputs must be disjoint.")
 run_child <- function(k,check) {
  tmp<-tempfile("c1c4_stage_",fileext=".dcf");on.exit(unlink(tmp))
  write.dcf(as.data.frame(as.list(resolved[[k]]),check.names=FALSE),file=tmp,width=100000L)
  row<-reg[match(k,reg$stage),,drop=FALSE]
  args<-c(file.path(root,"run.R"),"--worker",k,row$handler,tmp,if(check)"check" else "run")
  status<-system2(file.path(R.home("bin"),"Rscript"),vapply(args,shQuote,character(1)))
  if(!identical(as.integer(status),0L))die("Stage failed: ",k," (exit ",status,"). No later stage was run.")
 }
 if(opt$mode=="check-config") {
  for(k in selected)run_child(k,TRUE)
  cat("MODULE_METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
  return(invisible(NULL))
 }
 log_root<-resolve_path(file.path(expand(head[["Output-Root"]]),".task_logs",if(is.null(opt$task))opt$module else opt$task),base)
 if(is_within(log_root,root)||is_within(root,log_root))die("Log directory overlaps package.")
 dir.create(log_root,recursive=TRUE,showWarnings=FALSE)
 log_file<-file.path(log_root,paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid(),".tsv"))
 if(file.exists(log_file))die("Run log already exists.")
 writeLines("stage\tevent\ttime",log_file)
 for(k in selected) {
  dir.create(dirname(resolved[[k]][["output_dir"]]),recursive=TRUE,showWarnings=FALSE)
  cat(paste(k,"STARTED",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),sep="\t"),"\n",file=log_file,append=TRUE)
  run_child(k,TRUE);run_child(k,FALSE)
  cat(paste(k,"COMPLETED",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),sep="\t"),"\n",file=log_file,append=TRUE)
 }
 cat("MODULE_COMPLETED\n")
 invisible(NULL)
}
