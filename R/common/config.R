# Purpose: Strict configuration parsing, path resolution and stage bindings.
# Inputs: one module DCF or one internal worker DCF; no study contents are loaded.
# 配置只处理路径与阶段选择；模型、阈值及检验范围保留在分析模块。
die <- function(...) stop(..., call. = FALSE)
is_within <- function(path, root) {
 if (root == "/") return(startsWith(path, "/"))
 identical(path, root) || startsWith(path, paste0(root, "/"))
}
resolve_path <- function(path, base_dir) {
 path <- path.expand(path)
 if (!grepl("^/", path)) path <- file.path(base_dir, path)
 link <- Sys.readlink(path)
 if (!is.na(link) && nzchar(link) && !file.exists(path) && !dir.exists(path)) die("Dangling symbolic link: ",path)
 if (file.exists(path) || dir.exists(path)) return(normalizePath(path, winslash="/", mustWork=TRUE))
 parent <- dirname(path)
 if (identical(parent,path)) die("Cannot resolve path: ",path)
 pp <- resolve_path(parent,base_dir); leaf <- basename(path)
 if (leaf == ".") return(pp)
 if (leaf == "..") return(dirname(pp))
 file.path(pp,leaf)
}
read_records <- function(path) {
 if (!file_test("-f",path) || file.access(path,4L)!=0L || is.na(file.info(path)$size) || file.info(path)$size<=0) die("Readable nonempty configuration required.")
 lines <- readLines(path,warn=FALSE)
 lines <- lines[!grepl("^\\s*#",lines)]
 # read.dcf silently accepts duplicate fields: reject them before parsing.
 blocks <- split(lines,cumsum(!nzchar(trimws(lines))))
 for (b in blocks) {
  starts <- b[grepl("^[^[:space:]#:][^:]*:",b)]
  keys <- sub(":.*$","",starts)
  if (anyDuplicated(keys)) die("Repeated DCF field.")
 }
 con <- textConnection(lines); on.exit(close(con))
 x <- read.dcf(con)
 lapply(seq_len(nrow(x)),function(i) {
  ok <- !is.na(x[i,]); r <- setNames(as.character(x[i,ok]),colnames(x)[ok])
  r[] <- gsub("[\\r\\n]+[ \\t]*"," ",r,perl=TRUE); r
 })
}
read_stage_config <- function(path,fields) {
 records <- read_records(path)
 if (length(records)!=1L || !setequal(names(records[[1]]),fields) || length(records[[1]])!=length(fields)) die("Worker DCF fields differ from the stage contract.")
 cfg <- records[[1]][fields]
 if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) die("Replace all placeholders.")
 if ("edge_id" %in% fields && !cfg[["edge_id"]] %in% c("E1","E3")) die("edge_id must be E1 or E3.")
 cfg
}
bind_stage_helper <- function(fn,scope) {environment(fn) <- scope; fn}
prepare_stage <- function(cfg,context,check_only,scope) {
 fields <- names(cfg); path_fields <- setdiff(fields,"edge_id")
 paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(context$config))
 args <- c("--config",context$config,if(check_only)"--check-config")
 list2env(list(fields=fields,required_fields=fields,path_fields=path_fields,
  paths=paths,config_path=context$config,args=args,
  flags=if(check_only)"--check-config" else character()),envir=scope)
 # One retained ORA implementation explicitly looks up paths in .GlobalEnv.
 # 每个阶段均在新进程运行；此绑定保持原脚本的全局路径查询行为。
 assign("paths",paths,envir=.GlobalEnv)
 invisible(NULL)
}
