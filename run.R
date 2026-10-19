#!/usr/bin/env Rscript
# Purpose: Single entry point for 18 analysis modules and 8 configuration domains.
# Inputs: config/*.dcf copied and completed outside the checkout; see docs/reproduction.md.
# Outputs: stage-specific directories and logs; no outputs in this code checkout.
# 阶段在独立进程中执行，避免状态、随机种子及包加载互相影响。
options(stringsAsFactors=FALSE,warn=1)
self <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
if(length(self)!=1L)stop("Invoke with Rscript.")
root <- dirname(normalizePath(self,mustWork=TRUE))
for(file in c("config","io","statistics","plotting","runtime"))sys.source(file.path(root,"R/common",paste0(file,".R")),envir=.GlobalEnv)
cli_main(commandArgs(TRUE),root)
