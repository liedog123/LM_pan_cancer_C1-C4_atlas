#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 17_run_bassez_cellchat_rep1.R --config PATH [--check-config]\n       Rscript 17_run_bassez_cellchat_rep1.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("seurat_object", "balanced_cell_manifest", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("seurat_object", "balanced_cell_manifest")]
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

# ---- Retained source body; only documented path/read adapters ----

options(stringsAsFactors = FALSE, width = 180)

suppressPackageStartupMessages({
  library(data.table)
  library(SeuratObject)
  library(CellChat)
})

script_name <- "03_run_bassez_balanced_cellchat_smoke_v0.2.R"
start_time <- Sys.time()


out_dir <- unname(cfg_paths["output_dir"])
object_dir <- file.path(out_dir, "objects")
table_dir <- file.path(out_dir, "tables")
log_dir <- file.path(out_dir, "logs")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

seurat_rds <- unname(cfg_paths["seurat_object"])
manifest_path <- unname(cfg_paths["balanced_cell_manifest"])
states <- c("C1", "C2", "C3", "C4")
celltype_order <- c("Malignant", "T_cell", "Myeloid_cell", "Fibroblast")
smoke_replicate <- 1L
smoke_nboot <- 20L
base_seed <- 20260811L

log_event <- function(level, msg, ...) {
  extra <- list(...)
  extra_text <- if (length(extra)) paste(sprintf("%s=%s", names(extra), vapply(extra, as.character, character(1))), collapse = " ") else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, msg,
              if (nzchar(extra_text)) " " else "", extra_text))
  flush.console()
}

matrix_to_long <- function(mat, value_name) {
  out <- as.data.table(as.table(mat))
  setnames(out, c("source", "target", value_name))
  out
}

run_state <- function(obj, state_manifest, state, db_use) {
  selected_cells <- unique(state_manifest$cell_id)
  if (length(selected_cells) != 400L) stop(sprintf("%s does not contain exactly 400 selected cells", state))
  if (!all(selected_cells %in% colnames(obj))) stop(sprintf("%s has cells absent from Seurat object", state))

  expr <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")[, selected_cells, drop = FALSE]
  meta <- unique(state_manifest[, .(cell_id, cell_types_update)])
  meta <- meta[match(colnames(expr), cell_id)]
  if (anyNA(meta$cell_id)) stop(sprintf("%s metadata order failed", state))
  meta[, cell_types_update := factor(cell_types_update, levels = celltype_order)]
  meta_df <- as.data.frame(meta[, .(cell_types_update)])
  rownames(meta_df) <- meta$cell_id

  cellchat <- CellChat::createCellChat(object = expr, meta = meta_df, group.by = "cell_types_update", do.sparse = TRUE)
  cellchat@DB <- db_use
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = TRUE,
    population.size = FALSE,
    nboot = smoke_nboot,
    seed.use = base_seed + match(state, states)
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 20L)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  object_path <- file.path(object_dir, sprintf("CellChat_%s_rep001.rds", state))
  saveRDS(cellchat, object_path, compress = "xz")

  count_long <- matrix_to_long(cellchat@net$count, "interaction_count")
  strength_long <- matrix_to_long(cellchat@net$weight, "interaction_strength")
  route <- merge(count_long, strength_long, by = c("source", "target"), all = TRUE)
  route[, sample_C_subtype := state]
  fwrite(route, file.path(table_dir, sprintf("%s_route_network.tsv", state)), sep = "\t", quote = FALSE)

  pathway_flow <- data.table()
  if (length(cellchat@netP$prob)) {
    prob_arr <- cellchat@netP$prob
    flow <- apply(prob_arr, 3L, sum, na.rm = TRUE)
    pathway_flow <- data.table(sample_C_subtype = state, pathway = names(flow), information_flow = as.numeric(flow))
    fwrite(pathway_flow, file.path(table_dir, sprintf("%s_pathway_information_flow.tsv", state)), sep = "\t", quote = FALSE)
  }

  lr <- tryCatch(CellChat::subsetCommunication(cellchat, slot.name = "net"), error = function(e) data.frame())
  lr_dt <- as.data.table(lr)
  if (nrow(lr_dt)) lr_dt[, sample_C_subtype := state]
  fwrite(lr_dt, file.path(table_dir, sprintf("%s_lr_communication.tsv.gz", state)), sep = "\t", quote = FALSE, compress = "gzip")

  data.table(
    sample_C_subtype = state,
    n_cells = ncol(expr),
    n_celltypes = length(unique(meta$cell_types_update)),
    n_interactions = sum(cellchat@net$count, na.rm = TRUE),
    total_interaction_strength = sum(cellchat@net$weight, na.rm = TRUE),
    n_pathways = nrow(pathway_flow),
    n_lr_rows = nrow(lr_dt),
    object_path = object_path,
    object_size_bytes = file.info(object_path)$size
  )
}

main <- function() {
  log_event("START", "Bassez balanced CellChat four-state smoke run begin", replicate = smoke_replicate, nboot = smoke_nboot)
  if (!file.exists(seurat_rds) || !file.exists(manifest_path)) stop("Required Seurat object or balanced manifest is missing")

  manifest <- fread(manifest_path)
  manifest <- manifest[replicate == smoke_replicate & sample_C_subtype %in% states]
  if (nrow(manifest) != 1600L) stop("Smoke manifest must contain exactly 1,600 cells across four states")

  data("CellChatDB.human", package = "CellChat")
  db_use <- CellChatDB.human
  wanted_annotations <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
  db_use$interaction <- db_use$interaction[db_use$interaction$annotation %in% wanted_annotations, , drop = FALSE]
  if (!nrow(db_use$interaction)) stop("CellChatDB subset is empty")

  obj <- readRDS(seurat_rds)
  if (!("RNA" %in% names(obj@assays))) stop("RNA assay is absent")
  if (!("data" %in% SeuratObject::Layers(obj[["RNA"]]))) stop("RNA normalized data layer is absent")
  log_event("CHECKPOINT", "Seurat object loaded", step = "EB1_1", n_object_cells = ncol(obj))

  state_results <- vector("list", length(states))
  for (i in seq_along(states)) {
    state <- states[[i]]
    log_event("CHECKPOINT", "CellChat state run begin", step = paste0("EB1_", i + 1L), state = state)
    state_results[[i]] <- run_state(obj, manifest[sample_C_subtype == state], state, db_use)
    log_event("CHECKPOINT", "CellChat state run complete", step = paste0("EB1_", i + 1L), state = state,
              interactions = state_results[[i]]$n_interactions, pathways = state_results[[i]]$n_pathways)
  }
  summary <- rbindlist(state_results)
  fwrite(summary, file.path(out_dir, "cellchat_smoke_state_summary.tsv"), sep = "\t", quote = FALSE)

  checks <- data.table(
    check = c("all_four_states_completed", "all_states_400_cells", "all_states_four_celltypes", "all_states_nonzero_interactions", "all_states_nonzero_pathways", "all_objects_nonempty"),
    pass = c(
      nrow(summary) == 4L && setequal(summary$sample_C_subtype, states),
      all(summary$n_cells == 400L),
      all(summary$n_celltypes == 4L),
      all(summary$n_interactions > 0L),
      all(summary$n_pathways > 0L),
      all(summary$object_size_bytes > 0L)
    )
  )
  gate <- if (all(checks$pass)) "PASS_EB1_SMOKE_GO_10REP_PILOT" else "STOP_EB1_SMOKE"
  fwrite(checks, file.path(out_dir, "cellchat_smoke_checks.tsv"), sep = "\t", quote = FALSE)
  fwrite(data.table(gate = gate, replicate = smoke_replicate, nboot = smoke_nboot,
                    next_step = if (gate == "PASS_EB1_SMOKE_GO_10REP_PILOT") "Run 10 balanced CellChat replicates and aggregate consensus landscape" else "Stop and inspect failed state/network outputs"),
         file.path(out_dir, "cellchat_smoke_gate.tsv"), sep = "\t", quote = FALSE)
  capture.output(sessionInfo(), file = file.path(log_dir, "eb1_cellchat_smoke_session_info_v0.2.txt"))
  runtime <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  log_event("DONE", "Bassez balanced CellChat four-state smoke run complete", gate = gate, runtime_sec = runtime, output = out_dir)
}

tryCatch(
  main(),
  error = function(e) {
    fwrite(data.table(gate = "ERROR_EB1_SMOKE", message = conditionMessage(e)),
           file.path(out_dir, "cellchat_smoke_error.tsv"), sep = "\t", quote = FALSE)
    capture.output(sessionInfo(), file = file.path(log_dir, "eb1_cellchat_smoke_session_info_v0.2.txt"))
    log_event("ERROR", conditionMessage(e), output = out_dir)
    quit(status = 1)
  }
)

