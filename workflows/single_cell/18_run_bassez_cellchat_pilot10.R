#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 18_run_bassez_cellchat_pilot10.R --config PATH [--check-config]\n       Rscript 18_run_bassez_cellchat_pilot10.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("seurat_object", "balanced_cell_manifest", "rep1_objects_dir", "output_dir")
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
input_dirs <- unname(paths[c("rep1_objects_dir")])
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

# 普通执行只定位rep1的固定对象basename，保留源缺失对象stop。
rep1_input <- function(leaf) {
  if(length(leaf)!=1L || is.na(leaf) || !nzchar(leaf) || basename(leaf)!=leaf ||
     leaf %in% c(".", "..")) stop("Invalid replicate-1 object basename.")
  candidate <- resolve_path(file.path(cfg_paths["rep1_objects_dir"],leaf), dirname(config_path))
  if(!is_within(candidate,cfg_paths["rep1_objects_dir"]) ||
     identical(candidate,cfg_paths["rep1_objects_dir"])) stop("Replicate-1 object escapes configured directory.")
  if(!file.exists(candidate) && !dir.exists(candidate)) return(candidate)
  if(!file_test("-f",candidate) || file.access(candidate,4L)!=0L ||
     is.na(file.info(candidate)$size) || file.info(candidate)$size<=0)
    stop("Selected replicate-1 object must be a readable nonempty regular file.")
  candidate
}

# ---- Retained source body; only documented path/read adapters ----

options(stringsAsFactors = FALSE, width = 180)

suppressPackageStartupMessages({
  library(data.table)
  library(SeuratObject)
  library(CellChat)
})

script_name <- "04_run_bassez_balanced_cellchat_pilot10_v0.2.R"
start_time <- Sys.time()



out_dir <- unname(cfg_paths["output_dir"])
object_dir <- file.path(out_dir, "objects")
replicate_dir <- file.path(out_dir, "replicate_tables")
consensus_dir <- file.path(out_dir, "consensus_tables")
log_dir <- file.path(out_dir, "logs")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(consensus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

seurat_rds <- unname(cfg_paths["seurat_object"])
manifest_path <- unname(cfg_paths["balanced_cell_manifest"])
states <- c("C1", "C2", "C3", "C4")
celltype_order <- c("Malignant", "T_cell", "Myeloid_cell", "Fibroblast")
pilot_replicates <- 1:10
nboot <- 20L
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

extract_state_tables <- function(cellchat, state, rep_id) {
  route <- merge(
    matrix_to_long(cellchat@net$count, "interaction_count"),
    matrix_to_long(cellchat@net$weight, "interaction_strength"),
    by = c("source", "target"), all = TRUE
  )
  route[, `:=`(replicate = rep_id, sample_C_subtype = state)]

  pathway <- data.table()
  if (length(cellchat@netP$prob)) {
    flow <- apply(cellchat@netP$prob, 3L, sum, na.rm = TRUE)
    pathway <- data.table(replicate = rep_id, sample_C_subtype = state,
                          pathway = names(flow), information_flow = as.numeric(flow))
  }

  lr <- as.data.table(tryCatch(CellChat::subsetCommunication(cellchat, slot.name = "net"), error = function(e) data.frame()))
  if (nrow(lr)) lr[, `:=`(replicate = rep_id, sample_C_subtype = state)]

  summary <- data.table(
    replicate = rep_id,
    sample_C_subtype = state,
    n_interactions = sum(cellchat@net$count, na.rm = TRUE),
    total_interaction_strength = sum(cellchat@net$weight, na.rm = TRUE),
    n_pathways = nrow(pathway),
    n_lr_rows = nrow(lr)
  )
  list(route = route, pathway = pathway, lr = lr, summary = summary)
}

run_state <- function(obj, state_manifest, state, rep_id, db_use) {
  selected_cells <- unique(state_manifest$cell_id)
  if (length(selected_cells) != 400L) stop(sprintf("replicate %d %s does not contain 400 cells", rep_id, state))
  expr <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")[, selected_cells, drop = FALSE]
  meta <- unique(state_manifest[, .(cell_id, cell_types_update)])[match(colnames(expr), cell_id)]
  if (anyNA(meta$cell_id)) stop(sprintf("replicate %d %s metadata match failed", rep_id, state))
  meta[, cell_types_update := factor(cell_types_update, levels = celltype_order)]
  meta_df <- as.data.frame(meta[, .(cell_types_update)])
  rownames(meta_df) <- meta$cell_id

  cellchat <- CellChat::createCellChat(object = expr, meta = meta_df, group.by = "cell_types_update", do.sparse = TRUE)
  cellchat@DB <- db_use
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat, type = "triMean", raw.use = TRUE, population.size = FALSE,
    nboot = nboot, seed.use = base_seed + rep_id * 100L + match(state, states)
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 20L)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file.path(object_dir, sprintf("CellChat_%s_rep%03d.rds", state, rep_id)), compress = "gzip")
  extract_state_tables(cellchat, state, rep_id)
}

read_smoke_rep1 <- function() {
  objects <- lapply(states, function(state) {
    path <- rep1_input(sprintf("CellChat_%s_rep001.rds", state))
    if (!file.exists(path)) stop("Replicate-1 smoke object missing: ", path)
    extract_state_tables(readRDS(path), state, 1L)
  })
  objects
}

main <- function() {
  log_event("START", "Bassez balanced CellChat 10-replicate pilot begin", replicates = "1-10", nboot = nboot)
  manifest <- fread(manifest_path)
  manifest <- manifest[replicate %in% pilot_replicates & sample_C_subtype %in% states]
  if (nrow(manifest) != 16000L) stop("Pilot manifest must contain 16,000 rows")

  data("CellChatDB.human", package = "CellChat")
  wanted_annotations <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
  db_use <- CellChatDB.human
  db_use$interaction <- db_use$interaction[db_use$interaction$annotation %in% wanted_annotations, , drop = FALSE]

  result_list <- read_smoke_rep1()
  log_event("CHECKPOINT", "replicate-1 smoke objects imported", step = "EB1P_1")

  obj <- readRDS(seurat_rds)
  if (!("data" %in% SeuratObject::Layers(obj[["RNA"]]))) stop("RNA normalized data layer is absent")
  log_event("CHECKPOINT", "Seurat object loaded", step = "EB1P_2", n_object_cells = ncol(obj))

  idx <- length(result_list)
  for (rep_id in 2:10) {
    for (state in states) {
      log_event("CHECKPOINT", "pilot state run begin", step = "EB1P_RUN", replicate = rep_id, state = state)
      idx <- idx + 1L
      result_list[[idx]] <- run_state(obj, manifest[replicate == rep_id & sample_C_subtype == state], state, rep_id, db_use)
      log_event("CHECKPOINT", "pilot state run complete", step = "EB1P_RUN", replicate = rep_id, state = state,
                interactions = result_list[[idx]]$summary$n_interactions, pathways = result_list[[idx]]$summary$n_pathways)
    }
  }

  route_all <- rbindlist(lapply(result_list, `[[`, "route"), fill = TRUE)
  pathway_all <- rbindlist(lapply(result_list, `[[`, "pathway"), fill = TRUE)
  lr_all <- rbindlist(lapply(result_list, `[[`, "lr"), fill = TRUE)
  summary_all <- rbindlist(lapply(result_list, `[[`, "summary"), fill = TRUE)

  fwrite(route_all, file.path(replicate_dir, "route_network_all_replicates.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")
  fwrite(pathway_all, file.path(replicate_dir, "pathway_information_flow_all_replicates.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")
  fwrite(lr_all, file.path(replicate_dir, "lr_communication_all_replicates.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")
  fwrite(summary_all, file.path(replicate_dir, "state_network_summary_all_replicates.tsv"), sep = "\t", quote = FALSE)

  route_consensus <- route_all[, .(
    median_interaction_count = median(interaction_count, na.rm = TRUE),
    median_interaction_strength = median(interaction_strength, na.rm = TRUE),
    recovery_fraction = mean(interaction_count > 0, na.rm = TRUE)
  ), by = .(sample_C_subtype, source, target)]

  pathway_grid <- CJ(replicate = pilot_replicates, sample_C_subtype = states,
                     pathway = sort(unique(pathway_all$pathway)), unique = TRUE)
  pathway_complete <- merge(pathway_grid, pathway_all, by = c("replicate", "sample_C_subtype", "pathway"), all.x = TRUE)
  pathway_complete[is.na(information_flow), information_flow := 0]
  pathway_consensus <- pathway_complete[, .(
    median_information_flow = median(information_flow),
    mean_information_flow = mean(information_flow),
    recovery_fraction = mean(information_flow > 0)
  ), by = .(sample_C_subtype, pathway)]

  lr_id_cols <- intersect(c("source", "target", "ligand", "receptor", "interaction_name", "interaction_name_2", "pathway_name", "annotation"), names(lr_all))
  lr_present <- unique(lr_all[, c("replicate", "sample_C_subtype", "prob", lr_id_cols), with = FALSE])
  lr_keys <- unique(lr_present[, c("sample_C_subtype", lr_id_cols), with = FALSE])
  replicate_grid <- CJ(replicate = pilot_replicates, unique = TRUE)
  replicate_grid[, join_key := 1L]
  lr_keys[, join_key := 1L]
  lr_grid <- merge(replicate_grid, lr_keys, by = "join_key", allow.cartesian = TRUE)[, join_key := NULL]
  lr_complete <- merge(lr_grid, lr_present, by = c("replicate", "sample_C_subtype", lr_id_cols), all.x = TRUE)
  lr_complete[is.na(prob), prob := 0]
  lr_consensus <- lr_complete[, .(
    median_probability = median(prob),
    mean_probability = mean(prob),
    recovery_fraction = mean(prob > 0)
  ), by = c("sample_C_subtype", lr_id_cols)]

  state_consensus <- summary_all[, .(
    median_n_interactions = median(n_interactions),
    range_n_interactions = paste(range(n_interactions), collapse = "-"),
    median_total_strength = median(total_interaction_strength),
    range_total_strength = paste(signif(range(total_interaction_strength), 5), collapse = "-"),
    median_n_pathways = median(n_pathways),
    range_n_pathways = paste(range(n_pathways), collapse = "-")
  ), by = sample_C_subtype]

  fwrite(route_consensus, file.path(consensus_dir, "route_network_consensus.tsv"), sep = "\t", quote = FALSE)
  fwrite(pathway_consensus, file.path(consensus_dir, "pathway_information_flow_consensus.tsv"), sep = "\t", quote = FALSE)
  fwrite(lr_consensus, file.path(consensus_dir, "lr_communication_consensus.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")
  fwrite(state_consensus, file.path(consensus_dir, "state_network_consensus_summary.tsv"), sep = "\t", quote = FALSE)

  checks <- data.table(
    check = c("forty_state_runs", "all_runs_nonzero_interactions", "all_runs_nonzero_pathways", "all_consensus_tables_nonempty"),
    pass = c(nrow(summary_all) == 40L, all(summary_all$n_interactions > 0), all(summary_all$n_pathways > 0),
             nrow(route_consensus) > 0 && nrow(pathway_consensus) > 0 && nrow(lr_consensus) > 0)
  )
  gate <- if (all(checks$pass)) "PASS_EB1_PILOT10_READY_FOR_BIOLOGICAL_REVIEW" else "STOP_EB1_PILOT10"
  fwrite(checks, file.path(out_dir, "pilot10_checks.tsv"), sep = "\t", quote = FALSE)
  fwrite(data.table(gate = gate, n_replicates = 10L, n_state_runs = nrow(summary_all), nboot = nboot,
                    next_step = if (grepl("PASS", gate)) "Review recurrent pathways/routes/LR and design Fig6-E preview" else "Stop and inspect failed pilot outputs"),
         file.path(out_dir, "pilot10_gate.tsv"), sep = "\t", quote = FALSE)
  capture.output(sessionInfo(), file = file.path(log_dir, "eb1_cellchat_pilot10_session_info_v0.2.txt"))
  runtime <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  log_event("DONE", "Bassez balanced CellChat 10-replicate pilot complete", gate = gate, runtime_sec = runtime, output = out_dir)
}

tryCatch(
  main(),
  error = function(e) {
    fwrite(data.table(gate = "ERROR_EB1_PILOT10", message = conditionMessage(e)),
           file.path(out_dir, "pilot10_error.tsv"), sep = "\t", quote = FALSE)
    capture.output(sessionInfo(), file = file.path(log_dir, "eb1_cellchat_pilot10_session_info_v0.2.txt"))
    log_event("ERROR", conditionMessage(e), output = out_dir)
    quit(status = 1)
  }
)
