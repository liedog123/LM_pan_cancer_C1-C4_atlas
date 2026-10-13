#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 16_build_bassez_balanced_cellchat_input.R --config PATH [--check-config]\n       Rscript 16_build_bassez_balanced_cellchat_input.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("cell_metadata", "sample_patient", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("cell_metadata", "sample_patient")]
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
  library(digest)
})

script_name <- "02_freeze_bassez_balanced_cellchat_input_v0.2.R"
start_time <- Sys.time()

out_dir <- unname(cfg_paths["output_dir"])
log_dir <- file.path(out_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  cell_metadata = unname(cfg_paths["cell_metadata"]),
  sample_patient = unname(cfg_paths["sample_patient"])
)

major_types <- c("Malignant", "T_cell", "Myeloid_cell", "Fibroblast")
states <- c("C1", "C2", "C3", "C4")
min_cells_per_block <- 20L
patients_per_state <- 5L
n_replicates <- 30L
base_seed <- 20260811L

log_event <- function(level, msg, ...) {
  extra <- list(...)
  extra_text <- if (length(extra)) paste(sprintf("%s=%s", names(extra), vapply(extra, as.character, character(1))), collapse = " ") else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, msg,
              if (nzchar(extra_text)) " " else "", extra_text))
  flush.console()
}

write_tsv <- function(x, filename, compress = FALSE) {
  target <- file.path(out_dir, filename)
  fwrite(x, target, sep = "\t", quote = FALSE, na = "NA", compress = if (compress) "gzip" else "none")
  target
}

main <- function() {
  log_event("START", "Bassez balanced CellChat input freeze begin", output = out_dir)

  meta <- fread(paths$cell_metadata, select = c("cell_id", "study", "sample_unit_key", "cell_types_update"))
  bmeta <- meta[study == "Bassez" & cell_types_update %in% major_types]
  if (!nrow(bmeta)) stop("No Bassez major-cell-type metadata found")
  if (anyDuplicated(bmeta$cell_id)) stop("Duplicated Bassez cell IDs detected")

  sample_map <- fread(paths$sample_patient)
  required_map <- c("sample_unit_key", "patient_key", "study", "sample_C_subtype")
  if (!all(required_map %in% names(sample_map))) stop("Sample-patient mapping lacks required fields")
  bsample <- unique(sample_map[study == "Bassez", ..required_map])

  coverage_long <- bmeta[, .(n_cells = .N), by = .(sample_unit_key, cell_types_update)]
  coverage_wide <- dcast(coverage_long, sample_unit_key ~ cell_types_update, value.var = "n_cells", fill = 0L)
  for (ct in major_types) if (!(ct %in% names(coverage_wide))) coverage_wide[, (ct) := 0L]
  coverage_wide[, eligible_all_four_ge20 := Reduce(`&`, lapply(.SD, `>=`, min_cells_per_block)), .SDcols = major_types]
  ecosystem <- merge(bsample, coverage_wide, by = "sample_unit_key", all.x = TRUE)
  for (ct in major_types) set(ecosystem, which(is.na(ecosystem[[ct]])), ct, 0L)
  ecosystem[is.na(eligible_all_four_ge20), eligible_all_four_ge20 := FALSE]
  ecosystem <- ecosystem[sample_C_subtype %in% states]
  setorder(ecosystem, sample_C_subtype, patient_key, sample_unit_key)
  write_tsv(ecosystem, "eligible_ecosystem_coverage.tsv")

  eligible <- ecosystem[eligible_all_four_ge20 == TRUE]
  eligibility_summary <- eligible[, .(
    n_eligible_samples = uniqueN(sample_unit_key),
    n_eligible_patients = uniqueN(patient_key)
  ), by = sample_C_subtype][CJ(sample_C_subtype = states), on = "sample_C_subtype"]
  eligibility_summary[is.na(n_eligible_samples), `:=`(n_eligible_samples = 0L, n_eligible_patients = 0L)]
  eligibility_summary <- eligibility_summary[match(states, sample_C_subtype)]
  write_tsv(eligibility_summary, "eligibility_summary.tsv")

  if (any(eligibility_summary$n_eligible_patients < patients_per_state)) {
    stop("At least one state has fewer than five eligible independent patients")
  }
  log_event("CHECKPOINT", "eligible four-compartment ecosystems frozen", step = "EB0_1",
            state_patients = paste(eligibility_summary$n_eligible_patients, collapse = ","))

  # 每个replicate内先按患者选一个sample，再按sample×cell type固定抽20个细胞。
  design_rows <- vector("list", n_replicates * length(states) * patients_per_state * length(major_types))
  idx <- 0L
  for (rep_id in seq_len(n_replicates)) {
    for (state_id in seq_along(states)) {
      state <- states[[state_id]]
      state_eligible <- eligible[sample_C_subtype == state]
      patient_ids <- sort(unique(state_eligible$patient_key))
      set.seed(base_seed + rep_id * 1000L + state_id * 100L)
      selected_patients <- if (length(patient_ids) == patients_per_state) patient_ids else sample(patient_ids, patients_per_state, replace = FALSE)

      for (patient_pos in seq_along(selected_patients)) {
        patient_id <- selected_patients[[patient_pos]]
        sample_ids <- sort(unique(state_eligible[patient_key == patient_id, sample_unit_key]))
        set.seed(base_seed + rep_id * 100000L + state_id * 10000L + patient_pos * 100L)
        selected_sample <- if (length(sample_ids) == 1L) sample_ids[[1]] else sample(sample_ids, 1L)

        for (type_id in seq_along(major_types)) {
          cell_type <- major_types[[type_id]]
          block_cells <- bmeta[sample_unit_key == selected_sample & cell_types_update == cell_type, cell_id]
          if (length(block_cells) < min_cells_per_block) stop("Eligible block unexpectedly has fewer than 20 cells")
          block_seed <- base_seed + rep_id * 1000000L + state_id * 100000L + patient_pos * 1000L + type_id
          set.seed(block_seed)
          selected_cells <- sample(sort(block_cells), min_cells_per_block, replace = FALSE)
          idx <- idx + 1L
          design_rows[[idx]] <- data.table(
            replicate = rep_id,
            sample_C_subtype = state,
            patient_key = patient_id,
            sample_unit_key = selected_sample,
            cell_types_update = cell_type,
            block_seed = block_seed,
            cell_id = selected_cells
          )
        }
      }
    }
  }
  design <- rbindlist(design_rows[seq_len(idx)])
  design[, `:=`(
    state_order = match(sample_C_subtype, states),
    celltype_order = match(cell_types_update, major_types)
  )]
  setorder(design, replicate, state_order, patient_key, celltype_order, cell_id)
  design[, c("state_order", "celltype_order") := NULL]
  design_file <- write_tsv(design, "balanced_cell_manifest_30rep.tsv.gz", compress = TRUE)

  block_summary <- design[, .(
    n_cells = .N,
    n_unique_cells = uniqueN(cell_id)
  ), by = .(replicate, sample_C_subtype, patient_key, sample_unit_key, cell_types_update)]
  replicate_summary <- design[, .(
    n_cells = .N,
    n_patients = uniqueN(patient_key),
    n_samples = uniqueN(sample_unit_key),
    n_celltypes = uniqueN(cell_types_update)
  ), by = .(replicate, sample_C_subtype)]
  write_tsv(block_summary, "balanced_block_summary.tsv.gz", compress = TRUE)
  write_tsv(replicate_summary, "balanced_replicate_summary.tsv")

  expected_rows <- n_replicates * length(states) * patients_per_state * length(major_types) * min_cells_per_block
  checks <- data.table(
    check = c(
      "all_states_have_five_eligible_patients",
      "manifest_expected_rows",
      "each_block_exactly_20_cells",
      "each_block_unique_cells",
      "each_replicate_state_five_patients",
      "each_replicate_state_four_celltypes",
      "each_replicate_state_400_cells",
      "cell_ids_exist_in_metadata"
    ),
    pass = c(
      all(eligibility_summary$n_eligible_patients >= patients_per_state),
      nrow(design) == expected_rows,
      all(block_summary$n_cells == min_cells_per_block),
      all(block_summary$n_unique_cells == min_cells_per_block),
      all(replicate_summary$n_patients == patients_per_state),
      all(replicate_summary$n_celltypes == length(major_types)),
      all(replicate_summary$n_cells == patients_per_state * length(major_types) * min_cells_per_block),
      all(design$cell_id %chin% bmeta$cell_id)
    )
  )
  gate <- if (all(checks$pass)) "PASS_EB0_GO_EB1" else "STOP_EB0_DESIGN_OR_COVERAGE"
  write_tsv(checks, "balanced_input_checks.tsv")
  write_tsv(data.table(
    gate = gate,
    n_replicates = n_replicates,
    states = paste(states, collapse = ","),
    patients_per_state = patients_per_state,
    celltypes = paste(major_types, collapse = ","),
    cells_per_sample_celltype = min_cells_per_block,
    cells_per_state_replicate = patients_per_state * length(major_types) * min_cells_per_block,
    total_manifest_rows = nrow(design),
    manifest_sha256 = digest(file = design_file, algo = "sha256"),
    next_step = if (gate == "PASS_EB0_GO_EB1") "Run one CellChat smoke replicate for all four states" else "Stop and inspect failed balanced-input checks"
  ), "balanced_input_gate.tsv")

  capture.output(sessionInfo(), file = file.path(log_dir, "eb0_biology_first_session_info_v0.2.txt"))
  runtime <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  log_event("DONE", "Bassez balanced CellChat input freeze complete", gate = gate,
            rows = nrow(design), runtime_sec = runtime, output = out_dir)
}

tryCatch(
  main(),
  error = function(e) {
    fwrite(data.table(gate = "ERROR_EB0_SCRIPT", message = conditionMessage(e)),
           file.path(out_dir, "balanced_input_error.tsv"), sep = "\t", quote = FALSE)
    log_event("ERROR", conditionMessage(e), output = out_dir)
    quit(status = 1)
  }
)
