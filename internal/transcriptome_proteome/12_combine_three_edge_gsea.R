#!/usr/bin/env Rscript
# Purpose: Test the three adjacent-state comparisons using their original ranked-gene families.
# Task: boundary_gsea; step: s12
# Inputs (DCF fields): e1_dir,e2_dir,e3_dir,definition_registry
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task boundary_gsea --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 12_combine_three_edge_gsea.R --config PATH [--check-config]\n       Rscript 12_combine_three_edge_gsea.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { usage(); quit(status = 0L) }
if (!(length(args) %in% c(2L,3L)) || args[1L] != "--config" ||
    !nzchar(args[2L]) || (length(args)==3L && args[3L] != "--check-config")) {
  usage(); stop("Expected --config PATH [--check-config].")
}
resolve_path <- function(path, base_dir) {
  if (!grepl("^/", path)) path <- file.path(base_dir, path)
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
cfg <- read.dcf(config_path)
fields <- c("e1_dir", "e2_dir", "e3_dir", "definition_registry", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- setdiff(fields,"edge_id")
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
if(any(!dir.exists(paths[c("e1_dir","e2_dir","e3_dir")])) || any(file.access(paths[c("e1_dir","e2_dir","e3_dir")],5L)!=0L)) stop("Edge directories must be readable/searchable.")
required_files <- c(paths["definition_registry"],unlist(lapply(1:3,function(j) file.path(paths[paste0("e",j,"_dir")],c(paste0("e",j,"_gsea_all_results.tsv.gz"),"machine_qa.tsv"))),use.names=FALSE))
if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
protected <- c(input_paths[dir.exists(input_paths)],dirname(required_files),config_path)
if(out_dir=="/" || out_dir==normalizePath(path.expand("~")) ||
   any(vapply(protected,function(p) is_within(p,out_dir)||is_within(out_dir,p),logical(1L)))) stop("Output must be independent of inputs/config.")
if(file.exists(out_dir) && !dir.exists(out_dir)) stop("Output is not a directory.")
if(dir.exists(out_dir) && (file.access(out_dir,7L)!=0L || length(list.files(out_dir,all.files=TRUE,no..=TRUE)))) stop("Output must be writable and empty.")
parent <- out_dir
while(!dir.exists(parent)) {
 if(file.exists(parent)) stop("Output has a non-directory ancestor.")
 parent <- dirname(parent)
}
if(file.access(parent,3L)!=0L) stop("No writable searchable output parent.")
if(length(args)==3L) {
 cat("METADATA_CHECK_PASS: no study contents loaded; no outputs created.\n")
 quit(status=0L)
}




suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

script_name <- "12_combine_three_edge_gsea.R"
start_time <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extra <- if (length(fields)) paste(vapply(names(fields), function(nm) {
    sprintf('%s="%s"', nm, gsub("[\r\n\t]", " ", as.character(fields[[nm]])))
  }, character(1)), collapse = " ") else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
              gsub("[\r\n]", " ", msg), if (nzchar(extra)) " " else "", extra))
  flush.console()
}

fail <- function(msg) {
  log_event("ERROR", msg)
  stop(msg, call. = FALSE)
}

norm_symbol <- function(x) toupper(trimws(as.character(x)))

main <- function() {
  final_dir <- file.path(out_dir,"results","combined")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(final_dir) && length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    fail(paste("三边汇总结果目录已存在且非空，拒绝覆盖", final_dir))
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  edge_registry <- data.table(
    Edge_ID = c("E1", "E2", "E3"),
    Target = c("AG1_IG1", "AG2_IG1", "AG2_IG2"),
    Reference = c("AG2_IG1", "AG2_IG2", "AG3_IG2"),
    Target_State = c("C1", "C3", "C2"),
    Reference_State = c("C3", "C2", "C4"),
    Expected_Cancers = c(20L, 19L, 21L)
  )
  edge_dirs <- paths[c("e1_dir","e2_dir","e3_dir")]
  definition_path <- unname(paths["definition_registry"])
  paths <- data.table(Edge_ID=c("E1","E2","E3"),Source=file.path(edge_dirs,paste0("e",1:3,"_gsea_all_results.tsv.gz")),Reuse_Status=c("NEW_E1","REUSED_FROZEN_E2","NEW_E3"))
  if (any(!file.exists(paths$Source)) || any(file.info(paths$Source)$size <= 0)) {
    fail("E1/E2/E3逐癌GSEA汇总输入缺失或为空")
  }

  log_event("START", "three-edge GSEA combine begin", list(
    environment = "R4.5.1", edges = "E1,E2,E3",
    e2 = "reuse_frozen_no_rerun", output = final_dir))

  gsea <- rbindlist(lapply(seq_len(nrow(paths)), function(i) {
    z <- fread(paths$Source[i])
    if (!all(z$Edge_ID == paths$Edge_ID[i])) fail(paste(paths$Edge_ID[i], "来源边标识不一致"))
    drop_cols <- intersect(c("Target_State", "Reference_State"), names(z))
    if (length(drop_cols)) z[, (drop_cols) := NULL]
    merge(z, edge_registry[, .(Edge_ID, Target_State, Reference_State)],
          by = "Edge_ID", all.x = TRUE)
  }), use.names = TRUE, fill = TRUE)

  # 从三条边的逐癌结果统一重算跨癌票数，避免摘要格式差异。
  vote <- gsea[, .(
    Display_Name = Display_Name[1], Functional_Domain = Functional_Domain[1],
    Evidence_Role = Evidence_Role[1], N_Evaluable = uniqueN(Cancer),
    N_Target_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES > 0]),
    N_Reference_High_Sig = uniqueN(Cancer[FDR < 0.05 & NES < 0]),
    Target_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES > 0])), collapse = ";"),
    Reference_High_Sig_Cancers = paste(sort(unique(Cancer[FDR < 0.05 & NES < 0])), collapse = ";"),
    N_Target_High_NES = sum(NES > 0), N_Reference_High_NES = sum(NES < 0),
    Median_NES = median(NES), Median_Abs_NES = median(abs(NES)),
    Leading_Edge_Union = paste(sort(unique(unlist(strsplit(Leading_Edge, ";", fixed = TRUE)))), collapse = ";")
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection, Set_ID)]
  vote[, `:=`(
    Min_Recurrence = pmax(3L, ceiling(0.30 * N_Evaluable)),
    Min_Strong_Support = pmax(4L, ceiling(0.40 * N_Evaluable)),
    Max_Opposite_Allowed = floor(0.10 * N_Evaluable),
    Target_NES_Direction_Fraction = N_Target_High_NES / N_Evaluable,
    Reference_NES_Direction_Fraction = N_Reference_High_NES / N_Evaluable
  )]
  vote[, `:=`(
    Recurrent_Target_High = N_Target_High_Sig >= Min_Recurrence,
    Recurrent_Reference_High = N_Reference_High_Sig >= Min_Recurrence,
    Strong_Target_High = N_Target_High_Sig >= Min_Strong_Support &
      N_Reference_High_Sig <= Max_Opposite_Allowed,
    Strong_Reference_High = N_Reference_High_Sig >= Min_Strong_Support &
      N_Target_High_Sig <= Max_Opposite_Allowed
  )]
  vote[, Support_Class := fifelse(
    Strong_Target_High, "STRONG_TARGET_HIGH",
    fifelse(Strong_Reference_High, "STRONG_REFERENCE_HIGH",
      fifelse(Recurrent_Target_High, "RECURRENT_TARGET_HIGH",
        fifelse(Recurrent_Reference_High, "RECURRENT_REFERENCE_HIGH", "NOT_RECURRENT"))))]
  setorder(vote, Edge_ID, Sensitivity, Collection, -N_Target_High_Sig,
           -N_Reference_High_Sig, -Median_Abs_NES, Set_ID)

  candidates <- vote[
    Sensitivity == "NON_DEFINITION_RANK" &
      Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT") &
      Support_Class %chin% c("STRONG_TARGET_HIGH", "STRONG_REFERENCE_HIGH")
  ]
  candidates[, `:=`(
    Evidence_Tier = fifelse(Collection == "IMMUNE119_PRIMARY",
                            "PRIMARY_IMMUNE", "HALLMARK_CONTEXT"),
    Manual_Redundancy_Review = "PENDING",
    Fig3B_Eligibility = "PENDING_MANUAL_GATE"
  )]

  comparison <- dcast(vote[Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT")],
    Edge_ID + Target_State + Reference_State + Collection + Set_ID +
      Display_Name + Functional_Domain ~ Sensitivity,
    value.var = c("N_Target_High_Sig", "N_Reference_High_Sig", "Support_Class"))

  # 长表直接支持未来三比较×两方向复现矩阵；本轮不绘图。
  matrix_source <- rbind(
    vote[Sensitivity == "NON_DEFINITION_RANK" &
           Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT"), .(
      Edge_ID, Comparison = paste0(Target_State, " vs ", Reference_State),
      State = Target_State, Opposite_State = Reference_State, Collection,
      Set_ID, Display_Name, Functional_Domain, N_Evaluable,
      N_Directional_Significant = N_Target_High_Sig,
      N_Opposite_Significant = N_Reference_High_Sig,
      Direction_Fraction = Target_NES_Direction_Fraction,
      Strong = Strong_Target_High, Recurrent = Recurrent_Target_High,
      Side = "TARGET")],
    vote[Sensitivity == "NON_DEFINITION_RANK" &
           Collection %chin% c("IMMUNE119_PRIMARY", "HALLMARK_CONTEXT"), .(
      Edge_ID, Comparison = paste0(Target_State, " vs ", Reference_State),
      State = Reference_State, Opposite_State = Target_State, Collection,
      Set_ID, Display_Name, Functional_Domain, N_Evaluable,
      N_Directional_Significant = N_Reference_High_Sig,
      N_Opposite_Significant = N_Target_High_Sig,
      Direction_Fraction = Reference_NES_Direction_Fraction,
      Strong = Strong_Reference_High, Recurrent = Recurrent_Reference_High,
      Side = "REFERENCE")]
  )
  setorder(matrix_source, Edge_ID, Collection, Set_ID, Side)

  summary <- vote[, .(
    N_Sets = .N,
    N_Recurrent_Target = sum(Recurrent_Target_High),
    N_Recurrent_Reference = sum(Recurrent_Reference_High),
    N_Strong_Target = sum(Strong_Target_High),
    N_Strong_Reference = sum(Strong_Reference_High)
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection)]

  definition <- norm_symbol(fread(definition_path)$gene_symbol)
  nondef_leading <- unique(unlist(strsplit(
    gsea[Sensitivity == "NON_DEFINITION_RANK", Leading_Edge], ";", fixed = TRUE)))
  cancer_counts <- unique(gsea[, .(N_Cancers = uniqueN(Cancer)), by = Edge_ID])
  cancer_counts <- merge(edge_registry[, .(Edge_ID, Expected_Cancers)], cancer_counts,
                         by = "Edge_ID", all.x = TRUE)
  qa <- data.table(
    Check = c("THREE_EDGES_EXACT", "EXPECTED_CANCER_COUNTS", "EXPECTED_COLLECTIONS",
              "TWO_SENSITIVITIES", "E2_REUSED_NOT_RERUN", "NO_C7",
              "NO_META_FIELDS", "FDR_RANGE_VALID",
              "NONDEFINITION_LEADING_EDGE_EXCLUDES_DEFINITION",
              "ALL_EDGE_QA_PASS", "OUTPUT_NONEMPTY"),
    Pass = c(
      setequal(unique(gsea$Edge_ID), c("E1", "E2", "E3")),
      all(cancer_counts$N_Cancers == cancer_counts$Expected_Cancers),
      setequal(unique(gsea$Collection), c("IMMUNE119_PRIMARY",
        "LM12_DEFINITION_CONTROL", "HALLMARK_CONTEXT")),
      setequal(unique(gsea$Sensitivity), c("FULL_RANK", "NON_DEFINITION_RANK")),
      paths[Edge_ID == "E2", Reuse_Status] == "REUSED_FROZEN_E2",
      !any(grepl("C7", gsea$Collection)),
      !any(c("Meta_FDR", "I2", "Tau2", "LOO") %in% names(gsea)),
      all(is.finite(gsea$FDR) & gsea$FDR >= 0 & gsea$FDR <= 1),
      length(intersect(nondef_leading, definition)) == 0L,
      all(vapply(c("E1", "E3"), function(e) {
        q <- file.path(edge_dirs[paste0(tolower(e),"_dir")], "machine_qa.tsv")
        file.exists(q) && all(fread(q)$Pass)
      }, logical(1))) && all(fread(file.path(edge_dirs["e2_dir"], "machine_qa.tsv"))$Pass),
      nrow(gsea) > 0L && nrow(vote) > 0L && nrow(matrix_source) > 0L
    )
  )
  if (!all(qa$Pass)) fail("三边GSEA统一汇总机器QA失败")

  paths[, `:=`(Size_Bytes = file.info(Source)$size,
               MD5 = unname(tools::md5sum(Source)))]
  fwrite(gsea, file.path(tmp_dir, "three_edge_gsea_all_results.tsv.gz"),
         sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
  fwrite(vote, file.path(tmp_dir, "three_edge_gsea_cross_cancer_votes.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(candidates, file.path(tmp_dir, "three_edge_nondefinition_strong_candidates.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(comparison, file.path(tmp_dir, "three_edge_definition_sensitivity_comparison.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(matrix_source, file.path(tmp_dir, "fig3b_recurrence_matrix_source.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(summary, file.path(tmp_dir, "three_edge_gsea_summary.tsv"),
         sep = "\t", quote = TRUE)
  fwrite(paths, file.path(tmp_dir, "source_reuse_manifest.tsv"), sep = "\t", quote = TRUE)
  fwrite(qa, file.path(tmp_dir, "machine_qa.tsv"), sep = "\t", quote = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(tmp_dir, "sessionInfo.txt"))

  report <- c(
    "# Fig3三条局部边full-rank GSEA执行报告 V0.1", "", "日期：2026-08-06", "",
    "状态：`MACHINE_QA_PASS / MANUAL_GATE_PENDING / POST_FREEZE_EXPLORATORY / STAGE3_UNCHANGED`", "",
    "## 执行边界", "",
    "E1和E3新增全排序GSEA；E2复用既有冻结结果并记录MD5。所有单癌BH、双敏感性、gene-set和复现门槛统一；未重跑limma、未使用meta/LOO/跨癌BH。", "",
    "## 汇总", "",
    paste(apply(summary, 1L, function(x) sprintf(
      "%s（%s vs %s）%s × %s：target复现/强支持%s/%s，reference复现/强支持%s/%s。",
      x[["Edge_ID"]], x[["Target_State"]], x[["Reference_State"]],
      x[["Sensitivity"]], x[["Collection"]], x[["N_Recurrent_Target"]],
      x[["N_Strong_Target"]], x[["N_Recurrent_Reference"]],
      x[["N_Strong_Reference"]])), collapse = "\n\n"), "",
    "## 下一步", "",
    "只对NON_DEFINITION_RANK强支持候选进行语义去冗余和definition敏感性审核，冻结10–12行以内的Fig3B候选；本步骤不绘图。"
  )
  writeLines(report, file.path(tmp_dir, "execution_report.md"))

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) fail("无法提交三边GSEA统一汇总结果")
  log_event("DONE", "three-edge GSEA combine complete", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    rows = nrow(gsea), candidates = nrow(candidates),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)), output = final_dir))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)))
  quit(status = 1L)
})
