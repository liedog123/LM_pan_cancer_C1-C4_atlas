#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 11_run_e1_e3_full_rank_gsea.R --config PATH [--check-config]\n       Rscript 11_run_e1_e3_full_rank_gsea.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("model_dir", "definition_registry", "immune_registry", "immune_membership", "lm12_membership", "hallmark_gmt", "edge_id", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
if(!cfg["edge_id"] %in% c("E1","E3")) stop("edge_id must be E1 or E3.")
path_fields <- setdiff(fields,"edge_id")
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","THCA","UCEC")
if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("model_dir must be readable/searchable.")
required_files <- c(input_paths[names(input_paths)!="model_dir"],file.path(paths["model_dir"],"input_manifest.tsv"),file.path(paths["model_dir"],"individual",paste0(expected_cancers,"_mrna_boundary_models_v0.1_20260806.tsv.gz")))
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




suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
script_name <- "11_run_e1_e3_full_rank_gsea.R"
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

seed_from_keys <- function(...) {
  key <- paste(..., collapse = "|")
  as.integer((20260806 + sum(utf8ToInt(key))) %% .Machine$integer.max)
}

read_gmt <- function(path, collection, evidence_role) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines[nzchar(trimws(lines))], "\t", fixed = TRUE)
  if (any(lengths(parts) < 3L)) fail(paste(collection, "GMT存在少于3列的记录"))
  registry <- rbindlist(lapply(seq_along(parts), function(i) data.table(
    Set_ID = paste0(collection, "::", trimws(parts[[i]][1])),
    Display_Name = trimws(parts[[i]][1]), Collection = collection,
    Functional_Domain = "GENERAL_CONTEXT", Evidence_Role = evidence_role,
    Source_Gene_Count = length(unique(norm_symbol(parts[[i]][-c(1, 2)])))
  )))
  membership <- rbindlist(lapply(seq_along(parts), function(i) {
    genes <- unique(norm_symbol(parts[[i]][-c(1, 2)]))
    genes <- genes[nzchar(genes)]
    data.table(Set_ID = paste0(collection, "::", trimws(parts[[i]][1])), Gene = genes)
  }))
  list(registry = registry, membership = membership)
}

run_fgsea <- function(ranks, registry, membership, cancer, sensitivity, collection,
                      edge_meta) {
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  if (anyDuplicated(names(ranks))) fail(paste(cancer, sensitivity, "rank gene重复"))
  present <- membership[Gene %chin% names(ranks)]
  sizes <- present[, .(Mapped_Size = uniqueN(Gene)), by = Set_ID]
  eligible <- sizes[Mapped_Size >= 10L & Mapped_Size <= 500L]
  meta <- merge(registry, eligible, by = "Set_ID")
  pathways <- split(present[Set_ID %chin% meta$Set_ID, Gene],
                    present[Set_ID %chin% meta$Set_ID, Set_ID])
  pathways <- lapply(pathways, unique)
  if (!length(pathways)) return(data.table())

  set.seed(seed_from_keys(cancer, sensitivity, collection))
  fg <- suppressWarnings(as.data.table(fgseaMultilevel(
    pathways = pathways, stats = ranks, sampleSize = 101L,
    minSize = 10L, maxSize = 500L, eps = 0, scoreType = "std",
    nproc = 1L, nPermSimple = 10000L
  )))
  setnames(fg, "pathway", "Set_ID")
  fg[, GSEA_Method := "fgseaMultilevel"]

  failed <- fg[!is.finite(pval) | !is.finite(NES) | !is.finite(ES), Set_ID]
  if (length(failed)) {
    set.seed(seed_from_keys(cancer, sensitivity, collection, "fallback"))
    fb <- suppressWarnings(as.data.table(fgseaSimple(
      pathways = pathways[failed], stats = ranks, nperm = 100000L,
      minSize = 10L, maxSize = 500L, scoreType = "std", nproc = 1L
    )))
    setnames(fb, "pathway", "Set_ID")
    for (id in failed) {
      i <- which(fg$Set_ID == id)
      z <- fb[Set_ID == id]
      if (length(i) != 1L || nrow(z) != 1L) {
        fail(paste("GSEA fallback无法唯一匹配", cancer, sensitivity, id))
      }
      if (!is.finite(z$pval) || !is.finite(z$NES) || !is.finite(z$ES)) {
        if (!is.finite(fg$ES[i]) || fg$ES[i] == 0) {
          fail(paste("GSEA方向fallback缺少有限ES", cancer, sensitivity, id))
        }
        score_type <- ifelse(fg$ES[i] > 0, "pos", "neg")
        set.seed(seed_from_keys(cancer, sensitivity, collection, id, score_type))
        z <- suppressWarnings(as.data.table(fgseaSimple(
          pathways = pathways[id], stats = ranks, nperm = 100000L,
          minSize = 10L, maxSize = 500L, scoreType = score_type, nproc = 1L
        )))
        setnames(z, "pathway", "Set_ID")
        if (nrow(z) != 1L || !is.finite(z$pval) ||
            !is.finite(z$NES) || !is.finite(z$ES)) {
          fail(paste("GSEA方向fallback仍失败", cancer, sensitivity, id))
        }
        method_label <- paste0("fgseaSimple_nperm100000_", score_type,
                               "_directional_fallback")
      } else {
        method_label <- "fgseaSimple_nperm100000_fallback"
      }
      fg$pval[i] <- z$pval; fg$ES[i] <- z$ES; fg$NES[i] <- z$NES
      fg$size[i] <- z$size; fg$leadingEdge[[i]] <- z$leadingEdge[[1L]]
      fg$log2err[i] <- NA_real_; fg$GSEA_Method[i] <- method_label
    }
  }
  if (any(!is.finite(fg$pval)) || any(!is.finite(fg$NES)) || any(!is.finite(fg$ES))) {
    fail(paste("GSEA仍有非有限结果", cancer, sensitivity, collection))
  }
  fg[, Leading_Edge := vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))]
  fg[, leadingEdge := NULL]
  fg[, FDR := p.adjust(pval, method = "BH")]
  fg[, Direction := fifelse(NES > 0, edge_meta$Positive_Direction,
                            edge_meta$Negative_Direction)]
  fg <- merge(fg, meta, by = "Set_ID", all.x = TRUE)
  fg[, `:=`(Cancer = cancer, Edge_ID = edge_meta$Edge_ID,
            Sensitivity = sensitivity, Collection = collection,
            Target = edge_meta$Target, Reference = edge_meta$Reference,
            Target_State = edge_meta$Target_State,
            Reference_State = edge_meta$Reference_State)]
  setnames(fg, c("pval", "log2err", "size"), c("P_Value", "Log2_Error", "Set_Size_Used"))
  fg[]
}

main <- function() {
  edge_id <- unname(cfg["edge_id"])
  edge_registry <- data.table(
    Edge_ID = c("E1", "E3"),
    Target = c("AG1_IG1", "AG2_IG2"),
    Reference = c("AG2_IG1", "AG3_IG2"),
    Target_State = c("C1", "C2"),
    Reference_State = c("C3", "C4"),
    Positive_Direction = c("C1_LOCAL_HIGH", "C2_LOCAL_HIGH"),
    Negative_Direction = c("C3_LOCAL_HIGH", "C4_LOCAL_HIGH"),
    Expected_Cancers = c(20L, 21L)
  )
  edge_meta <- edge_registry[Edge_ID == edge_id]
  final_dir <- file.path(out_dir,"results","edges",edge_id)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if (dir.exists(final_dir) && length(list.files(final_dir, all.files = TRUE, no.. = TRUE))) {
    fail(paste("输出目录已存在且非空，拒绝覆盖", final_dir))
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(file.path(tmp_dir, "individual"), recursive = TRUE, showWarnings = FALSE)

  model_dir <- unname(paths["model_dir"])
  paths <- as.list(setNames(unname(paths[c("definition_registry","immune_registry","immune_membership","lm12_membership","hallmark_gmt")]),c("definition","immune_registry","immune_long","lm12_long","hallmark")))
  paths$model_manifest <- file.path(model_dir,"input_manifest.tsv")
  if (any(!file.exists(unlist(paths)))) fail(paste("缺少", edge_id, "GSEA冻结输入"))

  manifest <- fread(paths$model_manifest)
  # 只按完整癌种及固定basename重绑路径，不读取历史绝对路径。
  if(!all(c("Cancer","Output") %in% names(manifest)) || nrow(manifest)!=22L || anyNA(manifest$Cancer) || anyDuplicated(manifest$Cancer) || !setequal(manifest$Cancer,expected_cancers)) fail("Complete unique 22-cancer model manifest required.")
  expected_names <- paste0(manifest$Cancer,"_mrna_boundary_models_v0.1_20260806.tsv.gz")
  if(anyNA(manifest$Output) || !identical(basename(manifest$Output),expected_names)) fail("Model manifest basenames do not match Cancer.")
  manifest[, Output := file.path(model_dir,"individual",expected_names)]
  definition <- norm_symbol(fread(paths$definition)$gene_symbol)
  immune_registry <- fread(paths$immune_registry)[, .(
    Set_ID = set_id, Display_Name = display_name, Collection = "IMMUNE119_PRIMARY",
    Functional_Domain = functional_domain, Evidence_Role = evidence_role,
    Source_Gene_Count = source_gene_count
  )]
  immune_membership <- fread(paths$immune_long)[, .(Set_ID = set_id, Gene = norm_symbol(gene))]
  if (uniqueN(immune_registry$Set_ID) != 119L) fail("冻结免疫背景不是119条")

  lm12 <- fread(paths$lm12_long)
  lm12_registry <- lm12[, .(
    Display_Name = pathway[1],
    Collection = "LM12_DEFINITION_CONTROL", Functional_Domain = functional_group[1],
    Evidence_Role = "DEFINITION_POSITIVE_CONTROL",
    Source_Gene_Count = uniqueN(gene)
  ), by = .(Set_ID = feature_id)]
  lm12_membership <- unique(lm12[, .(Set_ID = feature_id, Gene = norm_symbol(gene))])
  if (uniqueN(lm12_registry$Set_ID) != 12L) fail("LM12定义对照不是12条")

  hallmark <- read_gmt(paths$hallmark, "HALLMARK_CONTEXT", "GENERAL_CONTEXT_ONLY")
  registry <- rbindlist(list(immune_registry, lm12_registry, hallmark$registry), use.names = TRUE)
  membership <- rbindlist(list(immune_membership, lm12_membership, hallmark$membership), use.names = TRUE)
  membership <- unique(membership[nzchar(Gene)])

  log_event("START", "single-edge full-rank GSEA begin", list(
    edge = edge_id, environment = "R4.5.1",
    model = paste0("reuse_existing_", edge_id, "_moderated_t"),
    cancers = nrow(manifest), collections = paste(unique(registry$Collection), collapse = ","),
    sensitivities = "FULL_RANK,NON_DEFINITION_RANK", output = final_dir,
    prohibited = "rerun_limma,purity,meta,LOO,cross_cancer_BH,C7"))

  all_results <- list()
  coverage <- list()
  out_manifest <- list()
  for (i in seq_len(nrow(manifest))) {
    cancer <- manifest$Cancer[i]
    x <- fread(manifest$Output[i])
    x <- x[Edge_ID == edge_id & Status == "ESTIMABLE" & is.finite(Statistic)]
    if (!nrow(x)) next
    x[, Gene := norm_symbol(Feature_ID)]
    if (anyDuplicated(x$Gene)) fail(paste(cancer, edge_id, "rank gene重复"))
    full_rank <- setNames(x$Statistic, x$Gene)
    rank_list <- list(
      FULL_RANK = full_rank,
      NON_DEFINITION_RANK = full_rank[!names(full_rank) %chin% definition]
    )
    cancer_rows <- list()
    for (sensitivity in names(rank_list)) {
      for (collection in unique(registry$Collection)) {
        reg0 <- registry[Collection == collection]
        mem0 <- membership[Set_ID %chin% reg0$Set_ID]
        z <- run_fgsea(rank_list[[sensitivity]], reg0, mem0, cancer,
                       sensitivity, collection, edge_meta)
        cancer_rows[[paste(sensitivity, collection)]] <- z
        coverage[[paste(cancer, sensitivity, collection)]] <- data.table(
          Cancer = cancer, Sensitivity = sensitivity, Collection = collection,
          Rank_N = length(rank_list[[sensitivity]]), N_Sets_Evaluable = nrow(z),
          N_FDR_005 = if (nrow(z)) sum(z$FDR < 0.05) else 0L
        )
      }
    }
    cancer_result <- rbindlist(cancer_rows, use.names = TRUE, fill = TRUE)
    out_path <- file.path(tmp_dir, "individual",
      paste0(cancer, "_", tolower(edge_id), "_full_rank_gsea_v0.1_20260806.tsv.gz"))
    fwrite(cancer_result, out_path, sep = "\t", quote = TRUE, na = "NA", compress = "gzip")
    if (!file.exists(out_path) || file.info(out_path)$size <= 0) fail(paste(cancer, "逐癌GSEA输出为空"))
    out_manifest[[cancer]] <- data.table(
      Cancer = cancer,
      Output = file.path(final_dir, "individual", basename(out_path)),
      N_Rows = nrow(cancer_result)
    )
    all_results[[cancer]] <- cancer_result
    log_event("CHECKPOINT", "cancer edge GSEA complete", list(
      edge = edge_id, cancer = cancer, full_rank_n = length(full_rank),
      nondefinition_rank_n = length(rank_list$NON_DEFINITION_RANK),
      result_rows = nrow(cancer_result), significant_rows = sum(cancer_result$FDR < 0.05)))
  }
  gsea <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  if (!nrow(gsea)) fail(paste("所有癌种", edge_id, "GSEA均为空"))
  prefix <- tolower(edge_id)
  fwrite(gsea, file.path(tmp_dir, paste0(prefix, "_gsea_all_results.tsv.gz")), sep = "\t",
         quote = TRUE, na = "NA", compress = "gzip")
  fwrite(rbindlist(coverage), file.path(tmp_dir, paste0(prefix, "_gsea_coverage.tsv")), sep = "\t", quote = TRUE)

  # 跨癌只做等权方向与显著票数，不合并P值。
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
  setorder(vote, Sensitivity, Collection, -N_Target_High_Sig,
           -N_Reference_High_Sig, -Median_Abs_NES, Set_ID)
  fwrite(vote, file.path(tmp_dir, paste0(prefix, "_gsea_cross_cancer_votes.tsv")), sep = "\t", quote = TRUE)

  candidates <- vote[
    Sensitivity == "NON_DEFINITION_RANK" & Collection == "IMMUNE119_PRIMARY" &
      Support_Class %chin% c("STRONG_TARGET_HIGH", "STRONG_REFERENCE_HIGH")
  ]
  candidates[, `:=`(Manual_Redundancy_Review = "PENDING",
                    Fig3B_Eligibility = "PENDING_MANUAL_GATE")]
  fwrite(candidates, file.path(tmp_dir, paste0(prefix, "_gsea_nondefinition_immune_candidates.tsv")),
         sep = "\t", quote = TRUE)

  # FULL与去定义层按同一集合并排，便于识别定义驱动。
  comparison <- dcast(vote[Collection == "IMMUNE119_PRIMARY"],
    Collection + Set_ID + Display_Name + Functional_Domain ~ Sensitivity,
    value.var = c("N_Target_High_Sig", "N_Reference_High_Sig", "Support_Class"))
  fwrite(comparison, file.path(tmp_dir, paste0(prefix, "_gsea_definition_sensitivity_comparison.tsv")),
         sep = "\t", quote = TRUE)

  qa <- data.table(
    Check = c("REQUESTED_EDGE_ONLY", "EXPECTED_COLLECTIONS", "TWO_SENSITIVITIES",
              "IMMUNE119_EXACT", "LM12_EXACT", "NO_C7", "NO_META_FIELDS",
              "FDR_RANGE_VALID", "NONDEFINITION_RANK_EXCLUDES_DEFINITION",
              "EXPECTED_EVALUABLE_CANCERS", "OUTPUT_NONEMPTY"),
    Pass = c(
      all(gsea$Edge_ID == edge_id),
      setequal(unique(gsea$Collection), c("IMMUNE119_PRIMARY", "LM12_DEFINITION_CONTROL", "HALLMARK_CONTEXT")),
      setequal(unique(gsea$Sensitivity), c("FULL_RANK", "NON_DEFINITION_RANK")),
      uniqueN(immune_registry$Set_ID) == 119L,
      uniqueN(lm12_registry$Set_ID) == 12L,
      !any(grepl("C7", gsea$Collection)),
      !any(c("Meta_FDR", "I2", "Tau2", "LOO") %in% names(gsea)),
      all(is.finite(gsea$FDR) & gsea$FDR >= 0 & gsea$FDR <= 1),
      length(intersect(
        unique(unlist(strsplit(
          gsea[Sensitivity == "NON_DEFINITION_RANK", Leading_Edge],
          ";", fixed = TRUE
        ))),
        definition
      )) == 0L,
      uniqueN(gsea$Cancer) == edge_meta$Expected_Cancers,
      nrow(gsea) > 0L && nrow(vote) > 0L
    )
  )
  fwrite(qa, file.path(tmp_dir, "machine_qa.tsv"), sep = "\t", quote = TRUE)
  if (!all(qa$Pass)) fail(paste(edge_id, "full-rank GSEA机器QA失败"))

  summary <- vote[, .(
    N_Sets = .N,
    N_Recurrent_Target = sum(Recurrent_Target_High),
    N_Recurrent_Reference = sum(Recurrent_Reference_High),
    N_Strong_Target = sum(Strong_Target_High),
    N_Strong_Reference = sum(Strong_Reference_High)
  ), by = .(Edge_ID, Target_State, Reference_State, Sensitivity, Collection)]
  fwrite(summary, file.path(tmp_dir, paste0(prefix, "_gsea_summary.tsv")), sep = "\t", quote = TRUE)
  fwrite(rbindlist(out_manifest), file.path(tmp_dir, "output_manifest.tsv"), sep = "\t", quote = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(tmp_dir, "sessionInfo.txt"))
  report <- c(
    paste0("# ", edge_id, " full-rank GSEA执行报告 V0.1"), "", "日期：2026-08-06", "",
    "状态：`MACHINE_QA_PASS / MANUAL_GATE_PENDING / STAGE3_UNCHANGED`", "",
    "## 方法边界", "",
    sprintf("只复用既有%s逐癌moderated t；正NES表示%s/%s局部侧高，负NES表示%s/%s局部侧高。每癌种独立GSEA；跨癌仅做等权方向与显著票数，不运行meta、LOO或跨癌BH。", edge_id, edge_meta$Target, edge_meta$Target_State, edge_meta$Reference, edge_meta$Reference_State), "",
    "FULL_RANK为描述层；NON_DEFINITION_RANK同时从排名和集合中去除全部Acid/LM12 definition genes。只有去定义后的冻结免疫119层可以进入新增免疫模块人工审核。", "",
    "## 数量摘要", "",
    paste(apply(summary, 1L, function(x) sprintf(
      "%s × %s：%s个可评价集合；target侧复现/强支持%s/%s，reference侧复现/强支持%s/%s。",
      x[["Sensitivity"]], x[["Collection"]], x[["N_Sets"]],
      x[["N_Recurrent_Target"]], x[["N_Strong_Target"]],
      x[["N_Recurrent_Reference"]], x[["N_Strong_Reference"]]
    )), collapse = "\n\n"), "",
    "## 下一步", "",
    "与既有E2结果统一汇总后，人工审核三条边NON_DEFINITION_RANK强支持term的语义冗余、definition敏感性和Fig3B展示资格。"
  )
  writeLines(report, file.path(tmp_dir, "execution_report.md"))

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) fail(paste("无法提交", edge_id, "full-rank GSEA结果目录"))
  log_event("DONE", "single-edge full-rank GSEA complete", list(
    edge = edge_id,
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1),
    cancers = uniqueN(gsea$Cancer), gsea_rows = nrow(gsea),
    nondefinition_immune_strong_candidates = nrow(candidates),
    qa = sprintf("%d/%d", sum(qa$Pass), nrow(qa)), output = final_dir,
    manual_gate = "PENDING"))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)))
  quit(status = 1L)
})
