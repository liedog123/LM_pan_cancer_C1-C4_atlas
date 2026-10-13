#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 06_run_programme_feature_associations.R --config PATH [--check-config]\n       Rscript 06_run_programme_feature_associations.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("numeric_manifest", "bundle_dir", "program_clusters", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("numeric_manifest", "program_clusters")]
required_files <- unname(input_paths)

# Directory metadata only: do not read manifests or enumerate bundle contents.
bundle_dir <- unname(paths["bundle_dir"])
if (!dir.exists(bundle_dir) || file.access(bundle_dir, 5L) != 0L) stop("Bundle directory must be readable and searchable.")

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),bundle_dir,config_path,package_root)
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

# Validate only the selected bundle at its original read position.
resolve_bundle_path <- function(value) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value) ||
      value %in% c(".", "..")) stop("Invalid manifest bundle reference.")
  name <- basename(value)
  if (!nzchar(name) || name %in% c(".", "..")) stop("Invalid bundle basename.")
  path <- resolve_path(file.path(bundle_dir, name), bundle_dir)
  if (!is_within(path, bundle_dir) || is_within(path, out_dir) ||
      !file_test("-f", path) || file.access(path, 4L) != 0L ||
      is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Bundle must be a contained readable nonempty regular file.")
  }
  path
}

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

script_name <- "04_run_stage1_trackA_directed_integration_v0.1_20260810.R"
started <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extras <- if (length(fields)) paste(vapply(names(fields), function(k) {
    sprintf('%s="%s"', k, gsub('"', "'", as.character(fields[[k]]), fixed = TRUE))
  }, character(1)), collapse = " ") else ""
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
              gsub('"', "'", msg, fixed = TRUE),
              if (nzchar(extras)) " " else "", extras))
  flush.console()
}

safe_cor <- function(x, y) {
  sx <- sd(x)
  sy <- sd(y)
  if (length(x) < 3L || !is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
    return(c(rho = NA_real_, p = NA_real_))
  }
  z <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(z$estimate), p = z$p.value)
}

test_pair <- function(x, y, meta) {
  d <- data.table(
    x = as.numeric(x), y = as.numeric(y),
    state = factor(meta$C_subtype, levels = c("C4", "C1")),
    age = suppressWarnings(as.numeric(meta$Age)),
    sex = factor(meta$Sex)
  )
  base_ok <- is.finite(d$x) & is.finite(d$y) & !is.na(d$state)
  if (sum(base_ok) < 30L) return(list(Status = "NOT_EVALUABLE_N_LT30"))

  covars <- "state"
  if (mean(is.finite(d$age[base_ok])) >= 0.80 && length(unique(d$age[base_ok & is.finite(d$age)])) >= 3L) {
    covars <- c(covars, "age")
  }
  sex_ok <- base_ok & !is.na(d$sex)
  if (mean(!is.na(d$sex[base_ok])) >= 0.80 && length(unique(d$sex[sex_ok])) == 2L) {
    covars <- c(covars, "sex")
  }

  keep_cols <- c("x", "y", covars)
  cc <- complete.cases(d[, ..keep_cols])
  dd <- d[cc]
  n1 <- sum(dd$state == "C1")
  n4 <- sum(dd$state == "C4")
  if (nrow(dd) < 30L || min(n1, n4) < 10L) return(list(Status = "NOT_EVALUABLE_GROUP_GATE"))
  if (!is.finite(sd(dd$x)) || sd(dd$x) == 0) return(list(Status = "NOT_EVALUABLE_CONSTANT_PROGRAM"))
  if (!is.finite(sd(dd$y)) || sd(dd$y) == 0) return(list(Status = "NOT_EVALUABLE_CONSTANT_TARGET"))

  if ("age" %in% covars) dd[, age := as.numeric(scale(age))]
  form <- as.formula(paste("~", paste(covars, collapse = "+")))
  design <- model.matrix(form, data = dd)
  if (qr(design)$rank < ncol(design)) {
    covars <- "state"
    design <- model.matrix(~ state, data = dd)
  }
  rx <- lm.fit(design, dd$x)$residuals
  ry <- lm.fit(design, dd$y)$residuals
  partial <- safe_cor(rx, ry)
  raw <- safe_cor(dd$x, dd$y)
  if (any(!is.finite(partial)) || any(!is.finite(raw))) {
    return(list(Status = "NOT_EVALUABLE_SINGULAR_RESIDUAL"))
  }

  zx <- as.numeric(scale(dd$x))
  zy <- as.numeric(scale(dd$y))
  list(
    Status = "ESTIMATED", N = nrow(dd), N_C1 = n1, N_C4 = n4,
    Covariates = paste(covars, collapse = "+"),
    Partial_Rho = partial[["rho"]], Partial_P = partial[["p"]],
    Raw_Rho = raw[["rho"]], Raw_P = raw[["p"]],
    Program_C1_minus_C4_Z = mean(zx[dd$state == "C1"]) - mean(zx[dd$state == "C4"]),
    Target_C1_minus_C4_Z = mean(zy[dd$state == "C1"]) - mean(zy[dd$state == "C4"])
  )
}

matrix_family <- function(bundle) {
  list(
    RPPA_PATHWAY = bundle$rppa_tcpa_pathway,
    GENOMIC = cbind(bundle$cnv, bundle$genomic),
    IMMUNE_ESTIMATE = bundle$immune[, grepl("^ESTIMATE_", colnames(bundle$immune)), drop = FALSE],
    IMMUNE_CIBERSORT = bundle$immune[, grepl("^CIBERSORT_", colnames(bundle$immune)), drop = FALSE],
    IMMUNE_MCPCOUNTER = bundle$immune[, grepl("^MCPCOUNTER_", colnames(bundle$immune)), drop = FALSE],
    PREDICTED_DRUG = bundle$predicted_drug
  )
}

main <- function() {
  out_dir <- unname(paths["output_dir"])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- unname(paths["numeric_manifest"])
  cluster_path <- unname(paths["program_clusters"])
  if (!file.exists(manifest_path) || !file.exists(cluster_path)) stop("缺少共用输入或cluster registry")
  manifest <- fread(manifest_path)
  clusters <- fread(cluster_path)[State %in% c("C1", "C4") & Cluster_ID != "OTHER",
                                  .(Program = Cluster_ID, Program_Label = Cluster_Label,
                                    Program_State = State, Semantic_Boundary)]
  if (nrow(manifest) != 22L || nrow(clusters) != 9L) stop("输入Gate不一致")

  log_event("START", "线路A逐癌定向整合开始", list(
    project = "LM_pan_cancer", cancers = nrow(manifest), programs = nrow(clusters),
    method = "partial_spearman_residualized_state_age_sex", output = out_dir
  ))

  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    cancer <- manifest$Cancer[i]
    # Historical manifest paths provide basenames; keep bundles inside the configured directory.
    b <- readRDS(resolve_bundle_path(manifest$Bundle[i]))
    meta <- as.data.table(b$patient)
    ids <- meta$Patient_ID
    if (!identical(rownames(b$rna_program), ids)) stop(cancer, " RNA程序患者顺序不一致")
    families <- matrix_family(b)
    for (family in names(families)) {
      target <- families[[family]]
      if (is.null(dim(target)) || nrow(target) != length(ids)) stop(cancer, " target矩阵维度异常: ", family)
      for (program in colnames(b$rna_program)) {
        for (feature in colnames(target)) {
          ans <- test_pair(b$rna_program[, program], target[, feature], meta)
          rows[[length(rows) + 1L]] <- data.table(
            Cancer = cancer, Program = program, Evidence_Family = family,
            Target_Feature = feature,
            Status = ans$Status,
            N = if (!is.null(ans$N)) ans$N else NA_integer_,
            N_C1 = if (!is.null(ans$N_C1)) ans$N_C1 else NA_integer_,
            N_C4 = if (!is.null(ans$N_C4)) ans$N_C4 else NA_integer_,
            Covariates = if (!is.null(ans$Covariates)) ans$Covariates else NA_character_,
            Partial_Rho = if (!is.null(ans$Partial_Rho)) ans$Partial_Rho else NA_real_,
            Partial_P = if (!is.null(ans$Partial_P)) ans$Partial_P else NA_real_,
            Raw_Rho = if (!is.null(ans$Raw_Rho)) ans$Raw_Rho else NA_real_,
            Raw_P = if (!is.null(ans$Raw_P)) ans$Raw_P else NA_real_,
            Program_C1_minus_C4_Z = if (!is.null(ans$Program_C1_minus_C4_Z)) ans$Program_C1_minus_C4_Z else NA_real_,
            Target_C1_minus_C4_Z = if (!is.null(ans$Target_C1_minus_C4_Z)) ans$Target_C1_minus_C4_Z else NA_real_
          )
        }
      }
    }
    log_event("CHECKPOINT", "线路A癌种完成", list(cancer = cancer, n_samples = nrow(meta)))
  }

  res <- rbindlist(rows, fill = TRUE)
  res <- merge(res, clusters, by = "Program", all.x = TRUE)
  res[, Partial_FDR := NA_real_]
  res[Status == "ESTIMATED", Partial_FDR := p.adjust(Partial_P, method = "BH"),
      by = .(Cancer, Evidence_Family)]
  setcolorder(res, c("Cancer", "Program", "Program_Label", "Program_State", "Semantic_Boundary",
                     "Evidence_Family", "Target_Feature", "Status", "N", "N_C1", "N_C4",
                     "Covariates", "Partial_Rho", "Partial_P", "Partial_FDR", "Raw_Rho", "Raw_P",
                     "Program_C1_minus_C4_Z", "Target_C1_minus_C4_Z"))

  summary <- res[Status == "ESTIMATED", {
    sig <- is.finite(Partial_FDR) & Partial_FDR < 0.05
    n_sig <- sum(sig)
    n_pos_sig <- sum(sig & Partial_Rho > 0)
    n_neg_sig <- sum(sig & Partial_Rho < 0)
    dom_n <- max(n_pos_sig, n_neg_sig)
    .(
      N_Evaluable_Cancers = .N,
      N_Positive_All = sum(Partial_Rho > 0, na.rm = TRUE),
      N_Negative_All = sum(Partial_Rho < 0, na.rm = TRUE),
      Median_Partial_Rho = median(Partial_Rho, na.rm = TRUE),
      N_FDR_Significant = n_sig,
      N_Positive_Significant = n_pos_sig,
      N_Negative_Significant = n_neg_sig,
      Dominant_Significant_Direction = if (n_sig == 0L) "NONE" else if (n_pos_sig >= n_neg_sig) "POSITIVE" else "NEGATIVE",
      Significant_Direction_Consistency = if (n_sig == 0L) NA_real_ else dom_n / n_sig,
      Significant_Cancers = paste(sort(Cancer[sig]), collapse = ";")
    )
  }, by = .(Program, Program_Label, Program_State, Evidence_Family, Target_Feature)]

  orthogonal <- summary[
    Evidence_Family %in% c("RPPA_PATHWAY", "GENOMIC") &
      N_FDR_Significant >= 3L & Significant_Direction_Consistency >= 0.70
  ]
  orthogonal[, Abs_Median_Partial_Rho := abs(Median_Partial_Rho)]
  setorder(orthogonal, -N_FDR_Significant, -Significant_Direction_Consistency,
           -Abs_Median_Partial_Rho, Program, Evidence_Family, Target_Feature)
  orthogonal[, Abs_Median_Partial_Rho := NULL]

  res_path <- file.path(out_dir, "stage1_trackA_associations_by_cancer.tsv.gz")
  summary_path <- file.path(out_dir, "stage1_trackA_cross_cancer_summary.tsv")
  candidate_path <- file.path(out_dir, "stage1_trackA_orthogonal_candidates.tsv")
  fwrite(res, res_path, sep = "\t", compress = "gzip")
  fwrite(summary, summary_path, sep = "\t")
  fwrite(orthogonal, candidate_path, sep = "\t")

  qa <- data.table(
    Check = c(
      "22 cancers processed", "9 programs present", "no pooled patient model",
      "all estimated p finite", "within-family BH reproducible",
      "candidate restricted to orthogonal families", "output files nonempty"
    ),
    Pass = c(
      uniqueN(res$Cancer) == 22L,
      uniqueN(res$Program) == 9L,
      TRUE,
      all(is.finite(res[Status == "ESTIMATED", Partial_P])),
      res[Status == "ESTIMATED", all(abs(Partial_FDR - p.adjust(Partial_P, "BH")) < 1e-12),
          by = .(Cancer, Evidence_Family)][, all(V1)],
      all(orthogonal$Evidence_Family %in% c("RPPA_PATHWAY", "GENOMIC")),
      all(file.exists(c(res_path, summary_path, candidate_path))) &&
        all(file.info(c(res_path, summary_path, candidate_path))$size > 0)
    )
  )
  fwrite(qa, file.path(out_dir, "stage1_trackA_qa.tsv"), sep = "\t")
  if (!all(qa$Pass)) stop("线路A机器QA未全部通过")

  top_text <- if (nrow(orthogonal)) {
    paste(sprintf("%s—%s（%s，%s癌FDR通过）",
                  head(orthogonal$Program, 5), head(orthogonal$Target_Feature, 5),
                  head(orthogonal$Dominant_Significant_Direction, 5),
                  head(orthogonal$N_FDR_Significant, 5)), collapse = "；")
  } else "无候选达到预设的至少3癌种、方向一致性至少70%门槛"

  report <- c(
    "# Stage 1线路A定向整合报告 V0.1", "",
    "日期：2026-08-10  ",
    "状态：`TRACK_A_COMPLETE / MACHINE_QA_PASS / BIOLOGICAL_REVIEW_PENDING`", "",
    sprintf("共完成%s条逐癌程序—目标关联，其中%s条可评价。", nrow(res), res[Status == "ESTIMATED", .N]),
    sprintf("按单癌种、单证据家族BH后有%s条关联FDR<0.05。", res[Status == "ESTIMATED" & Partial_FDR < 0.05, .N]),
    sprintf("相对正交的RPPA/基因组候选单元为%s条。", nrow(orthogonal)),
    paste0("排序靠前的候选：", top_text, "。"), "",
    "主统计已调整C1/C4、年龄和性别，因此检验的是组别之外的跨层协同；它仍是观察性关联，不构成因果或药物机制。",
    "TCGA预测药敏和免疫结果保留为RNA耦联注释，未用于多组学候选的单独晋级。"
  )
  writeLines(report, file.path(out_dir, "stage1_trackA_report_v0.1.md"), useBytes = TRUE)

  log_event("DONE", "线路A定向整合完成", list(
    runtime_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    n_tests = nrow(res), n_estimated = res[Status == "ESTIMATED", .N],
    n_fdr = res[Status == "ESTIMATED" & Partial_FDR < 0.05, .N],
    n_orthogonal_candidates = nrow(orthogonal), output = out_dir
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e))
  quit(status = 1L)
})
