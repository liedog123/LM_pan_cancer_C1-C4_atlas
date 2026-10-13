#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 04_run_shared_axis_sensitivity.R --config PATH [--check-config]\n       Rscript 04_run_shared_axis_sensitivity.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("phase1_matches", "phase1_context", "drug_complete", "drug_final_grade", "shared_patient", "shared_cancer", "prediction_audit", "mrna", "rppa", "mutation", "tmb", "mantis", "msi", "cnv_burden", "cnv_highlevel", "prediction_dir", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("phase1_matches", "phase1_context", "drug_complete", "drug_final_grade", "shared_patient", "shared_cancer", "prediction_audit", "mrna", "rppa", "mutation", "tmb", "mantis", "msi", "cnv_burden", "cnv_highlevel")]
required_files <- unname(input_paths)

# 仅验证预测目录元信息，不读取audit内容或枚举预测文件。
prediction_dir <- unname(paths["prediction_dir"])
if (!dir.exists(prediction_dir) || file.access(prediction_dir, 5L) != 0L) stop("Prediction directory must be readable and searchable.")

if(any(!file_test("-f",required_files)) || any(file.access(required_files,4L)!=0L) ||
   any(is.na(file.info(required_files)$size)|file.info(required_files)$size<=0)) stop("Inputs must be readable nonempty regular files.")
required_files <- vapply(required_files,normalizePath,character(1L),mustWork=TRUE)
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[1L]),mustWork=TRUE)
package_root <- normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
protected <- c(dirname(required_files),prediction_dir,config_path,package_root)
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

options(stringsAsFactors = FALSE)
output_dir <- out_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  read.delim(
    con, sep = "\t", header = TRUE, quote = "\"", comment.char = "",
    check.names = FALSE, na.strings = c("NA", "")
  )
}

write_tsv <- function(x, name) {
  write.table(
    x, file.path(output_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

collapse_unique <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0L) "" else paste(x, collapse = ";")
}

cliffs_delta <- function(x, y) {
  n1 <- length(x)
  n0 <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u1 <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n0) - 1
}

signal_class <- function(direction, tested_n) {
  direction <- unique(direction[!is.na(direction) & nzchar(direction)])
  if (tested_n == 0L) return("NO_APPLICABLE_FEATURE")
  if (length(direction) == 0L) return("NO_SIGNIFICANT_SIGNAL")
  has_c1 <- "C1" %in% direction
  has_c4 <- "C4" %in% direction
  if (has_c1 && has_c4) return("MIXED_DIRECTION_SIGNAL")
  if (has_c1) return("C1_SIGNAL_ONLY")
  if (has_c4) return("C4_SIGNAL_ONLY")
  "NO_SIGNIFICANT_SIGNAL"
}

state_assessment <- function(signal, sensitive_state) {
  if (signal %in% c("NO_APPLICABLE_FEATURE", "NO_SIGNIFICANT_SIGNAL")) {
    return(signal)
  }
  if (signal == "MIXED_DIRECTION_SIGNAL") return("MIXED_STATE_ASSOCIATION")
  signal_state <- sub("_SIGNAL_ONLY$", "", signal)
  if (signal_state == sensitive_state) {
    "SAME_STATE_COMPATIBLE_ASSOCIATION"
  } else {
    "OPPOSITE_STATE_ASSOCIATION"
  }
}

manifest <- data.frame(
  input_id = names(input_paths),
  path = unname(input_paths),
  exists = file.exists(input_paths),
  nonempty = file.exists(input_paths) & file.info(input_paths)$size > 0L,
  size_bytes = ifelse(file.exists(input_paths), file.info(input_paths)$size, NA_real_),
  md5 = ifelse(file.exists(input_paths), unname(tools::md5sum(input_paths)), NA_character_)
)
write_tsv(manifest, "input_manifest.tsv")
if (!all(manifest$nonempty)) stop("Required Phase 1 or Stage 3-5 input is missing/empty.")

d <- lapply(input_paths, read_tsv)
strict_drugs <- d$drug_final_grade$drug[d$drug_final_grade$evidence_grade == "STRICT_SUPPORTED"]
if (!setequal(strict_drugs, c("Vinorelbine", "Linsitinib", "Palbociclib"))) {
  stop("Strict drug freeze differs from the Phase 2 contract.")
}

raw_drug <- d$drug_complete[
  d$drug_complete$drug %in% strict_drugs & d$drug_complete$status == "ok", ,
  drop = FALSE
]
formal_cancers <- sort(unique(raw_drug$Cancer[raw_drug$significant %in% TRUE]))
raw_drug <- raw_drug[raw_drug$Cancer %in% formal_cancers, , drop = FALSE]
if (length(formal_cancers) != 18L || nrow(raw_drug) != 54L) {
  stop("Expected 18 formal cancers and 54 Cancer x Drug units.")
}

# 每个药物使用排除自身后的285药物PC1，避免当前药物直接进入共享轴。
residual_rows <- list()
rr <- 0L
for (cancer in formal_cancers) {
  audit_row <- d$prediction_audit[d$prediction_audit$Cancer == cancer, , drop = FALSE]
  if (nrow(audit_row) != 1L) stop("Prediction source unresolved for ", cancer)
  # 内嵌历史路径只提供文件名；在配置预测目录内解析并保护输出边界。
  source_name <- basename(as.character(audit_row$source_file))
  if (length(source_name) != 1L || is.na(source_name) || !nzchar(source_name) ||
      source_name %in% c(".", "..") || grepl("[/\\\\]", source_name)) {
    stop("Prediction basename unresolved for ", cancer)
  }
  prediction_path <- resolve_path(file.path(prediction_dir, source_name), prediction_dir)
  if (!is_within(prediction_path, prediction_dir) || is_within(prediction_path, out_dir) ||
      !file_test("-f", prediction_path) || file.access(prediction_path, 4L) != 0L ||
      is.na(file.info(prediction_path)$size) || file.info(prediction_path)$size <= 0L) {
    stop("Prediction source must be a protected readable nonempty file for ", cancer)
  }
  pred <- read.csv(
    prediction_path, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sample_col <- names(pred)[1L]
  pred$Sample15 <- substr(as.character(pred[[sample_col]]), 1L, 15L)
  shared <- d$shared_patient[d$shared_patient$Cancer == cancer, , drop = FALSE]
  idx <- match(shared$Sample15, pred$Sample15)
  if (anyNA(idx)) stop("Shared-score patients missing from prediction matrix: ", cancer)
  pred <- pred[idx, , drop = FALSE]
  if (!identical(pred$Sample15, shared$Sample15)) stop("Sample order mismatch: ", cancer)

  drug_cols <- setdiff(names(pred), c(sample_col, "Sample15"))
  if (!all(strict_drugs %in% drug_cols)) stop("Strict drug column missing: ", cancer)

  cancer_rows <- list()
  for (drug in strict_drugs) {
    other_cols <- setdiff(drug_cols, drug)
    other_mat <- as.matrix(data.frame(lapply(pred[other_cols], as.numeric), check.names = FALSE))
    keep <- apply(other_mat, 2L, function(x) all(is.finite(x)) && stats::sd(x) > 0)
    other_mat <- other_mat[, keep, drop = FALSE]
    if (ncol(other_mat) < 280L) stop("Too few nonconstant leave-one-drug-out columns: ", cancer)
    pca <- prcomp(other_mat, center = TRUE, scale. = TRUE, rank. = 1L)
    pc1 <- pca$x[, 1L]
    pc1_var <- pca$sdev[1L]^2 / sum(pca$sdev^2)

    y <- as.numeric(pred[[drug]])
    fit <- lm(y ~ pc1)
    resid_y <- residuals(fit)
    is_c1 <- shared$C_subtype == "C1"
    x <- resid_y[is_c1]
    y4 <- resid_y[!is_c1]
    wt <- suppressWarnings(wilcox.test(x, y4, exact = FALSE))
    delta <- cliffs_delta(x, y4)
    median_diff <- median(x) - median(y4)
    raw_row <- raw_drug[raw_drug$Cancer == cancer & raw_drug$drug == drug, , drop = FALSE]

    rr <- rr + 1L
    cancer_rows[[drug]] <- data.frame(
      drug = drug,
      Cancer = cancer,
      n_C1 = sum(is_c1),
      n_C4 = sum(!is_c1),
      n_other_drugs_for_PC1 = ncol(other_mat),
      loo_PC1_variance_fraction = pc1_var,
      abs_spearman_drug_vs_loo_PC1 = abs(suppressWarnings(cor(
        as.numeric(pred[[drug]]), pc1, method = "spearman"
      ))),
      raw_cliffs_delta_C1_minus_C4 = raw_row$cliffs_delta_C1_minus_C4,
      raw_more_sensitive_group = raw_row$more_sensitive_group,
      raw_FDR = raw_row$FDR,
      residual_median_diff_C1_minus_C4 = median_diff,
      residual_cliffs_delta_C1_minus_C4 = delta,
      residual_p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
  cancer_df <- do.call(rbind, cancer_rows)
  cancer_df$residual_FDR_within_cancer_3_drugs <- p.adjust(
    cancer_df$residual_p_value, method = "BH"
  )
  cancer_df$residual_more_sensitive_group <- ifelse(
    cancer_df$residual_median_diff_C1_minus_C4 < 0, "C1",
    ifelse(cancer_df$residual_median_diff_C1_minus_C4 > 0, "C4", "tie")
  )
  cancer_df$residual_significant <- (
    cancer_df$residual_FDR_within_cancer_3_drugs < 0.05 &
      abs(cancer_df$residual_cliffs_delta_C1_minus_C4) >= 0.20
  )
  cancer_df$direction_preserved_after_residual <- (
    cancer_df$residual_more_sensitive_group == cancer_df$raw_more_sensitive_group
  )
  residual_rows[[cancer]] <- cancer_df
}
residual <- do.call(rbind, residual_rows)
rownames(residual) <- NULL
write_tsv(residual, "strict_drug_shared_axis_residual_by_cancer.tsv")

# Phase 1冻结基因/feature是唯一允许的逐癌查询集合。
phase1 <- d$phase1_matches[d$phase1_matches$drug %in% strict_drugs, , drop = FALSE]
phase1_gene <- unique(phase1[phase1$gene != "MODULE", c("drug", "gene", "target_role")])

evidence <- data.frame()
append_evidence <- function(x) {
  evidence <<- rbind(evidence, x)
}

# mRNA：只作与预测模型相关的表达语境。
for (i in seq_len(nrow(phase1_gene))) {
  z <- d$mrna[
    d$mrna$Gene == phase1_gene$gene[i] &
      d$mrna$Design_ID == "C1_vs_C4" &
      d$mrna$Cancer %in% formal_cancers, ,
    drop = FALSE
  ]
  if (nrow(z) == 0L) next
  append_evidence(data.frame(
    drug = phase1_gene$drug[i], Cancer = z$Cancer,
    gene = z$Gene, feature = z$Feature_ID,
    target_role = phase1_gene$target_role[i],
    source_domain = "mRNA",
    evidence_role = "MODEL_INPUT_RELATED_EXPRESSION_CONTEXT",
    tested = z$Status == "ESTIMABLE",
    significant = z$Significant_In_Direction %in% TRUE,
    signal_state = ifelse(z$Estimate > 0, "C1", ifelse(z$Estimate < 0, "C4", NA)),
    effect = z$Estimate, within_cancer_FDR = z$FDR,
    detail = paste0("direction=", z$Direction),
    stringsAsFactors = FALSE
  ))
}

# RPPA：只使用Phase 1已匹配的具体feature。
rppa_phase1 <- unique(phase1[phase1$source_domain == "RPPA", c("drug", "gene", "target_role", "source_detail")])
if (nrow(rppa_phase1) > 0L) {
  rppa_phase1$feature <- sub(";.*$", "", rppa_phase1$source_detail)
  for (i in seq_len(nrow(rppa_phase1))) {
    z <- d$rppa[
      d$rppa$Feature_ID == rppa_phase1$feature[i] &
        d$rppa$Cancer %in% formal_cancers, ,
      drop = FALSE
    ]
    if (nrow(z) == 0L) next
    append_evidence(data.frame(
      drug = rppa_phase1$drug[i], Cancer = z$Cancer,
      gene = rppa_phase1$gene[i], feature = z$Feature_ID,
      target_role = rppa_phase1$target_role[i],
      source_domain = "RPPA",
      evidence_role = "ORTHOGONAL_PROTEIN_CONTEXT",
      tested = z$Status == "ESTIMABLE",
      significant = z$Single_Cancer_Direction %in% c("C1_HIGH", "C4_HIGH"),
      signal_state = ifelse(
        z$Single_Cancer_Direction == "C1_HIGH", "C1",
        ifelse(z$Single_Cancer_Direction == "C4_HIGH", "C4", NA)
      ),
      effect = z$Estimate, within_cancer_FDR = z$FDR,
      detail = paste0(z$Feature_Class, ";", z$Single_Cancer_Direction),
      stringsAsFactors = FALSE
    ))
  }
}

# Mutation和high-level CNV只查询Phase 1中已有匹配的基因。
mutation_phase1 <- unique(phase1[phase1$source_domain == "mutation", c("drug", "gene", "target_role")])
if (nrow(mutation_phase1) > 0L) {
  for (i in seq_len(nrow(mutation_phase1))) {
    z <- d$mutation[
      d$mutation$Gene == mutation_phase1$gene[i] &
        d$mutation$Cancer %in% formal_cancers, ,
      drop = FALSE
    ]
    if (nrow(z) == 0L) next
    append_evidence(data.frame(
      drug = mutation_phase1$drug[i], Cancer = z$Cancer,
      gene = z$Gene, feature = z$Gene,
      target_role = mutation_phase1$target_role[i],
      source_domain = "mutation",
      evidence_role = "GENOMIC_BACKGROUND_CONTEXT",
      tested = TRUE,
      significant = z$single_cancer_significant %in% TRUE,
      signal_state = ifelse(z$direction == "C1_higher", "C1", "C4"),
      effect = z$frequency_difference_C1_minus_C4,
      within_cancer_FDR = z$FDR_within_cancer,
      detail = "nonsilent_callable_mutation",
      stringsAsFactors = FALSE
    ))
  }
}

cnv_phase1 <- unique(phase1[phase1$source_domain == "CNV_highlevel", c("drug", "gene", "target_role")])
if (nrow(cnv_phase1) > 0L) {
  for (i in seq_len(nrow(cnv_phase1))) {
    z <- d$cnv_highlevel[
      d$cnv_highlevel$Gene == cnv_phase1$gene[i] &
        d$cnv_highlevel$Cancer %in% formal_cancers, ,
      drop = FALSE
    ]
    if (nrow(z) == 0L) next
    append_evidence(data.frame(
      drug = cnv_phase1$drug[i], Cancer = z$Cancer,
      gene = z$Gene, feature = paste(z$Gene, z$event_direction, sep = "::"),
      target_role = cnv_phase1$target_role[i],
      source_domain = "CNV_highlevel",
      evidence_role = "GENOMIC_BACKGROUND_CONTEXT",
      tested = TRUE,
      significant = z$significant %in% TRUE,
      signal_state = ifelse(z$direction == "C1_higher", "C1", "C4"),
      effect = z$effect_c1_minus_c4,
      within_cancer_FDR = z$q_value,
      detail = z$event_direction,
      stringsAsFactors = FALSE
    ))
  }
}

evidence <- unique(evidence)
write_tsv(evidence, "strict_drug_gene_feature_evidence_long.tsv")

# 全局基因组背景只作Cancer-level context，不作药物特异支持。
get_context <- function(x, prefix, effect_col = "effect_c1_minus_c4") {
  z <- x[x$Cancer %in% formal_cancers, , drop = FALSE]
  data.frame(
    Cancer = z$Cancer,
    value = z[[effect_col]],
    significant = z$significant,
    direction = z$direction,
    stringsAsFactors = FALSE
  )
}
context <- data.frame(Cancer = formal_cancers)
for (spec in list(
  list("TMB", d$tmb),
  list("MANTIS", d$mantis),
  list("MSI", d$msi),
  list("CNV_burden", d$cnv_burden[d$cnv_burden$formal_primary %in% TRUE, , drop = FALSE])
)) {
  z <- get_context(spec[[2]], spec[[1]])
  names(z)[-1L] <- paste0(tolower(spec[[1]]), c("_effect", "_significant", "_direction"))
  context <- merge(context, z, by = "Cancer", all.x = TRUE, sort = FALSE)
}
context <- context[match(formal_cancers, context$Cancer), , drop = FALSE]
write_tsv(context, "genomic_context_by_cancer.tsv")

# 汇总每个Cancer × Drug的表达语境和相对正交状态背景。
cross_rows <- list()
counter <- 0L
for (drug in strict_drugs) {
  for (cancer in formal_cancers) {
    counter <- counter + 1L
    z <- evidence[evidence$drug == drug & evidence$Cancer == cancer, , drop = FALSE]
    summarize_domain <- function(domain) {
      q <- z[z$source_domain == domain, , drop = FALSE]
      tested_n <- sum(q$tested %in% TRUE)
      sig <- q[q$tested %in% TRUE & q$significant %in% TRUE, , drop = FALSE]
      list(
        class = signal_class(sig$signal_state, tested_n),
        tested_n = tested_n,
        sig_n = nrow(sig),
        genes = collapse_unique(sig$gene),
        features = collapse_unique(sig$feature)
      )
    }
    m <- summarize_domain("mRNA")
    p <- summarize_domain("RPPA")
    u <- summarize_domain("mutation")
    c <- summarize_domain("CNV_highlevel")
    orth <- z[
      z$source_domain %in% c("RPPA", "mutation", "CNV_highlevel") &
        z$tested %in% TRUE & z$significant %in% TRUE, ,
      drop = FALSE
    ]
    orth_tested <- sum(
      z$source_domain %in% c("RPPA", "mutation", "CNV_highlevel") &
        z$tested %in% TRUE
    )
    orth_class <- signal_class(orth$signal_state, orth_tested)
    rr0 <- residual[residual$drug == drug & residual$Cancer == cancer, , drop = FALSE]
    cross_rows[[counter]] <- data.frame(
      drug = drug, Cancer = cancer,
      raw_sensitive_state = rr0$raw_more_sensitive_group,
      raw_cliffs_delta = rr0$raw_cliffs_delta_C1_minus_C4,
      residual_sensitive_state = rr0$residual_more_sensitive_group,
      residual_cliffs_delta = rr0$residual_cliffs_delta_C1_minus_C4,
      residual_FDR = rr0$residual_FDR_within_cancer_3_drugs,
      residual_significant = rr0$residual_significant,
      residual_direction_preserved = rr0$direction_preserved_after_residual,
      mRNA_signal_class = m$class,
      mRNA_n_tested = m$tested_n,
      mRNA_n_significant = m$sig_n,
      mRNA_significant_genes = m$genes,
      mRNA_state_assessment = state_assessment(m$class, rr0$raw_more_sensitive_group),
      RPPA_signal_class = p$class,
      RPPA_n_tested = p$tested_n,
      RPPA_n_significant = p$sig_n,
      RPPA_significant_features = p$features,
      mutation_signal_class = u$class,
      CNV_highlevel_signal_class = c$class,
      orthogonal_signal_class = orth_class,
      orthogonal_significant_genes = collapse_unique(orth$gene),
      orthogonal_state_assessment = state_assessment(orth_class, rr0$raw_more_sensitive_group),
      stringsAsFactors = FALSE
    )
  }
}
cross <- do.call(rbind, cross_rows)
cross <- merge(cross, context, by = "Cancer", all.x = TRUE, sort = FALSE)
cross <- merge(
  cross,
  d$shared_cancer[, c(
    "Cancer", "PC1_variance_fraction", "significant_fraction",
    "shared_axis_evidence"
  )],
  by = "Cancer", all.x = TRUE, sort = FALSE
)
cross <- cross[order(match(cross$drug, strict_drugs), match(cross$Cancer, formal_cancers)), ]
write_tsv(cross, "strict_drug_crosslayer_by_cancer.tsv")

# 仅统计癌种数，不计算跨癌P值。
summary_rows <- lapply(strict_drugs, function(drug) {
  z <- cross[cross$drug == drug, , drop = FALSE]
  n_resid <- sum(z$residual_significant %in% TRUE)
  n_preserved <- sum(z$residual_significant %in% TRUE & z$residual_direction_preserved %in% TRUE)
  preserved_fraction <- if (n_resid > 0L) n_preserved / n_resid else NA_real_
  shared_class <- if (n_resid >= 9L && !is.na(preserved_fraction) && preserved_fraction >= 0.80) {
    "PERSISTENT_BEYOND_SHARED_AXIS"
  } else if (n_resid >= 4L) {
    "PARTIAL_PERSISTENCE_CONTEXT_DEPENDENT"
  } else {
    "LARGELY_ATTENUATED_BY_SHARED_AXIS"
  }
  next_step <- switch(
    shared_class,
    PERSISTENT_BEYOND_SHARED_AXIS = "RETAIN_FOR_TARGETED_VALIDATION",
    PARTIAL_PERSISTENCE_CONTEXT_DEPENDENT = "RETAIN_CONTEXT_DEPENDENT",
    LARGELY_ATTENUATED_BY_SHARED_AXIS = "DOWNGRADE_NO_MECHANISTIC_EXTENSION"
  )
  data.frame(
    drug = drug,
    n_formal_cancers = nrow(z),
    n_raw_significant = 18L,
    n_residual_significant = n_resid,
    n_residual_C1_sensitive = sum(z$residual_significant & z$residual_sensitive_state == "C1"),
    n_residual_C4_sensitive = sum(z$residual_significant & z$residual_sensitive_state == "C4"),
    n_residual_direction_preserved = n_preserved,
    residual_direction_preserved_fraction = preserved_fraction,
    median_abs_raw_cliffs_delta = median(abs(z$raw_cliffs_delta)),
    median_abs_residual_cliffs_delta = median(abs(z$residual_cliffs_delta)),
    n_mRNA_same_state = sum(z$mRNA_state_assessment == "SAME_STATE_COMPATIBLE_ASSOCIATION"),
    n_mRNA_opposite_state = sum(z$mRNA_state_assessment == "OPPOSITE_STATE_ASSOCIATION"),
    n_orthogonal_same_state = sum(z$orthogonal_state_assessment == "SAME_STATE_COMPATIBLE_ASSOCIATION"),
    n_orthogonal_opposite_state = sum(z$orthogonal_state_assessment == "OPPOSITE_STATE_ASSOCIATION"),
    n_orthogonal_mixed = sum(z$orthogonal_state_assessment == "MIXED_STATE_ASSOCIATION"),
    n_orthogonal_no_signal = sum(z$orthogonal_state_assessment %in% c(
      "NO_SIGNIFICANT_SIGNAL", "NO_APPLICABLE_FEATURE"
    )),
    shared_axis_diagnostic_class = shared_class,
    phase2_next_step = next_step,
    interpretation_boundary = paste0(
      "Residual diagnostic and same-state context only; mRNA is model-input-related; ",
      "RPPA/genomics do not validate drug response."
    ),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
write_tsv(summary, "strict_drug_phase2_summary.tsv")

qa <- data.frame(
  check_id = c(
    "QA01_ALL_INPUTS_EXIST_NONEMPTY",
    "QA02_EXACTLY_THREE_FROZEN_STRICT_DRUGS",
    "QA03_EXACTLY_18_FORMAL_CANCERS",
    "QA04_RAW_UNITS_EQUAL_54",
    "QA05_RESIDUAL_UNITS_EQUAL_54",
    "QA06_THREE_DRUGS_PER_CANCER",
    "QA07_LOO_PC1_EXCLUDES_CURRENT_DRUG",
    "QA08_RESIDUAL_STATISTICS_FINITE",
    "QA09_RESIDUAL_FDR_IN_RANGE",
    "QA10_PHASE1_GENE_RESTRICTION",
    "QA11_CROSSLAYER_UNITS_EQUAL_54",
    "QA12_SUMMARY_HAS_THREE_DRUGS",
    "QA13_ALL_SHARED_AXIS_EVIDENCE_RETAINED",
    "QA14_NO_POOLED_OR_CROSS_CANCER_P_VALUE",
    "QA15_MRNA_STRICTLY_C1_VS_C4"
  ),
  pass = c(
    all(manifest$nonempty),
    length(strict_drugs) == 3L,
    length(formal_cancers) == 18L,
    nrow(raw_drug) == 54L,
    nrow(residual) == 54L,
    all(table(residual$Cancer) == 3L),
    all(residual$n_other_drugs_for_PC1 >= 285L),
    all(is.finite(residual$residual_cliffs_delta_C1_minus_C4)) &
      all(is.finite(residual$residual_p_value)),
    all(residual$residual_FDR_within_cancer_3_drugs >= 0 &
          residual$residual_FDR_within_cancer_3_drugs <= 1),
    all(unique(evidence[, c("drug", "gene")])$gene %in% phase1_gene$gene),
    nrow(cross) == 54L,
    nrow(summary) == 3L & setequal(summary$drug, strict_drugs),
    all(cross$shared_axis_evidence %in% TRUE),
    !any(grepl("pooled|meta_p|cross_cancer_p", names(summary), ignore.case = TRUE)),
    all(evidence$source_domain != "mRNA" |
          grepl("^direction=(UP|DOWN)$", evidence$detail)) &
      !anyDuplicated(evidence[
        evidence$source_domain == "mRNA",
        c("drug", "Cancer", "gene", "feature")
      ])
  ),
  detail = c(
    paste0(sum(manifest$nonempty), "/", nrow(manifest)),
    collapse_unique(strict_drugs),
    paste0(length(formal_cancers), " cancers"),
    paste0(nrow(raw_drug), " units"),
    paste0(nrow(residual), " units"),
    collapse_unique(paste(names(table(residual$Cancer)), table(residual$Cancer), sep = ":")),
    paste0("min_other_drugs=", min(residual$n_other_drugs_for_PC1)),
    "all residual effects and P values finite",
    paste0("range=", signif(min(residual$residual_FDR_within_cancer_3_drugs), 4),
           "-", signif(max(residual$residual_FDR_within_cancer_3_drugs), 4)),
    paste0(nrow(phase1_gene), " frozen drug-gene rows"),
    paste0(nrow(cross), " units"),
    collapse_unique(summary$drug),
    paste0(sum(cross$shared_axis_evidence %in% TRUE), "/", nrow(cross)),
    "cross-cancer integration uses counts only",
    "mRNA source filtered to Design_ID=C1_vs_C4 with unique drug-Cancer-gene-feature"
  )
)
write_tsv(qa, "machine_qa.tsv")

status <- data.frame(
  phase = "STAGE6_PHASE2_STRICT_DRUG_BY_CANCER",
  status = if (all(qa$pass)) "PASS_PHASE2_STRICT_DRUG_DIRECTION_COMPLETE" else "HOLD_PHASE2_QA_FAILED",
  gate6 = "OPEN_NOT_EVALUATED",
  final_figure = "NOT_RUN",
  ppi = "NOT_RUN",
  molecular_docking = "NOT_RUN",
  single_cell_virtual_perturbation = "NOT_RUN"
)
write_tsv(status, "stage6_phase2_status.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

if (!all(qa$pass)) stop("Stage 6 Phase 2 QA failed.")
message(
  "Stage 6 Phase 2 v0.2 completed: 54 residual units; ",
  paste(summary$drug, summary$n_residual_significant, sep = "=", collapse = ", "),
  " residual-significant cancers."
)
