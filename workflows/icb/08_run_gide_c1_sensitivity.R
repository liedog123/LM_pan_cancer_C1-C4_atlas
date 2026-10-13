#!/usr/bin/env Rscript
# Base-R metadata checks precede dependency loading and study reads.
usage <- function() cat("Usage: Rscript 08_run_gide_c1_sensitivity.R --config PATH [--check-config]\n       Rscript 08_run_gide_c1_sensitivity.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("patient_scores","frozen_program_members","lm12_members","output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
cfg_out_dir <- unname(paths["output_dir"])
input_paths <- paths[c("patient_scores","frozen_program_members","lm12_members")]
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
# 普通运行才创建独立空输出；检查模式已退出。
dir.create(unname(cfg_paths["output_dir"]),recursive=TRUE,showWarnings=FALSE)

# ---- Retained complete source body; documented path/message substitutions only ----

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# Private project root replaced by explicit DCF inputs.
# Module paths are supplied by the checked DCF.
src <- unname(cfg_paths["patient_scores"])
out <- unname(cfg_paths["output_dir"])
res <- file.path(out, "results")
fig <- file.path(out, "figures")
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

log_event <- function(level, msg) {
  cat(sprintf("[%s] %s script=08_run_gide_c1_sensitivity.R msg=\"%s\"\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg))
  flush.console()
}

log_event("START", paste0(
  "post-hoc sensitivity in Gide C1; input=", src,
  "; target=regimen distribution and Cytotoxic-LM12 axis relation",
  "; Firth profile 95% CI; Spearman patient bootstrap B=5000 seed=20260923",
  "; outputs=", out))

d <- fread(src)
x <- d[cohort == "Gide" & C_state == "C1" & !is.na(strict_response)]
if (nrow(x) != 25L || uniqueN(x$patient_id) != 25L || uniqueN(x$sample_id) != 25L ||
    sum(x$strict_response == "R") != 13L || sum(x$strict_response == "NR") != 12L) {
  stop("Gide C1 differs from the fixed 25-patient input summary; review inclusion before analysis")
}
if (anyNA(x$regimen) || any(!x$strict_response %in% c("R", "NR")) ||
    any(!is.finite(x$Cytotoxic)) || any(!is.finite(x$LM12_inflammatory_axis))) {
  stop("Gide C1 has missing or invalid values in required outcome, regimen, or score fields")
}

# 沿原始治疗标签归并为单药/联合；保留原始 regimen 供核查。
regimen_lower <- tolower(trimws(x$regimen))
is_combo <- grepl("ctla|ipilimumab|combination|combined", regimen_lower)
is_mono <- grepl("pd.?1|nivolumab|pembrolizumab", regimen_lower) & !is_combo
if (any(!is_combo & !is_mono)) {
  stop(paste0("Unresolved Gide C1 regimen label(s): ",
              paste(unique(x$regimen[!is_combo & !is_mono]), collapse = " | ")))
}
x[, regimen_group := ifelse(is_combo, "Combination", "Monotherapy")]
x[, regimen_group := factor(regimen_group, levels = c("Monotherapy", "Combination"))]
x[, response_label := factor(strict_response, levels = c("NR", "R"),
                              labels = c("SD/PD", "CR/PR"))]

# 列出每种方案内响应率，并同时显示每种方案占响应者/非响应者的比例。
regimen_summary <- x[, .(
  n = .N,
  n_R = sum(strict_response == "R"),
  n_NR = sum(strict_response == "NR"),
  regimen_source = paste(sort(unique(regimen)), collapse = " | ")
), by = .(regimen_group)]
regimen_summary[, response_rate := n_R / n]
regimen_summary[, `:=`(
  response_ci_low = mapply(function(r, n) binom.test(r, n)$conf.int[1], n_R, n),
  response_ci_high = mapply(function(r, n) binom.test(r, n)$conf.int[2], n_R, n),
  share_among_R = n_R / sum(x$strict_response == "R"),
  share_among_NR = n_NR / sum(x$strict_response == "NR")
)]
setorder(regimen_summary, regimen_group)
fwrite(regimen_summary, file.path(res, "Gide_C1_regimen_response_distribution.tsv"),
       sep = "\t", na = "NA")

# 只作小样本描述：展示各治疗/响应格子的冻结 Cytotoxic 分数分布。
cell_grid <- CJ(regimen_group = levels(x$regimen_group),
                strict_response = c("R", "NR"), unique = TRUE)
cyt_summary <- x[, .(
  n = .N,
  mean = mean(Cytotoxic),
  median = median(Cytotoxic),
  min = min(Cytotoxic),
  max = max(Cytotoxic),
  q1 = unname(quantile(Cytotoxic, .25)),
  q3 = unname(quantile(Cytotoxic, .75))
), by = .(regimen_group = as.character(regimen_group), strict_response)]
cyt_summary <- merge(cell_grid, cyt_summary,
                     by = c("regimen_group", "strict_response"), all.x = TRUE, sort = FALSE)
cyt_summary[is.na(n), n := 0L]
setorder(cyt_summary, regimen_group, strict_response)
fwrite(cyt_summary, file.path(res, "Gide_C1_cytotoxic_by_regimen_response.tsv"),
       sep = "\t", na = "NA")

# 在同一完整病例中比较原始与 regimen 调整的 Firth logistic OR/区间，不作第二套显著性筛选。
cc <- x[complete.cases(x[, .(strict_response, regimen_group, Cytotoxic)])]
cc[, response_binary := as.integer(strict_response == "R")]
cc[, Cytotoxic_cohort_SD := Cytotoxic]

firth_profile_ci <- function(X, y, index, profile = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(y)
  if (qr(X)$rank < ncol(X)) stop("Firth design matrix is rank deficient")
  penalized_loglik <- function(beta) {
    eta <- drop(X %*% beta)
    prob <- plogis(eta)
    weight <- pmax(prob * (1 - prob), 1e-15)
    info <- crossprod(X * sqrt(weight))
    det_info <- determinant(info, logarithm = TRUE)
    if (det_info$sign <= 0 || !is.finite(det_info$modulus)) return(-1e50)
    loglik <- sum(y * eta - ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta))))
    loglik + 0.5 * as.numeric(det_info$modulus)
  }
  full <- optim(rep(0, ncol(X)), function(beta) -penalized_loglik(beta),
                method = "BFGS", control = list(maxit = 1500, reltol = 1e-10))
  if (full$convergence != 0L) stop("Firth full model optimization failed")
  beta_hat <- full$par
  if (!profile) return(c(OR = exp(beta_hat[index]), logOR = beta_hat[index]))
  ll_top <- penalized_loglik(beta_hat)
  remaining <- setdiff(seq_along(beta_hat), index)
  profile_ll <- function(value) {
    beta <- beta_hat
    beta[index] <- value
    if (!length(remaining)) return(penalized_loglik(beta))
    opt <- optim(beta_hat[remaining], function(rest) {
      candidate <- beta
      candidate[remaining] <- rest
      -penalized_loglik(candidate)
    }, method = "BFGS", control = list(maxit = 1500, reltol = 1e-9))
    if (opt$convergence != 0L) stop("Firth profile optimization failed")
    -opt$value
  }
  target <- ll_top - qchisq(.95, df = 1) / 2
  find_bound <- function(side) {
    step <- .5
    bound <- beta_hat[index] + side * step
    while (profile_ll(bound) > target && step < 32) {
      step <- step * 2
      bound <- beta_hat[index] + side * step
    }
    if (profile_ll(bound) > target) return(if (side < 0) -Inf else Inf)
    interval <- if (side < 0) c(bound, beta_hat[index]) else c(beta_hat[index], bound)
    uniroot(function(value) profile_ll(value) - target, interval, tol = 1e-5)$root
  }
  c(OR = exp(beta_hat[index]),
    CI_low = exp(find_bound(-1)),
    CI_high = exp(find_bound(1)))
}

fit_one <- function(dat, adjusted) {
  if (adjusted && uniqueN(dat$regimen_group) > 1L) {
    design_data <- data.frame(Cytotoxic_cohort_SD = dat$Cytotoxic_cohort_SD,
                              regimen_group = droplevels(dat$regimen_group))
    X <- model.matrix(~ Cytotoxic_cohort_SD + regimen_group, design_data)
    term <- which(colnames(X) == "Cytotoxic_cohort_SD")
    adjustment <- "regimen group"
  } else {
    design_data <- data.frame(Cytotoxic_cohort_SD = dat$Cytotoxic_cohort_SD)
    X <- model.matrix(~ Cytotoxic_cohort_SD, design_data)
    term <- which(colnames(X) == "Cytotoxic_cohort_SD")
    adjustment <- if (adjusted) "regimen has one level; omitted" else "none"
  }
  est <- firth_profile_ci(X, dat$response_binary, term)
  c1_sd <- sd(dat$Cytotoxic_cohort_SD)
  data.table(model = if (adjusted) "Regimen-adjusted" else "Unadjusted",
             adjustment = adjustment, n_complete = nrow(dat),
             n_R = sum(dat$response_binary == 1), n_NR = sum(dat$response_binary == 0),
             score_scale = "one SD in all Gide baseline patients",
             OR_per_SD = est["OR"], CI_low = est["CI_low"], CI_high = est["CI_high"],
             C1_25_patient_SD_in_cohort_SD_units = c1_sd,
             OR_per_C1_25_patient_SD = exp(log(est["OR"]) * c1_sd),
             CI_low_per_C1_25_patient_SD = est["CI_low"] ^ c1_sd,
             CI_high_per_C1_25_patient_SD = est["CI_high"] ^ c1_sd,
             interval = "95% profile penalized-likelihood CI")
}
firth_rows <- rbindlist(list(fit_one(cc, FALSE), fit_one(cc, TRUE)))
fwrite(firth_rows, file.path(res, "Gide_C1_cytotoxic_Firth_regimen_sensitivity.tsv"),
       sep = "\t", na = "NA")

# 散点中可见高 Cytotoxic 响应者；逐一剔除患者检查 regimen 调整点估计是否由单人支配。
loo_rows <- rbindlist(lapply(seq_len(nrow(cc)), function(i) {
  keep <- cc[-i]
  design_data <- data.frame(Cytotoxic_cohort_SD = keep$Cytotoxic_cohort_SD,
                            regimen_group = droplevels(keep$regimen_group))
  X <- model.matrix(~ Cytotoxic_cohort_SD + regimen_group, design_data)
  term <- which(colnames(X) == "Cytotoxic_cohort_SD")
  est <- firth_profile_ci(X, keep$response_binary, term, profile = FALSE)
  data.table(omitted_patient_id = cc$patient_id[i],
             omitted_response = cc$strict_response[i],
             omitted_regimen = as.character(cc$regimen_group[i]),
             omitted_Cytotoxic = cc$Cytotoxic[i],
             n_remaining = nrow(keep),
             OR_per_cohort_SD = est["OR"],
             OR_per_full_C1_SD = exp(est["logOR"] * sd(cc$Cytotoxic)))
}))
fwrite(loo_rows, file.path(res, "Gide_C1_cytotoxic_regimen_adjusted_leave_one_out.tsv"),
       sep = "\t", na = "NA")
lo_idx <- which.min(loo_rows$OR_per_cohort_SD)
hi_idx <- which.max(loo_rows$OR_per_cohort_SD)
loo_summary <- data.table(
  n_leave_one_out_fits = nrow(loo_rows),
  full_data_OR_per_cohort_SD = firth_rows[model == "Regimen-adjusted", OR_per_SD],
  leave_one_out_OR_min = loo_rows$OR_per_cohort_SD[lo_idx],
  omitted_patient_at_min = loo_rows$omitted_patient_id[lo_idx],
  leave_one_out_OR_max = loo_rows$OR_per_cohort_SD[hi_idx],
  omitted_patient_at_max = loo_rows$omitted_patient_id[hi_idx],
  full_C1_SD_in_cohort_SD_units = sd(cc$Cytotoxic),
  interpretation = "point-estimate influence check only; no p values or intervals"
)
fwrite(loo_summary, file.path(res, "Gide_C1_cytotoxic_leave_one_out_summary.tsv"),
       sep = "\t", na = "NA")

# 以患者为单位 bootstrap Spearman 区间；该轴是既有 11 个非 AMP LM12 ssGSEA 分数均值。
set.seed(20260923)
n_boot <- 5000L
rho <- suppressWarnings(cor(x$Cytotoxic, x$LM12_inflammatory_axis, method = "spearman"))
boot_rho <- replicate(n_boot, {
  idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
  suppressWarnings(cor(x$Cytotoxic[idx], x$LM12_inflammatory_axis[idx], method = "spearman"))
})
rho_ci <- quantile(boot_rho[is.finite(boot_rho)], c(.025, .975), names = FALSE)
corr_result <- data.table(
  cohort = "Gide", C_state = "C1", n = nrow(x),
  x_score = "Cytotoxic effector programme score",
  y_score = "LM12 inflammatory axis (mean of 11 non-AMP ssGSEA scores)",
  spearman_rho = rho, CI_low = rho_ci[1], CI_high = rho_ci[2],
  interval = paste0("95% percentile patient bootstrap; B=", n_boot),
  seed = 20260923L, finite_bootstrap = sum(is.finite(boot_rho)),
  p_value = NA_real_, multiplicity_adjustment = "not performed; descriptive association only"
)
fwrite(corr_result, file.path(res, "Gide_C1_cytotoxic_LM12_axis_correlation.tsv"),
       sep = "\t", na = "NA")

# 检查 Cytotoxic 六基因是否直接参与 LM12 轴成员定义，以限定相关性的解释。
program_members <- fread(unname(cfg_paths["frozen_program_members"]))
if (!all(c("program", "gene") %in% names(program_members))) {
  stop("Unexpected frozen programme member table schema")
}
cytogenes <- unique(toupper(program_members[program == "Cytotoxic", gene]))
lm_long <- fread(unname(cfg_paths["lm12_members"]))
axis_genes <- unique(toupper(lm_long[feature_id != "LM12_Antimicrobial_Peptides", gene]))
overlap <- data.table(gene = sort(cytogenes),
                      member_of_LM12_inflammatory_axis = sort(cytogenes) %in% axis_genes,
                      axis_definition = "11 non-AMP LM12 gene sets")
fwrite(overlap, file.path(res, "Gide_C1_cytotoxic_LM12_member_overlap.tsv"),
       sep = "\t", na = "NA")

# 保留本次实际纳入的 25 名患者及用于本探索的冻结分数/临床字段。
patient_subset <- x[, .(patient_id, sample_id, strict_response, regimen_source = regimen,
                         regimen_group, Cytotoxic, LM12_inflammatory_axis, Acid_Score)]
setorder(patient_subset, patient_id)
fwrite(patient_subset, file.path(res, "Gide_C1_patients_used.tsv"), sep = "\t", na = "NA")

# 单张诊断图：散点同时编码疗效与方案，旁图呈现各方案下的 Cytotoxic 分布。
response_cols <- c("SD/PD" = "#5279a6", "CR/PR" = "#b85a50")
regimen_shapes <- c("Monotherapy" = 16, "Combination" = 17)
theme_explore <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9),
        legend.position = "bottom", legend.title = element_text(face = "bold"))
p_axis <- ggplot(x, aes(LM12_inflammatory_axis, Cytotoxic,
                        colour = response_label, shape = regimen_group)) +
  geom_point(size = 2.8, alpha = .9) +
  scale_colour_manual(values = response_cols, drop = FALSE) +
  scale_shape_manual(values = regimen_shapes, drop = FALSE) +
  labs(title = "Gide C1: cytotoxicity and inflammation scores",
       subtitle = sprintf("n=25; Spearman rho=%.2f (95%% bootstrap CI %.2f to %.2f)",
                          rho, rho_ci[1], rho_ci[2]),
       x = "LM12 inflammatory axis (mean of 11 ssGSEA scores)",
       y = "Cytotoxic programme score (cohort SD)",
       colour = "Best response", shape = "Pretreatment regimen") +
  theme_explore

p_regimen <- ggplot(x, aes(regimen_group, Cytotoxic,
                           colour = response_label)) +
  geom_boxplot(aes(group = regimen_group), width = .45,
               outlier.shape = NA, alpha = .15, colour = "#555555") +
  geom_point(position = position_jitter(width = .09, height = 0, seed = 20260923),
             size = 2.5, alpha = .9) +
  scale_colour_manual(values = response_cols, drop = FALSE) +
  scale_x_discrete(labels = c("Monotherapy" = paste0("Anti-PD-1\nmonotherapy\nn=", 
                                                       regimen_summary[regimen_group == "Monotherapy", n]),
                              "Combination" = paste0("Anti-PD-1 +\nanti-CTLA-4\nn=",
                                                      regimen_summary[regimen_group == "Combination", n])),
                   drop = FALSE) +
  labs(title = "Cytotoxic scores by regimen",
       subtitle = "Each point is one pretreatment C1 patient",
       x = NULL, y = "Cytotoxic programme score (cohort SD)",
       colour = "Best response") + theme_explore

draw_diagnostic <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.15, 1), "null"))))
  print(p_axis, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_regimen, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  popViewport()
}
png(file.path(fig, "Gide_C1_regimen_axis_diagnostic.png"),
    width = 2600, height = 1350, res = 180)
draw_diagnostic(); dev.off()
pdf(file.path(fig, "Gide_C1_regimen_axis_diagnostic.pdf"), width = 13, height = 6.75)
draw_diagnostic(); dev.off()

expected <- c(
  "Gide_C1_regimen_response_distribution.tsv",
  "Gide_C1_cytotoxic_by_regimen_response.tsv",
  "Gide_C1_cytotoxic_Firth_regimen_sensitivity.tsv",
  "Gide_C1_cytotoxic_regimen_adjusted_leave_one_out.tsv",
  "Gide_C1_cytotoxic_leave_one_out_summary.tsv",
  "Gide_C1_cytotoxic_LM12_axis_correlation.tsv",
  "Gide_C1_cytotoxic_LM12_member_overlap.tsv",
  "Gide_C1_patients_used.tsv")
if (any(!file.exists(file.path(res, expected))) ||
    any(file.info(file.path(res, expected))$size <= 0) ||
    any(file.info(file.path(fig, c("Gide_C1_regimen_axis_diagnostic.png",
                                  "Gide_C1_regimen_axis_diagnostic.pdf")))$size <= 0)) {
  stop("One or more exploratory outputs are missing or empty")
}
log_event("CHECKPOINT", sprintf(
  "Gide C1 n=%d; R=%d NR=%d; regimen levels=%s; cytotoxic-LM12 rho=%.3f",
  nrow(x), sum(x$strict_response == "R"), sum(x$strict_response == "NR"),
  paste(levels(droplevels(x$regimen_group)), collapse = "/"), rho))
log_event("DONE", paste0("outputs written under ", out,
                         "; exploratory sensitivity only; fixed C1 three-program results unchanged"))
