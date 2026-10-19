#!/usr/bin/env Rscript
# Purpose: Compare C1-C4 clinical features and overall-survival models.
# Task: c1c4_clinical; step: s01
# Inputs (DCF fields): patient_manifest,canonical_labels,clinical_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_clinical --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
usage <- function() {
  cat("Usage: Rscript 01_run_stage1_clinical_models.R --config PATH [--check-config]\n",
      "       Rscript 01_run_stage1_clinical_models.R --help\n\n",
      "DCF fields: patient_manifest, canonical_labels, clinical_data, output_dir.\n",
      "Relative paths resolve from the configuration file's directory.\n",
      "--check-config checks file metadata only: no input loading, dependencies or output creation.\n",
      "An ordinary configured invocation runs clinical tests, Cox models, meta-analysis and plots.\n",
      sep = "")
}
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  usage()
  quit(status = 0L)
}
if (!(length(args) %in% c(2L, 3L)) || args[1L] != "--config" || !nzchar(args[2L]) ||
    (length(args) == 3L && args[3L] != "--check-config")) {
  usage()
  stop("Expected --config PATH, optionally followed by --check-config.", call. = FALSE)
}
flags <- if (length(args) == 3L) args[3L] else character()

# Resolve DCF-relative paths through existing ancestors and symlinks.
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
config_path <- normalizePath(args[2L], winslash = "/", mustWork = TRUE)
if (!file_test("-f", config_path) || file.access(config_path, mode = 4L) != 0L ||
    is.na(file.info(config_path)$size) || file.info(config_path)$size <= 0) {
  stop("Configuration must be a readable nonempty regular file.", call. = FALSE)
}
cfg <- read.dcf(config_path)
required_fields <- c("patient_manifest", "canonical_labels", "clinical_data", "output_dir")
if (nrow(cfg) != 1L || !setequal(colnames(cfg), required_fields) ||
    ncol(cfg) != length(required_fields)) {
  stop("Configuration must contain one DCF record with exactly: ",
       paste(required_fields, collapse = ", "), call. = FALSE)
}
cfg <- setNames(as.character(cfg[1L, required_fields]), required_fields)
if (anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME", cfg, fixed = TRUE))) {
  stop("Replace all example placeholders with explicit paths.", call. = FALSE)
}
paths <- vapply(cfg, resolve_path, character(1L), base_dir = dirname(config_path))
input_file <- unname(paths["patient_manifest"])
canonical_file <- unname(paths["canonical_labels"])
clinical_file <- unname(paths["clinical_data"])
out_dir <- unname(paths["output_dir"])
required_files <- c(input_file, canonical_file, clinical_file)
# Resolve each input's final target before checking output independence.
required_files <- vapply(required_files, normalizePath, character(1L),
                        winslash = "/", mustWork = TRUE)

# Protect inputs and existing results; checks create no directories.
protected_paths <- c(dirname(required_files), config_path)
if (out_dir == "/" || out_dir == normalizePath(path.expand("~")) ||
    any(vapply(protected_paths, function(p) {
      is_within(p, out_dir) || is_within(out_dir, p)
    }, logical(1L)))) {
  stop("output_dir must be independent of inputs and the configuration file.", call. = FALSE)
}
if (file.exists(out_dir) && !dir.exists(out_dir)) {
  stop("output_dir exists but is not a directory.", call. = FALSE)
}
if (dir.exists(out_dir)) {
  if (file.access(out_dir, mode = 7L) != 0L ||
      length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    stop("output_dir must be accessible, writable and empty.", call. = FALSE)
  }
}
if (any(!file_test("-f", required_files)) ||
    any(file.access(required_files, mode = 4L) != 0L) ||
    any(is.na(file.info(required_files)$size) | file.info(required_files)$size <= 0)) {
  stop("All inputs must be readable nonempty regular files.", call. = FALSE)
}
existing_parent <- out_dir
while (!dir.exists(existing_parent)) {
  if (file.exists(existing_parent)) {
    stop("output_dir has a non-directory ancestor.", call. = FALSE)
  }
  existing_parent <- dirname(existing_parent)
}
if (file.access(existing_parent, mode = 3L) != 0L) {
  stop("output_dir has no writable searchable existing parent.", call. = FALSE)
}
if ("--check-config" %in% flags) {
  cat("Configuration paths passed metadata checks. No inputs loaded; no analysis or outputs created.\n")
  quit(status = 0L)
}

set.seed(20260722)
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})
script_name <- "01_run_stage1_clinical_models.R"
result_dir <- file.path(out_dir, "results")
figure_dir <- file.path(out_dir, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

version_tag <- "v0.1_20260722"
start_time <- Sys.time()

log_event <- function(level, msg, ...) {
  extra <- list(...)
  extra_text <- if (length(extra) == 0L) "" else paste(
    vapply(names(extra), function(nm) sprintf("%s=\"%s\"", nm, as.character(extra[[nm]])), character(1)),
    collapse = " "
  )
  cat(sprintf("[%s] %s script=%s msg=\"%s\"%s%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name, msg,
              ifelse(nzchar(extra_text), " ", ""), extra_text))
  flush.console()
}

fail <- function(msg) {
  log_event("ERROR", msg)
  quit(status = 1)
}

is_missing_text <- function(x) {
  z <- toupper(trimws(as.character(x)))
  is.na(z) | z %in% c("", "NA", "N/A", "NULL", "[NOT AVAILABLE]", "NOT REPORTED", "[UNKNOWN]")
}

clean_event <- function(x) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  z[!(z %in% c(0, 1))] <- NA_real_
  z
}

clean_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_gender <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[!(z %in% c("FEMALE", "MALE"))] <- NA_character_
  factor(z, levels = c("FEMALE", "MALE"))
}

clean_race_bin <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[is_missing_text(z)] <- NA_character_
  out <- ifelse(is.na(z), NA_character_, ifelse(z == "WHITE", "WHITE", "NON_WHITE"))
  factor(out, levels = c("WHITE", "NON_WHITE"))
}

normalize_stage4_raw <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[is_missing_text(z)] <- NA_character_
  z <- gsub("[[:space:]]+", "", z)
  z <- gsub("^STAGE", "", z)
  z <- gsub("[\\[\\]()]", "", z)
  z <- gsub("NOS$", "", z)
  out <- rep(NA_character_, length(z))
  out[!is.na(z) & grepl("^I", z)] <- "I"
  out[!is.na(z) & grepl("^II", z)] <- "II"
  out[!is.na(z) & grepl("^III", z)] <- "III"
  out[!is.na(z) & grepl("^IV", z)] <- "IV"
  out
}

standardize_safe <- function(x) {
  x <- clean_numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

stage_fixed_merge <- function(x) {
  z <- as.character(x)
  z[!(z %in% c("I", "II", "III", "IV"))] <- NA_character_
  tab <- table(z, useNA = "no")
  present <- names(tab)[tab > 0]
  if (length(present) < 2L) return(list(value = factor(rep(NA_character_, length(z))), df = 0L, rule = "NOT_USABLE"))
  if (any(tab[present] < 10L)) {
    out <- ifelse(z %in% c("I", "II"), "EARLY", ifelse(z %in% c("III", "IV"), "ADVANCED", NA_character_))
    out <- factor(out, levels = c("EARLY", "ADVANCED"))
    if (length(unique(na.omit(as.character(out)))) < 2L) {
      return(list(value = factor(rep(NA_character_, length(z))), df = 0L, rule = "NOT_USABLE_AFTER_FIXED_MERGE"))
    }
    return(list(value = out, df = 1L, rule = "FIXED_I_II_VS_III_IV"))
  }
  out <- factor(z, levels = c("I", "II", "III", "IV"))
  list(value = droplevels(out), df = length(unique(na.omit(z))) - 1L, rule = "FOUR_LEVEL_I_IV")
}

safe_coxph <- function(formula_obj, data_df) {
  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(
      coxph(formula_obj, data = data_df, ties = "efron", x = TRUE, model = TRUE,
            control = coxph.control(iter.max = 300)),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(fit, "error")) {
    return(list(fit = NULL, error = conditionMessage(fit), warnings = unique(warnings), valid = FALSE))
  }
  cf <- coef(fit)
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  bad_warning <- any(grepl("did not converge|infinite|singular", warnings, ignore.case = TRUE))
  valid <- length(cf) > 0L && all(is.finite(cf)) && !is.null(vc) && all(is.finite(diag(vc))) && !bad_warning
  list(fit = fit, error = ifelse(valid, NA_character_, "NONFINITE_OR_NONCONVERGENT_MODEL"),
       warnings = unique(warnings), valid = valid)
}

safe_zph_global <- function(fit) {
  z <- tryCatch(cox.zph(fit), error = function(e) NULL)
  if (is.null(z)) return(NA_real_)
  tab <- z$table
  if ("GLOBAL" %in% rownames(tab)) as.numeric(tab["GLOBAL", "p"]) else NA_real_
}

model_df_count <- function(design_df, design_type, extra_covariates) {
  design_df <- as.data.table(design_df)
  design_df_count <- if (design_type == "omnibus") 3L else 1L
  extra_df <- 0L
  if ("age10" %in% extra_covariates) extra_df <- extra_df + 1L
  if ("sex" %in% extra_covariates) extra_df <- extra_df + 1L
  if ("purity_z" %in% extra_covariates) extra_df <- extra_df + 1L
  if ("race_bin" %in% extra_covariates) extra_df <- extra_df + 1L
  if ("stage_model" %in% extra_covariates) extra_df <- extra_df + max(0L, length(unique(na.omit(as.character(design_df$stage_model)))) - 1L)
  design_df_count + extra_df
}

design_specs <- list(
  C1C4_OMNIBUS = list(type = "omnibus", variable = "C_subtype", levels = c("C4", "C1", "C2", "C3"), target = NA_character_, reference = "C4"),
  AG1_IG1_vs_AG2_IG1 = list(type = "pair", variable = "AGIG", levels = c("AG2_IG1", "AG1_IG1"), target = "AG1_IG1", reference = "AG2_IG1"),
  AG2_IG1_vs_AG2_IG2 = list(type = "pair", variable = "AGIG", levels = c("AG2_IG2", "AG2_IG1"), target = "AG2_IG1", reference = "AG2_IG2"),
  AG2_IG2_vs_AG3_IG2 = list(type = "pair", variable = "AGIG", levels = c("AG3_IG2", "AG2_IG2"), target = "AG2_IG2", reference = "AG3_IG2"),
  C1_vs_C4 = list(type = "pair", variable = "C_subtype", levels = c("C4", "C1"), target = "C1", reference = "C4"),
  CONTINUOUS_ACID = list(type = "continuous", variable = "Acid_z", levels = NULL, target = "per_1SD", reference = NA_character_),
  CONTINUOUS_IG = list(type = "continuous", variable = "IG_z", levels = NULL, target = "per_1SD", reference = NA_character_)
)

endpoint_specs <- list(
  OS = c("OS", "OS_time"),
  PFI = c("PFI", "PFI_time"),
  DSS = c("DSS", "DSS_time"),
  DFI = c("DFI", "DFI_time")
)

prepare_design_data <- function(dt, spec, endpoint) {
  event_col <- endpoint_specs[[endpoint]][1]
  time_col <- endpoint_specs[[endpoint]][2]
  x <- copy(dt[!is.na(get(event_col)) & !is.na(get(time_col)) & get(time_col) >= 0])
  setnames(x, c(event_col, time_col), c("event", "time"), skip_absent = TRUE)
  if (spec$type %in% c("omnibus", "pair")) {
    x <- x[get(spec$variable) %in% spec$levels]
    x[, exposure := factor(get(spec$variable), levels = spec$levels)]
  } else {
    x <- x[!is.na(get(spec$variable))]
    x[, exposure := get(spec$variable)]
  }
  x
}

check_unadjusted_threshold <- function(x, spec) {
  if (spec$type == "omnibus") {
    n_tab <- table(x$exposure)
    e_tab <- tapply(x$event, x$exposure, sum, na.rm = TRUE)
    ok <- length(n_tab) == 4L && all(n_tab >= 10L) && all(e_tab >= 3L) && sum(x$event) >= 30L
    reason <- sprintf("state_n=%s;state_events=%s;total_events=%d",
                      paste(names(n_tab), n_tab, sep = ":", collapse = "|"),
                      paste(names(e_tab), e_tab, sep = ":", collapse = "|"), sum(x$event))
  } else if (spec$type == "pair") {
    n_tab <- table(x$exposure)
    e_tab <- tapply(x$event, x$exposure, sum, na.rm = TRUE)
    ok <- length(n_tab) == 2L && all(n_tab >= 20L) && all(e_tab >= 5L) && sum(x$event) >= 20L
    reason <- sprintf("group_n=%s;group_events=%s;total_events=%d",
                      paste(names(n_tab), n_tab, sep = ":", collapse = "|"),
                      paste(names(e_tab), e_tab, sep = ":", collapse = "|"), sum(x$event))
  } else {
    ok <- nrow(x) >= 80L && sum(x$event) >= 30L && sd(x$exposure, na.rm = TRUE) > 0
    reason <- sprintf("n=%d;total_events=%d;exposure_sd=%.4f", nrow(x), sum(x$event), sd(x$exposure, na.rm = TRUE))
  }
  list(ok = isTRUE(ok), detail = reason)
}

extract_fit_rows <- function(fit_result, cancer, endpoint, design_id, spec, model_type,
                             analysis_set, n_input, n_used, events, covariates, threshold_detail,
                             stage_rule = NA_character_) {
  if (!fit_result$valid) return(data.table())
  fit <- fit_result$fit
  cf <- coef(fit)
  vc <- vcov(fit)
  if (spec$type == "omnibus") {
    effect_names <- grep("^exposure", names(cf), value = TRUE)
    global_p <- unname(summary(fit)$logtest["pvalue"])
  } else {
    effect_names <- if (spec$type == "pair") grep("^exposure", names(cf), value = TRUE) else "exposure"
    global_p <- NA_real_
  }
  if (length(effect_names) == 0L) return(data.table())
  rows <- rbindlist(lapply(effect_names, function(nm) {
    beta <- unname(cf[nm])
    se <- sqrt(vc[nm, nm])
    effect_level <- if (spec$type == "omnibus") sub("^exposure", "", nm) else spec$target
    data.table(
      Cancer = cancer, Endpoint = endpoint, Design_ID = design_id, Design_Type = spec$type,
      Model_Type = model_type, Analysis_Set = analysis_set, Effect_Level = effect_level,
      Reference_Level = spec$reference, Beta = beta, SE = se, HR = exp(beta),
      CI_Lower = exp(beta - 1.96 * se), CI_Upper = exp(beta + 1.96 * se),
      P_Value = 2 * pnorm(abs(beta / se), lower.tail = FALSE), Global_P = global_p,
      N_Input = n_input, N_Used = n_used, Events = events, Model_DF = length(cf),
      Covariates = covariates, Threshold_Detail = threshold_detail, Stage_Rule = stage_rule,
      PH_Global_P = safe_zph_global(fit), Converged = TRUE,
      Warnings = ifelse(length(fit_result$warnings) == 0L, NA_character_, paste(fit_result$warnings, collapse = " | ")),
      Status = "ESTIMATED", Failure_Reason = NA_character_
    )
  }))
  rows
}

empty_check_row <- function(cancer, endpoint, design_id, spec, model_type, analysis_set,
                            x, status, reason, threshold_detail, completeness = NA_real_, events_complete = NA_integer_,
                            stage_rule = NA_character_) {
  data.table(
    Cancer = cancer, Endpoint = endpoint, Design_ID = design_id, Design_Type = spec$type,
    Model_Type = model_type, Analysis_Set = analysis_set, N_Input = nrow(x),
    Events_Input = sum(x$event, na.rm = TRUE), Complete_Rate = completeness,
    Complete_Events = events_complete, Status = status, Reason = reason,
    Threshold_Detail = threshold_detail, Stage_Rule = stage_rule
  )
}

run_one_model <- function(x, cancer, endpoint, design_id, spec, model_type, analysis_set, threshold_detail) {
  extra <- character()
  min_complete_rate <- NA_real_
  stage_rule <- NA_character_
  if (model_type == "CORE_AGE_SEX") {
    extra <- c("age10", "sex")
    min_complete_rate <- 0.70
  } else if (model_type == "CORE_STAGE") {
    stage_obj <- stage_fixed_merge(x$stage4)
    x[, stage_model := stage_obj$value]
    stage_rule <- stage_obj$rule
    extra <- c("age10", "sex", "stage_model")
    min_complete_rate <- 0.50
    if (stage_obj$df < 1L) {
      return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
        analysis_set, x, "NOT_ESTIMABLE_BY_CONTRACT", "STAGE_NOT_USABLE", threshold_detail,
        stage_rule = stage_rule)))
    }
  } else if (model_type == "CORE_PURITY") {
    extra <- c("age10", "sex", "purity_z")
    min_complete_rate <- 0.70
  } else if (model_type == "CORE_RACE") {
    extra <- c("age10", "sex", "race_bin")
    min_complete_rate <- 0.70
  }

  needed <- c("time", "event", "exposure", extra)
  complete_idx <- complete.cases(x[, ..needed])
  comp <- x[complete_idx]
  complete_rate <- if (nrow(x) == 0L) 0 else nrow(comp) / nrow(x)
  complete_events <- sum(comp$event, na.rm = TRUE)

  if (length(extra) > 0L && complete_rate < min_complete_rate) {
    return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
      analysis_set, x, "NOT_ESTIMABLE_BY_CONTRACT", "COMPLETE_RATE_BELOW_THRESHOLD", threshold_detail,
      complete_rate, complete_events, stage_rule)))
  }
  if (model_type == "CORE_AGE_SEX" && complete_events < 40L) {
    return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
      analysis_set, x, "NOT_ESTIMABLE_BY_CONTRACT", "CORE_COMPLETE_EVENTS_LT40", threshold_detail,
      complete_rate, complete_events, stage_rule)))
  }
  if ("sex" %in% extra && length(unique(na.omit(as.character(comp$sex)))) < 2L) {
    return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
      analysis_set, x, "NOT_ESTIMABLE_BY_CONTRACT", "SEX_SINGLE_LEVEL", threshold_detail,
      complete_rate, complete_events, stage_rule)))
  }
  if ("race_bin" %in% extra) {
    race_tab <- table(comp$race_bin)
    if (length(race_tab) < 2L || any(race_tab < 10L)) {
      return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
        analysis_set, x, "NOT_ESTIMABLE_BY_CONTRACT", "RACE_RETAINED_LEVEL_LT10", threshold_detail,
        complete_rate, complete_events, stage_rule)))
    }
  }
  df_count <- model_df_count(comp, spec$type, extra)
  if (length(extra) > 0L && complete_events < 10L * df_count) {
    return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
      analysis_set, x, "NOT_ESTIMABLE_BY_CONTRACT", "EVENTS_LT_10_PER_DF", threshold_detail,
      complete_rate, complete_events, stage_rule)))
  }

  rhs <- paste(c("exposure", extra), collapse = " + ")
  formula_obj <- as.formula(paste0("Surv(time, event) ~ ", rhs))
  fit_result <- safe_coxph(formula_obj, comp)
  if (!fit_result$valid) {
    reason <- paste(na.omit(c(fit_result$error, fit_result$warnings)), collapse = " | ")
    return(list(result = data.table(), check = empty_check_row(cancer, endpoint, design_id, spec, model_type,
      analysis_set, x, "MODEL_FAILED", reason, threshold_detail, complete_rate, complete_events, stage_rule)))
  }
  result <- extract_fit_rows(fit_result, cancer, endpoint, design_id, spec, model_type,
                             analysis_set, nrow(x), nrow(comp), complete_events,
                             ifelse(length(extra) == 0L, "exposure", paste(c("exposure", extra), collapse = "+")),
                             threshold_detail, stage_rule)
  check <- empty_check_row(cancer, endpoint, design_id, spec, model_type, analysis_set, x,
                           "ESTIMATED", NA_character_, threshold_detail, complete_rate, complete_events, stage_rule)
  list(result = result, check = check)
}

run_design <- function(dt, cancer, endpoint, design_id, spec, analysis_set) {
  x <- prepare_design_data(dt, spec, endpoint)
  threshold <- check_unadjusted_threshold(x, spec)
  model_types <- c("UNADJUSTED", "CORE_AGE_SEX", "CORE_STAGE", "CORE_PURITY", "CORE_RACE")
  if (analysis_set == "EXCLUDE_DISCORDANT") model_types <- c("UNADJUSTED", "CORE_AGE_SEX")
  if (!threshold$ok) {
    checks <- rbindlist(lapply(model_types, function(mt) {
      empty_check_row(cancer, endpoint, design_id, spec, mt, analysis_set, x,
                      "NOT_ESTIMABLE_BY_CONTRACT", "UNADJUSTED_THRESHOLD_NOT_MET", threshold$detail)
    }))
    return(list(result = data.table(), check = checks))
  }
  out <- lapply(model_types, function(mt) run_one_model(x, cancer, endpoint, design_id, spec, mt, analysis_set, threshold$detail))
  list(result = rbindlist(lapply(out, `[[`, "result"), fill = TRUE),
       check = rbindlist(lapply(out, `[[`, "check"), fill = TRUE))
}

reml_meta <- function(yi, sei) {
  ok <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[ok]
  vi <- sei[ok]^2
  k <- length(yi)
  if (k < 2L) return(NULL)
  objective <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  upper <- max(10, var(yi) * 20, max(vi) * 20)
  opt <- optimize(objective, interval = c(0, upper))
  tau2 <- if (objective(0) <= opt$objective + 1e-10) 0 else opt$minimum
  w <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  se_mu <- sqrt(1 / sum(w))
  z <- mu / se_mu
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  w_fixed <- 1 / vi
  mu_fixed <- sum(w_fixed * yi) / sum(w_fixed)
  q <- sum(w_fixed * (yi - mu_fixed)^2)
  q_df <- k - 1L
  q_p <- pchisq(q, df = q_df, lower.tail = FALSE)
  i2 <- ifelse(q <= 0, 0, max(0, (q - q_df) / q) * 100)
  list(k = k, beta = mu, se = se_mu, p = p, tau2 = tau2, Q = q, Q_df = q_df,
       Q_p = q_p, I2 = i2, pred_low = mu - 1.96 * sqrt(tau2 + se_mu^2),
       pred_high = mu + 1.96 * sqrt(tau2 + se_mu^2))
}

run_meta_tables <- function(effect_dt) {
  meta_input <- effect_dt[Design_Type %in% c("pair", "continuous") & Status == "ESTIMATED" & is.finite(Beta) & is.finite(SE)]
  if (nrow(meta_input) == 0L) return(list(meta = data.table(), loo = data.table()))
  groups <- unique(meta_input[, .(Endpoint, Design_ID, Model_Type, Analysis_Set)])
  meta_rows <- list()
  loo_rows <- list()
  idx <- 0L
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- meta_input[g, on = .(Endpoint, Design_ID, Model_Type, Analysis_Set)]
    m <- reml_meta(z$Beta, z$SE)
    if (is.null(m) || m$k < 3L) next
    idx <- idx + 1L
    direction_prop <- mean(sign(z$Beta) == sign(m$beta))
    meta_rows[[idx]] <- data.table(
      Endpoint = g$Endpoint, Design_ID = g$Design_ID, Model_Type = g$Model_Type,
      Analysis_Set = g$Analysis_Set, K_Cancers = m$k, Beta_REML = m$beta, SE_REML = m$se,
      HR_REML = exp(m$beta), CI_Lower = exp(m$beta - 1.96 * m$se),
      CI_Upper = exp(m$beta + 1.96 * m$se), P_Value = m$p, Tau2 = m$tau2,
      Q = m$Q, Q_DF = m$Q_df, Q_P = m$Q_p, I2 = m$I2,
      Prediction_HR_Lower = exp(m$pred_low), Prediction_HR_Upper = exp(m$pred_high),
      Direction_Consistency = direction_prop,
      Cancers = paste(sort(z$Cancer), collapse = ";")
    )
    for (left_out in z$Cancer) {
      zl <- z[Cancer != left_out]
      ml <- reml_meta(zl$Beta, zl$SE)
      if (is.null(ml)) next
      loo_rows[[length(loo_rows) + 1L]] <- data.table(
        Endpoint = g$Endpoint, Design_ID = g$Design_ID, Model_Type = g$Model_Type,
        Analysis_Set = g$Analysis_Set, Left_Out_Cancer = left_out, K_Cancers = ml$k,
        Beta_REML = ml$beta, SE_REML = ml$se, HR_REML = exp(ml$beta),
        CI_Lower = exp(ml$beta - 1.96 * ml$se), CI_Upper = exp(ml$beta + 1.96 * ml$se),
        P_Value = ml$p, I2 = ml$I2
      )
    }
  }
  meta <- rbindlist(meta_rows, fill = TRUE)
  loo <- rbindlist(loo_rows, fill = TRUE)
  if (nrow(meta) > 0L) {
    meta[, Meta_FDR := p.adjust(P_Value, method = "BH"), by = .(Endpoint, Model_Type, Analysis_Set)]
    meta[, LOO_Direction_Stable := {
      ll <- loo[.BY, on = .(Endpoint, Design_ID, Model_Type, Analysis_Set)]
      if (nrow(ll) == 0L) FALSE else all(sign(ll$Beta_REML) == sign(Beta_REML))
    }, by = .(Endpoint, Design_ID, Model_Type, Analysis_Set)]
    meta[, LOO_All_Nominal_Significant := {
      ll <- loo[.BY, on = .(Endpoint, Design_ID, Model_Type, Analysis_Set)]
      if (nrow(ll) == 0L) FALSE else all(ll$P_Value < 0.05)
    }, by = .(Endpoint, Design_ID, Model_Type, Analysis_Set)]
  }
  list(meta = meta, loo = loo)
}

clinical_test <- function(dt, cancer, design_id, spec, feature, analysis_set) {
  x <- copy(dt)
  if (analysis_set == "EXCLUDE_DISCORDANT") x <- x[discordant_multi_sample_state == FALSE]
  if (spec$type %in% c("omnibus", "pair")) {
    x <- x[get(spec$variable) %in% spec$levels]
    x[, exposure := factor(get(spec$variable), levels = spec$levels)]
  } else {
    x <- x[!is.na(get(spec$variable))]
    x[, exposure := get(spec$variable)]
  }
  x <- x[!is.na(get(feature))]
  if (nrow(x) < 20L) return(data.table())
  categorical <- feature %in% c("sex", "stage4", "race_bin")
  if (!categorical) {
    if (spec$type == "omnibus") {
      test <- tryCatch(kruskal.test(x[[feature]] ~ x$exposure), error = function(e) NULL)
      method <- "KRUSKAL_WALLIS"
      estimate <- NA_real_
    } else if (spec$type == "pair") {
      test <- tryCatch(wilcox.test(x[[feature]] ~ x$exposure, exact = FALSE), error = function(e) NULL)
      method <- "WILCOXON_RANK_SUM"
      target_med <- median(x[exposure == spec$target][[feature]], na.rm = TRUE)
      ref_med <- median(x[exposure == spec$reference][[feature]], na.rm = TRUE)
      estimate <- target_med - ref_med
    } else {
      test <- tryCatch(cor.test(x[[feature]], x$exposure, method = "spearman", exact = FALSE), error = function(e) NULL)
      method <- "SPEARMAN"
      estimate <- if (is.null(test)) NA_real_ else unname(test$estimate)
    }
  } else {
    if (spec$type %in% c("omnibus", "pair")) {
      tab <- table(x$exposure, x[[feature]])
      if (nrow(tab) < 2L || ncol(tab) < 2L) return(data.table())
      expected <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
      if (all(expected >= 5)) {
        test <- tryCatch(chisq.test(tab, correct = FALSE), error = function(e) NULL)
        method <- "PEARSON_CHISQ"
      } else {
        test <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 10000), error = function(e) NULL)
        method <- "FISHER_SIMULATED"
      }
      estimate <- NA_real_
    } else {
      y <- x[[spec$variable]]
      grp <- droplevels(factor(x[[feature]]))
      if (nlevels(grp) < 2L) return(data.table())
      test <- tryCatch(kruskal.test(y ~ grp), error = function(e) NULL)
      method <- "KRUSKAL_AXIS_BY_CATEGORY"
      estimate <- NA_real_
    }
  }
  if (is.null(test)) return(data.table())
  data.table(Cancer = cancer, Design_ID = design_id, Design_Type = spec$type,
             Analysis_Set = analysis_set, Feature = feature, Method = method,
             N = nrow(x), Estimate = estimate, P_Value = test$p.value)
}

plot_save <- function(p, stem, width_mm, height_mm) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), p, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), p, width = width_mm, height = height_mm, units = "mm", dpi = 600, compression = "lzw")
}

main <- function() {
  log_event("START", "Stage 1 clinical modeling begin", input = input_file, output = result_dir,
            parameters = "patient_unit;OS_primary;PFI_DSS_DFI_secondary;fixed_designs;REML")
  if (!file.exists(input_file)) fail("patient manifest missing")
  if (!file.exists(canonical_file)) fail("canonical label file missing")
  if (!file.exists(clinical_file)) fail("clinical source object missing")
  dt <- fread(input_file, na.strings = c("NA", "", "NaN"))
  required <- c("Patient_ID", "Sample", "Cancer", "C_subtype", "AGIG", "Acid_Score",
                "LM12_inflammatory_axis", "discordant_multi_sample_state", "age", "gender_clean",
                "race_clean", "stage4", "TumorPurity", "OS", "OS_time", "PFI", "PFI_time",
                "DSS", "DSS_time", "DFI", "DFI_time")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0L) fail(paste("missing columns", paste(missing_cols, collapse = ",")))
  if (anyDuplicated(dt$Patient_ID)) fail("patient manifest is not unique by Patient_ID")
  if (!all(dt$C_subtype %in% c("C1", "C2", "C3", "C4"))) fail("invalid C_subtype found")
  canonical <- fread(canonical_file, select = c("Sample", "Cancer", "C_subtype", "AGIG"))
  if (anyDuplicated(canonical$Sample)) fail("canonical label file is not unique by Sample")
  canonical_check <- canonical[dt, on = "Sample", nomatch = 0L]
  if (nrow(canonical_check) != nrow(dt)) fail("patient representative samples are not fully covered by canonical labels")
  if (any(canonical_check$Cancer != canonical_check$i.Cancer |
          canonical_check$C_subtype != canonical_check$i.C_subtype |
          canonical_check$AGIG != canonical_check$i.AGIG)) {
    fail("patient manifest labels do not match canonical labels")
  }
  clinical_env <- new.env(parent = emptyenv())
  loaded_objects <- load(clinical_file, envir = clinical_env)
  if (!("pancancer_data_list" %in% loaded_objects)) fail("pancancer_data_list not found in clinical source object")
  clinical_list <- clinical_env[["pancancer_data_list"]]
  stage_rows <- rbindlist(lapply(names(clinical_list), function(cancer) {
    cl <- as.data.table(clinical_list[[cancer]][["clinical_data"]])
    needed_stage <- c("sample", "_PATIENT", "ajcc_pathologic_tumor_stage", "clinical_stage")
    if (!all(needed_stage %in% names(cl))) return(data.table())
    ajcc <- normalize_stage4_raw(cl[["ajcc_pathologic_tumor_stage"]])
    clinical <- normalize_stage4_raw(cl[["clinical_stage"]])
    data.table(
      Cancer = cancer,
      Sample = substr(as.character(cl[["sample"]]), 1, 16),
      Patient_ID = substr(as.character(cl[["_PATIENT"]]), 1, 12),
      stage4_rebuilt = ifelse(!is.na(ajcc), ajcc, clinical),
      stage4_source = ifelse(!is.na(ajcc), "ajcc_pathologic_tumor_stage",
                             ifelse(!is.na(clinical), "clinical_stage", NA_character_))
    )
  }), fill = TRUE)
  stage_rows <- unique(stage_rows, by = c("Cancer", "Sample"))
  stage_key <- paste(stage_rows$Cancer, stage_rows$Sample, sep = "|")
  patient_stage_key <- paste(dt$Cancer, dt$Sample, sep = "|")
  stage_idx <- match(patient_stage_key, stage_key)
  dt[, stage4_stage0_interface := stage4]
  dt[, stage4 := stage_rows$stage4_rebuilt[stage_idx]]
  dt[, stage4_source := stage_rows$stage4_source[stage_idx]]
  dt[, age10 := clean_numeric(age) / 10]
  dt[, sex := clean_gender(gender_clean)]
  dt[, race_bin := clean_race_bin(race_clean)]
  dt[, purity_z := standardize_safe(TumorPurity), by = Cancer]
  dt[, Acid_z := standardize_safe(Acid_Score), by = Cancer]
  dt[, IG_z := standardize_safe(LM12_inflammatory_axis), by = Cancer]
  dt[, discordant_multi_sample_state := as.logical(discordant_multi_sample_state)]
  for (ep in names(endpoint_specs)) {
    dt[, (endpoint_specs[[ep]][1]) := clean_event(get(endpoint_specs[[ep]][1]))]
    dt[, (endpoint_specs[[ep]][2]) := clean_numeric(get(endpoint_specs[[ep]][2]))]
  }
  fwrite(dt, file.path(result_dir, paste0("stage1_patient_analysis_manifest_", version_tag, ".tsv")), sep = "\t", na = "NA")
  log_event("CHECKPOINT", "patient interface validated", n_patients = nrow(dt), n_cancers = uniqueN(dt$Cancer),
            discordant_patients = sum(dt$discordant_multi_sample_state), stage_complete = sum(!is.na(dt$stage4)))

  model_results <- list()
  model_checks <- list()
  counter <- 0L
  cancers <- sort(unique(dt$Cancer))
  for (cancer in cancers) {
    cancer_dt <- dt[Cancer == cancer]
    for (endpoint in names(endpoint_specs)) {
      for (design_id in names(design_specs)) {
        spec <- design_specs[[design_id]]
        for (analysis_set in c("REPRESENTATIVE_MAIN", "EXCLUDE_DISCORDANT")) {
          use_dt <- if (analysis_set == "REPRESENTATIVE_MAIN") cancer_dt else cancer_dt[discordant_multi_sample_state == FALSE]
          z <- run_design(use_dt, cancer, endpoint, design_id, spec, analysis_set)
          counter <- counter + 1L
          model_results[[counter]] <- z$result
          model_checks[[counter]] <- z$check
        }
      }
    }
    log_event("CHECKPOINT", "cancer models complete", cancer = cancer)
  }
  effects <- rbindlist(model_results, fill = TRUE)
  checks <- rbindlist(model_checks, fill = TRUE)
  if (nrow(effects) == 0L) fail("no Cox models estimated")
  effects[, FDR := p.adjust(P_Value, method = "BH"), by = .(Endpoint, Design_ID, Model_Type, Analysis_Set, Effect_Level)]
  omnibus_global <- unique(effects[Design_Type == "omnibus",
    .(Cancer, Endpoint, Design_ID, Model_Type, Analysis_Set, Global_P, N_Used, Events, PH_Global_P)])
  omnibus_global[, Global_FDR := p.adjust(Global_P, method = "BH"),
                 by = .(Endpoint, Design_ID, Model_Type, Analysis_Set)]
  fwrite(effects, file.path(result_dir, paste0("cox_effects_by_cancer_", version_tag, ".tsv")), sep = "\t", na = "NA")
  fwrite(omnibus_global, file.path(result_dir, paste0("cox_omnibus_global_by_cancer_", version_tag, ".tsv")), sep = "\t", na = "NA")
  fwrite(checks, file.path(result_dir, paste0("cox_model_evaluability_and_status_", version_tag, ".tsv")), sep = "\t", na = "NA")
  log_event("CHECKPOINT", "cancer-level Cox models written", n_effect_rows = nrow(effects), n_status_rows = nrow(checks))

  meta_out <- run_meta_tables(effects)
  fwrite(meta_out$meta, file.path(result_dir, paste0("random_effects_meta_", version_tag, ".tsv")), sep = "\t", na = "NA")
  fwrite(meta_out$loo, file.path(result_dir, paste0("leave_one_cancer_out_meta_", version_tag, ".tsv")), sep = "\t", na = "NA")
  log_event("CHECKPOINT", "random-effects meta and leave-one-out complete", n_meta = nrow(meta_out$meta), n_loo = nrow(meta_out$loo))

  clinical_rows <- list()
  k <- 0L
  clinical_features <- c("age", "sex", "stage4", "race_bin", "TumorPurity")
  for (cancer in cancers) {
    zdt <- dt[Cancer == cancer]
    for (design_id in names(design_specs)) {
      for (feature in clinical_features) {
        for (analysis_set in c("REPRESENTATIVE_MAIN", "EXCLUDE_DISCORDANT")) {
          k <- k + 1L
          clinical_rows[[k]] <- clinical_test(zdt, cancer, design_id, design_specs[[design_id]], feature, analysis_set)
        }
      }
    }
  }
  clinical <- rbindlist(clinical_rows, fill = TRUE)
  clinical[, FDR := p.adjust(P_Value, method = "BH"), by = .(Feature, Design_ID, Method, Analysis_Set)]
  fwrite(clinical, file.path(result_dir, paste0("clinical_feature_associations_", version_tag, ".tsv")), sep = "\t", na = "NA")

  endpoint_summary <- rbindlist(lapply(names(endpoint_specs), function(ep) {
    e_col <- endpoint_specs[[ep]][1]
    t_col <- endpoint_specs[[ep]][2]
    dt[!is.na(get(e_col)) & !is.na(get(t_col)) & get(t_col) >= 0,
       .(N = .N, Events = sum(get(e_col))), by = .(Cancer, C_subtype)][, Endpoint := ep]
  }), fill = TRUE)
  setcolorder(endpoint_summary, c("Cancer", "Endpoint", "C_subtype", "N", "Events"))
  fwrite(endpoint_summary, file.path(result_dir, paste0("endpoint_state_distribution_", version_tag, ".tsv")), sep = "\t", na = "NA")

  if (nrow(meta_out$meta) > 0L) {
    forest_src <- meta_out$meta[Endpoint == "OS" & Model_Type %in% c("UNADJUSTED", "CORE_AGE_SEX") & Analysis_Set == "REPRESENTATIVE_MAIN"]
    fwrite(forest_src, file.path(result_dir, paste0("figure_source_os_meta_forest_", version_tag, ".tsv")), sep = "\t", na = "NA")
    if (nrow(forest_src) > 0L) {
      design_labels <- c(
        AG1_IG1_vs_AG2_IG1 = "AG1_IG1 vs AG2_IG1",
        AG2_IG1_vs_AG2_IG2 = "AG2_IG1 vs AG2_IG2",
        AG2_IG2_vs_AG3_IG2 = "AG2_IG2 vs AG3_IG2",
        C1_vs_C4 = "C1 vs C4",
        CONTINUOUS_ACID = "Acid axis (per SD)",
        CONTINUOUS_IG = "LM12 inflammatory-remodeling axis (per SD)"
      )
      forest_src[, Design_Label := factor(design_labels[Design_ID], levels = rev(unname(design_labels)))]
      forest_src[, Model_Label := factor(Model_Type, levels = c("UNADJUSTED", "CORE_AGE_SEX"), labels = c("Unadjusted", "Age + sex"))]
      p <- ggplot(forest_src, aes(x = HR_REML, y = Design_Label, colour = Model_Label)) +
        geom_vline(xintercept = 1, linetype = 2, colour = "#777777", linewidth = 0.4) +
        geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), orientation = "y", width = 0,
                      position = position_dodge(width = 0.45), linewidth = 0.55) +
        geom_point(position = position_dodge(width = 0.45), size = 2.1) +
        scale_x_log10() +
        scale_colour_manual(values = c("Unadjusted" = "#2C7FB8", "Age + sex" = "#D95F0E")) +
        labs(x = "Random-effects hazard ratio (95% CI)", y = NULL, colour = NULL) +
        theme_classic(base_size = 8, base_family = "sans") +
        theme(legend.position = "top", axis.text.y = element_text(colour = "black"))
      plot_save(p, paste0("os_random_effects_meta_forest_", version_tag), 125, 78)
    }
  }

  heat_src <- effects[Endpoint == "OS" & Model_Type == "UNADJUSTED" & Analysis_Set == "REPRESENTATIVE_MAIN" &
                        Design_Type %in% c("pair", "continuous"),
                      .(Cancer, Design_ID, Beta, FDR, HR, CI_Lower, CI_Upper)]
  fwrite(heat_src, file.path(result_dir, paste0("figure_source_os_cancer_effect_heatmap_", version_tag, ".tsv")), sep = "\t", na = "NA")
  if (nrow(heat_src) > 0L) {
    heat_src[, Sig := FDR <= 0.05]
    p2 <- ggplot(heat_src, aes(x = Design_ID, y = Cancer, fill = pmax(-2, pmin(2, Beta)))) +
      geom_tile(colour = "white", linewidth = 0.15) +
      geom_point(data = heat_src[Sig == TRUE], shape = 21, size = 1.2, fill = "white", colour = "black", stroke = 0.25) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           limits = c(-2, 2), name = "log(HR)") +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 7, base_family = "sans") +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = "black"), legend.position = "right")
    plot_save(p2, paste0("os_cancer_level_effect_heatmap_", version_tag), 145, 105)
  }

  figure_manifest <- data.table(
    Figure_Stem = c(paste0("os_random_effects_meta_forest_", version_tag), paste0("os_cancer_level_effect_heatmap_", version_tag)),
    PDF = file.path(figure_dir, paste0(c("os_random_effects_meta_forest_", "os_cancer_level_effect_heatmap_"), version_tag, ".pdf")),
    TIFF = file.path(figure_dir, paste0(c("os_random_effects_meta_forest_", "os_cancer_level_effect_heatmap_"), version_tag, ".tiff")),
    Source_Data = file.path(result_dir, paste0(c("figure_source_os_meta_forest_", "figure_source_os_cancer_effect_heatmap_"), version_tag, ".tsv")),
    Role = c("OS cross-cancer random-effects summary", "OS cancer-specific effect overview")
  )
  fwrite(figure_manifest, file.path(result_dir, paste0("figure_manifest_", version_tag, ".tsv")), sep = "\t", na = "NA")

  runtime <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  log_event("DONE", "Stage 1 clinical modeling complete", runtime_sec = runtime,
            n_effect_rows = nrow(effects), n_meta_rows = nrow(meta_out$meta), output = result_dir)
}

tryCatch(main(), error = function(e) fail(conditionMessage(e)))
