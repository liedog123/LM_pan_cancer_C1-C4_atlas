# Purpose: clinical survival analysis stages.
# Configuration: config/clinical.dcf; stage contracts: config/stages.tsv.
# Inputs are read-only; each selected stage runs in a fresh Rscript process.
# 保留原始统计范围、种子、输入输出及阶段边界；共用函数按阶段绑定作用域。

# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/ag_survival/00_prepare_ag_survival_inputs.R
stage_ag_survival_00_prepare_ag_survival_inputs <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Task: ag_survival; step: s00
# Inputs (DCF fields): projects_data,clinical_groups
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task ag_survival --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
projects_file <- unname(paths["projects_data"])
clinical_file <- unname(paths["clinical_groups"])
out_dir <- unname(paths["output_dir"])
required_files <- unname(paths[setdiff(required_fields, "output_dir")])
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。
`load_named` <- bind_stage_helper(shared_load_named_025, environment())
projects <- load_named(projects_file, "projects")
clinical_data_cluster_list <- load_named(clinical_file, "clinical_data_cluster_list")
# 保持原始AG1/AG3分组与行顺序；正式KM及Cox在后续步骤计算。
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
survival_data_list <- list()
for (project in projects) {
  survival_data <- clinical_data_cluster_list[[project]]
  survival_data <- survival_data[survival_data$Group %in% c("High", "Low"), ]
  survival_data_list[[project]] <- survival_data
}
save(survival_data_list, file = file.path(out_dir, "tow_group_clinical_data_list.rdata"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
}

# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/ag_survival/01_run_ag_cox.R
stage_ag_survival_01_run_ag_cox <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Task: ag_survival; step: s01
# Inputs (DCF fields): projects_data,survival_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task ag_survival --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
projects_file <- unname(paths["projects_data"])
survival_file <- unname(paths["survival_data"])
out_dir <- unname(paths["output_dir"])
required_files <- unname(paths[setdiff(required_fields, "output_dir")])
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival) })
`load_named` <- bind_stage_helper(shared_load_named_025, environment())
projects <- load_named(projects_file, "projects")
survival_data_list <- load_named(survival_file, "survival_data_list")
# 原函数保留未使用的gene_list参数；不加载与模型无关的基因对象。
gene_list <- NULL
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 定义用于 Cox 回归分析的函数
run_cox_analysis <- function(project, gene_list) {
  # 准备临床数据
  clinical_data <- survival_data_list[[project]]
  clinical_data$Group <- as.factor(clinical_data$Group)
  clinical_data$Group <- relevel(clinical_data$Group, ref = "Low")

  # 单因素 Cox 回归
  univ_model <- coxph(Surv(OS.time, OS) ~ Group, data = clinical_data, control = coxph.control(iter.max = 300))
  univ_summary <- summary(univ_model)

  # 多因素 Cox 回归数据处理
  # 处理 gender 变量，并转化为因子
  if ("gender" %in% colnames(clinical_data)) {
    unique_genders <- unique(clinical_data$gender)
    if (!("MALE" %in% unique_genders && "FEMALE" %in% unique_genders)) {
      clinical_data <- clinical_data %>% select(-gender)
    } else {
      clinical_data$gender <- as.factor(clinical_data$gender)
      #设置female为参考值
      clinical_data$gender <- relevel(factor(clinical_data$gender), ref = "FEMALE")
    }
  }

  # 处理 race 变量
  if ("race" %in% colnames(clinical_data)) {
    clinical_data <- clinical_data %>%
      mutate(
        race = case_when(
          race == "WHITE" ~ "WHITE",
          TRUE ~ "other"
        )
      )
    # 转换为因子并设置 "white" 为参考值
    clinical_data$race <- factor(clinical_data$race, levels = c("WHITE", "other"))
    clinical_data$race <- relevel(clinical_data$race, ref = "WHITE")
    # 检查 race 是否只包含一个水平
    if (length(unique(clinical_data$race)) == 1) {
      clinical_data <- clinical_data %>% select(-race)
    }
  }


  # 删除 age_at_initial_pathologic_diagnosis 中包含 NA 的行
  if ("age_at_initial_pathologic_diagnosis" %in% colnames(clinical_data)) {
    clinical_data <- clinical_data %>%
      filter(!is.na(age_at_initial_pathologic_diagnosis))
  }

  summary <- summary(clinical_data)

  # 定义协变量
  available_vars <- c("Group", "race", "gender", "age_at_initial_pathologic_diagnosis")
  available_covariates <- intersect(available_vars, colnames(clinical_data))

  # 动态构建 Cox 回归公式
  cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(available_covariates, collapse = " + ")))

  # 拟合模型
  multiv_model <- coxph(cox_formula, data = clinical_data, control = coxph.control(iter.max = 300))
  multiv_summary <- summary(multiv_model)

  # 返回结果
  list(
    univariate = univ_summary,
    multivariate = multiv_summary,
    data_summary = summary
  )
}


# 保存所有结果的列表
results_list <- list()
# 保存每个项目的警告信息
warnings_list <- list()

for (project in projects) {
  warnings <- NULL
  result <- withCallingHandlers(
    tryCatch(
      run_cox_analysis(project, gene_list),
      error = function(e) {
        message(paste("Error in project", project, ":", e$message))
        NULL
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")  # 取消后续显示同一警告
    }
  )

  # 将结果存入结果列表
  results_list[[project]] <- result
  # 将警告信息存入警告列表
  if (!is.null(warnings)) {
    warnings_list[[project]] <- warnings
  }
}


# 保存结果到文件中
save(results_list, file = file.path(out_dir, "Cox_analysis_results.RData"))
save(warnings_list, file = file.path(out_dir, "Cox_warnings.RData"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
}

# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/ag_survival/02_finalize_ag_cox.R
stage_ag_survival_02_finalize_ag_cox <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Task: ag_survival; step: s02
# Inputs (DCF fields): cox_results,cancer_order
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task ag_survival --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
cox_file <- unname(paths["cox_results"])
cancer_order_file <- unname(paths["cancer_order"])
out_dir <- unname(paths["output_dir"])
required_files <- unname(paths[setdiff(required_fields, "output_dir")])
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(readr) })
`load_named` <- bind_stage_helper(shared_load_named_025, environment())
results_list <- load_named(cox_file, "results_list")
panel_dir <- file.path(out_dir, "panel_exports")
source_data_dir <- file.path(out_dir, "source_data")
source_dir <- source_data_dir
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
# 保留完整22癌/两模型家族；图中KM癌种选择不进入Cox/BH。
if (length(results_list) != 22L || any(vapply(results_list, is.null, logical(1L)))) {
  stop("Formal finalization requires all 22 source cancer model results; do not shorten the BH family.")
}

result_cancers <- sub("^TCGA-", "", names(results_list))
shared_order_df <- read_csv(cancer_order_file, show_col_types = FALSE)
cancer_order <- shared_order_df$Cancer[shared_order_df$Cancer %in% result_cancers]
# 仅验证完整输入家族身份/显示顺序；不改变源BH范围。
if (is.null(names(results_list)) || anyNA(names(results_list)) ||
    any(!nzchar(names(results_list))) || anyDuplicated(result_cancers) ||
    anyNA(shared_order_df$Cancer) || anyDuplicated(shared_order_df$Cancer) ||
    length(cancer_order) != 22L || !setequal(cancer_order, result_cancers)) {
  stop("Result cancer names must be unique and fully covered by a unique frozen cancer order.")
}


cox_df <- bind_rows(lapply(names(results_list), function(project) {
  result <- results_list[[project]]
  if (is.null(result)) {
    return(NULL)
  }

  univ_result <- result$univariate
  multiv_result <- result$multivariate

  if ("GroupHigh" %in% rownames(multiv_result$coefficients)) {
    multiv_hr <- exp(multiv_result$coefficients["GroupHigh", "coef"])
    multiv_p <- multiv_result$coefficients["GroupHigh", "Pr(>|z|)"]
  } else {
    multiv_hr <- NA_real_
    multiv_p <- NA_real_
  }

  tibble(
    Cancer_Type = sub("^TCGA-", "", project),
    Model_Type = c("Univariate", "Multivariate"),
    Hazard_Ratio = c(
      univ_result$coefficients[1, "exp(coef)"],
      multiv_hr
    ),
    P_Value = c(
      univ_result$coefficients[1, "Pr(>|z|)"],
      multiv_p
    )
  )
}))

cox_df <- cox_df %>%
  mutate(
    FDR = ave(P_Value, Model_Type, FUN = function(x) p.adjust(x, method = "fdr")),
    Significance = ifelse(FDR <= 0.05, "FDR <= 0.05", "FDR > 0.05"),
    Hazard_Ratio_Discrete = cut(
      Hazard_Ratio,
      breaks = c(-Inf, 0.15, 1.0, 2.0, 3.0, Inf),
      labels = c("<0.15", "0.15-1.0", "1.0-2.0", "2.0-3.0", ">3.0"),
      include.lowest = TRUE
    ),
    Cancer_Type = factor(Cancer_Type, levels = cancer_order),
    Model_Type = factor(Model_Type, levels = c("Multivariate", "Univariate"))
  )

write_csv(
  cox_df,
  file.path(source_data_dir, "Fig1E_cox_summary_values.csv")
)

hr_colors <- c(
  "<0.15" = "#5A88C9",
  "0.15-1.0" = "#B7D8E8",
  "1.0-2.0" = "#F4B6C2",
  "2.0-3.0" = "#F47C7C",
  ">3.0" = "#D7301F"
)

p <- ggplot(
  cox_df,
  aes(x = Cancer_Type, y = Model_Type)
) +
  geom_point(
    aes(fill = Hazard_Ratio_Discrete, size = Significance),
    shape = 22,
    colour = "black",
    stroke = 0.45,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = hr_colors,
    drop = FALSE,
    name = "Hazard ratio"
  ) +
  scale_size_manual(
    values = c("FDR > 0.05" = 4.1, "FDR <= 0.05" = 7.0),
    name = "FDR",
    labels = c("FDR > 0.05" = ">0.05", "FDR <= 0.05" = "<=0.05"),
    guide = guide_legend(
      override.aes = list(shape = 22, fill = "white", colour = "black")
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.6, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(3, 4, 1, 4)
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1, override.aes = list(size = 5.5)),
    size = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(shape = 22, fill = "white", colour = "black"),
      label.hjust = 0.5
    )
  )

ggsave(
  file.path(panel_dir, "Fig1D_cox_summary_v1.pdf"),
  p,
  width = 178 / 25.4,
  height = 55 / 25.4,
  units = "in",
  dpi = 300
)

ggsave(
  file.path(panel_dir, "Fig1D_cox_summary_v1.png"),
  p,
  width = 178 / 25.4,
  height = 55 / 25.4,
  units = "in",
  dpi = 300
)

message("Fig1D Cox summary exported successfully.")

extract_model <- function(project, result, model_type) {
  fit_summary <- if (model_type == "Univariate") result$univariate else result$multivariate
  coef_rows <- rownames(fit_summary$coefficients)
  group_row <- if ("GroupHigh" %in% coef_rows) "GroupHigh" else coef_rows[[1]]
  ci <- fit_summary$conf.int[group_row, , drop = FALSE]
  lower <- unname(ci[1, "lower .95"])
  upper <- unname(ci[1, "upper .95"])
  events <- unname(fit_summary$nevent)
  stability <- case_when(
    !is.finite(lower) | !is.finite(upper) | lower <= 0 ~ "NONFINITE_CI_OR_SEPARATION",
    events < 10 ~ "LOW_EVENT_COUNT_LT10",
    upper / lower > 100 ~ "VERY_WIDE_CI",
    TRUE ~ "OK"
  )
  tibble(
    Cancer_Type = sub("^TCGA-", "", project),
    Model_Type = model_type,
    Extracted_Hazard_Ratio = exp(fit_summary$coefficients[group_row, "coef"]),
    CI95_Lower = lower,
    CI95_Upper = upper,
    Extracted_P_Value = fit_summary$coefficients[group_row, "Pr(>|z|)"],
    Model_N = unname(fit_summary$n),
    Event_N = events,
    Model_Terms = paste(coef_rows, collapse = ";"),
    Stability_Flag = stability
  )
}

d_extract <- bind_rows(lapply(names(results_list), function(project) {
  result <- results_list[[project]]
  if (is.null(result)) return(NULL)
  bind_rows(
    extract_model(project, result, "Univariate"),
    extract_model(project, result, "Multivariate")
  )
}))

d_manifest <- cox_df %>%
  left_join(d_extract, by = c("Cancer_Type", "Model_Type")) %>%
  mutate(
    HR_Match = abs(Hazard_Ratio - Extracted_Hazard_Ratio) < 1e-8 * pmax(1, abs(Hazard_Ratio)),
    P_Match = abs(P_Value - Extracted_P_Value) < 1e-10 * pmax(1, abs(P_Value)),
    Contrast = "High Acid versus Low Acid",
    Endpoint = "overall survival"
  )
stopifnot(nrow(d_manifest) == 44L, all(d_manifest$HR_Match), all(d_manifest$P_Match))
write_tsv(d_manifest, file.path(source_dir, "Fig1E_cox_model_manifest.tsv"))


writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
}

# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/ag_survival/03_export_formal_ag_km.R
stage_ag_survival_03_export_formal_ag_km <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Reproduce AG1-versus-AG3 survival inputs, Cox models and formal KM source data.
# Task: ag_survival; step: s03
# Inputs (DCF fields): survival_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task ag_survival --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
input_rdata <- unname(paths["survival_data"])
out_dir <- unname(paths["output_dir"])
required_files <- unname(paths[setdiff(required_fields, "output_dir")])
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

suppressPackageStartupMessages({ library(survival) })
`load_named` <- bind_stage_helper(shared_load_named_025, environment())
survival_data_list <- load_named(input_rdata, "survival_data_list")
if (!is.list(survival_data_list)) stop("survival_data_list is not a list")
qa_dir <- file.path(out_dir, "qa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

expected_cancers <- c("BLCA", "BRCA", "CESC", "HNSC", "KIRC", "KIRP",
                      "LGG", "LIHC", "LUAD", "SARC", "SKCM")
expected_n <- c(BLCA = 248L, BRCA = 418L, CESC = 173L, HNSC = 353L,
                KIRC = 286L, KIRP = 118L, LGG = 336L, LIHC = 263L,
                LUAD = 287L, SARC = 178L, SKCM = 344L)
expected_low <- c(BLCA = 172L, BRCA = 213L, CESC = 95L, HNSC = 91L,
                  KIRC = 129L, KIRP = 67L, LGG = 181L, LIHC = 214L,
                  LUAD = 246L, SARC = 117L, SKCM = 220L)
expected_high <- c(BLCA = 76L, BRCA = 205L, CESC = 78L, HNSC = 262L,
                   KIRC = 157L, KIRP = 51L, LGG = 155L, LIHC = 49L,
                   LUAD = 41L, SARC = 61L, SKCM = 124L)
expected_events <- c(BLCA = 103L, BRCA = 65L, CESC = 39L, HNSC = 162L,
                     KIRC = 107L, KIRP = 22L, LGG = 100L, LIHC = 85L,
                     LUAD = 92L, SARC = 71L, SKCM = 159L)
expected_p <- c(BLCA = 0.0231412, BRCA = 0.000311019, CESC = 0.000247714,
                HNSC = 0.000851538, KIRC = 0.0196141, KIRP = 0.0175173,
                LGG = 1.28265e-06, LIHC = 8.43172e-05, LUAD = 4.66134e-05,
                SARC = 0.00816214, SKCM = 0.00604262)

resolve_name <- function(cancer) {
  candidates <- c(cancer, paste0("TCGA-", cancer))
  hit <- candidates[candidates %in% names(survival_data_list)]
  if (length(hit) != 1L) {
    stop("Cannot uniquely resolve cancer name: ", cancer,
         "; available matches: ", paste(hit, collapse = ","))
  }
  hit[[1L]]
}

as_numeric_safe <- function(x, field, label) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(out)) && any(!is.na(x))) {
    stop("Cannot coerce ", field, " to numeric for ", label)
  }
  out
}

prepare_cancer <- function(cancer) {
  source_name <- resolve_name(cancer)
  d <- survival_data_list[[source_name]]
  required <- c("Group", "OS", "OS.time")
  missing <- setdiff(required, colnames(d))
  if (length(missing) > 0L) stop("Missing fields for ", source_name, ": ", paste(missing, collapse = ","))
  d <- d[d$Group %in% c("High", "Low"), , drop = FALSE]
  d$Group <- factor(as.character(d$Group), levels = c("Low", "High"))
  d$OS <- as_numeric_safe(d$OS, "OS", source_name)
  d$OS.time <- as_numeric_safe(d$OS.time, "OS.time", source_name)
  d <- d[complete.cases(d[, c("Group", "OS", "OS.time")]), , drop = FALSE]
  d$Cancer <- cancer
  d$panel_id <- sprintf("KM%02d_%s", match(cancer, expected_cancers), cancer)
  d$panel_order <- match(cancer, expected_cancers)
  d$display_label <- cancer
  d$source_scope <- "within_cancer"
  d$sample_id <- rownames(d)
  if (is.null(d$sample_id) || anyNA(d$sample_id) || any(d$sample_id == "")) {
    d$sample_id <- paste0(d$Cancer, "_row_", seq_len(nrow(d)))
  }
  d[, c("panel_order", "panel_id", "display_label", "Cancer", "source_scope",
        "sample_id", "Group", "OS", "OS.time")]
}

get_logrank <- function(d) {
  fit <- survdiff(Surv(OS.time, OS) ~ Group, data = d)
  1 - pchisq(fit$chisq, df = length(fit$n) - 1L)
}

within_data <- do.call(rbind, lapply(expected_cancers, prepare_cancer))
rownames(within_data) <- NULL

# 严格复用 KM_total_samples.R：七个生存字段整体 complete-case 后再拟合 OS。
all_parts <- lapply(names(survival_data_list), function(source_name) {
  d <- survival_data_list[[source_name]]
  required <- c("Group", "OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time")
  missing <- setdiff(required, colnames(d))
  if (length(missing) > 0L) stop("Missing pan-cancer fields for ", source_name, ": ", paste(missing, collapse = ","))
  keep <- d$Group %in% c("High", "Low")
  d <- d[keep, required, drop = FALSE]
  d$Group <- factor(as.character(d$Group), levels = c("Low", "High"))
  for (field in setdiff(required, "Group")) d[[field]] <- as_numeric_safe(d[[field]], field, source_name)
  d$Cancer <- source_name
  d$sample_id <- rownames(d)
  if (is.null(d$sample_id) || anyNA(d$sample_id) || any(d$sample_id == "")) {
    d$sample_id <- paste0(source_name, "_row_", seq_len(nrow(d)))
  }
  d
})
all_survival <- do.call(rbind, all_parts)
all_survival <- all_survival[complete.cases(all_survival[, c("Group", "OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time")]), , drop = FALSE]
pan_data <- all_survival[, c("Group", "OS", "OS.time", "Cancer", "sample_id"), drop = FALSE]
pan_data$panel_order <- 12L
pan_data$panel_id <- "KM12_Pan-cancer"
pan_data$display_label <- "Pan-cancer"
pan_data$source_scope <- "pan_cancer_complete_case_all_endpoints"
pan_data <- pan_data[, c("panel_order", "panel_id", "display_label", "Cancer", "source_scope",
                         "sample_id", "Group", "OS", "OS.time")]

plot_data <- rbind(within_data, pan_data)
plot_data$panel_order <- as.integer(plot_data$panel_order)
plot_data <- plot_data[order(plot_data$panel_order, plot_data$sample_id), , drop = FALSE]
rownames(plot_data) <- NULL

summary_rows <- lapply(seq_len(12L), function(i) {
  panel <- plot_data[plot_data$panel_order == i, , drop = FALSE]
  p <- get_logrank(panel)
  data.frame(
    panel_order = i,
    panel_id = unique(panel$panel_id),
    display_label = unique(panel$display_label),
    source_scope = unique(panel$source_scope),
    n = nrow(panel),
    n_low = sum(panel$Group == "Low"),
    n_high = sum(panel$Group == "High"),
    events = sum(panel$OS == 1),
    logrank_p = p,
    significance_rule = if (i < 12L) "original KM.R log-rank p < 0.05" else "descriptive pan-cancer OS log-rank",
    significant = if (i < 12L) p < 0.05 else NA,
    complete_case_rule = if (i < 12L) "Group + OS + OS.time complete" else "Group + OS + OS.time + PFI + PFI.time + DSS + DSS.time complete",
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)

write.table(plot_data, file.path(path.expand(out_dir), "km_12panel_plot_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_df, file.path(path.expand(out_dir), "km_12panel_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_df[, c("panel_order", "panel_id", "display_label")],
            file.path(path.expand(out_dir), "km_12panel_panel_order.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(path.expand(out_dir), "sessionInfo.txt"))

qa <- data.frame(check_id = character(), value = character(), pass = logical(), stringsAsFactors = FALSE)
add_qa <- function(check_id, value, pass) {
  qa <<- rbind(qa, data.frame(check_id = check_id, value = as.character(value), pass = isTRUE(pass), stringsAsFactors = FALSE))
}

add_qa("input_exists", file.exists(path.expand(input_rdata)), file.exists(path.expand(input_rdata)))
add_qa("panel_count", paste(sort(unique(plot_data$panel_order)), collapse = ","), identical(sort(unique(plot_data$panel_order)), 1:12))
add_qa("within_cancer_order", paste(unique(plot_data$display_label[plot_data$panel_order < 12]), collapse = ","), identical(unique(plot_data$display_label[plot_data$panel_order < 12]), expected_cancers))
add_qa("all_groups_valid", paste(sort(unique(as.character(plot_data$Group))), collapse = ","), all(as.character(plot_data$Group) %in% c("Low", "High")))
add_qa("plot_data_no_missing", sum(!complete.cases(plot_data[, c("panel_order", "panel_id", "Group", "OS", "OS.time")])), all(complete.cases(plot_data[, c("panel_order", "panel_id", "Group", "OS", "OS.time")])))

within_summary <- summary_df[summary_df$panel_order < 12, , drop = FALSE]
within_match <- all(within_summary$n == unname(expected_n[within_summary$display_label])) &&
  all(within_summary$n_low == unname(expected_low[within_summary$display_label])) &&
  all(within_summary$n_high == unname(expected_high[within_summary$display_label])) &&
  all(within_summary$events == unname(expected_events[within_summary$display_label]))
add_qa("within_counts", paste(within_summary$n, collapse = ","), within_match)
add_qa("within_significant_count", sum(within_summary$significant), sum(within_summary$significant) == 11L)
add_qa("within_logrank_p_reproducible", paste(format(within_summary$logrank_p, scientific = TRUE, digits = 6), collapse = ","), all(abs(within_summary$logrank_p - unname(expected_p[within_summary$display_label])) < 1e-07))

pan_summary <- summary_df[summary_df$panel_order == 12, , drop = FALSE]
pan_match <- pan_summary$n == 5188L && pan_summary$n_low == 2639L && pan_summary$n_high == 2549L &&
  pan_summary$events == 1508L && abs(pan_summary$logrank_p - 1.165430382e-04) < 1e-10
add_qa("pan_cancer_counts_and_p", paste(pan_summary$n, pan_summary$n_low, pan_summary$n_high, pan_summary$events, format(pan_summary$logrank_p, scientific = TRUE), sep = ";"), pan_match)

out_files <- c("km_12panel_plot_data.tsv", "km_12panel_summary.tsv", "km_12panel_panel_order.tsv", "sessionInfo.txt")
nonempty <- all(file.info(file.path(path.expand(out_dir), out_files))$size > 0)
add_qa("output_files_nonempty", paste(out_files, collapse = ","), nonempty)
add_qa("all_pass", "see pass column", all(qa$pass))
write.table(qa, file.path(path.expand(qa_dir), "km_12panel_machine_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sprintf("all_qa_pass=%s\noutput_dir=%s\nqa_dir=%s", all(qa$pass), path.expand(out_dir), path.expand(qa_dir)),
           file.path(path.expand(qa_dir), "km_12panel_qa_summary.txt"))

manifest <- data.frame(
  file = out_files,
  bytes = file.info(file.path(path.expand(out_dir), out_files))$size,
  md5 = unname(tools::md5sum(file.path(path.expand(out_dir), out_files))),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(path.expand(qa_dir), "km_12panel_output_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!all(qa$pass)) stop("KM 12-panel machine QA failed; inspect ", path.expand(qa_dir))
message("KM 12-panel figure-ready export and QA completed: ", path.expand(out_dir))
}

# Purpose: Compare C1-C4 clinical features and overall-survival models.
# Source at 41670141faa1902e623bebcb725ec79e33b0b563: internal/clinical_outcomes/01_run_stage1_clinical_models.R
stage_clinical_outcomes_01_run_stage1_clinical_models <- function(cfg, context, check_only = FALSE) {
  prepare_stage(cfg, context, check_only, environment())
#!/usr/bin/env Rscript
# Purpose: Compare C1-C4 clinical features and overall-survival models.
# Task: c1c4_clinical; step: s01
# Inputs (DCF fields): patient_manifest,canonical_labels,clinical_data
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task c1c4_clinical --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

options(stringsAsFactors = FALSE, warn = 1)
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
  if (check_only) return(invisible(TRUE))
  # 以下保留研究计算正文；不在加载模块或检查配置时执行。

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
`clean_numeric` <- bind_stage_helper(shared_clean_numeric_038, environment())

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
}
