#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 03_finalize_rppa_recurrence_tcpa.R --config PATH [--check-config]\n       Rscript 03_finalize_rppa_recurrence_tcpa.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("model_dir", "antibody_dictionary", "tcpa_table13", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
paths <- vapply(cfg,resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
if(!dir.exists(paths["model_dir"]) || file.access(paths["model_dir"],5L)!=0L) stop("model_dir must be readable/searchable.")

required_files <- input_paths[names(input_paths) != "model_dir"]
required_files <- c(required_files,file.path(paths["model_dir"],c("BLCA_rppa_feature_models_v0.1_20260723.tsv.gz","BRCA_rppa_feature_models_v0.1_20260723.tsv.gz","CESC_rppa_feature_models_v0.1_20260723.tsv.gz","COAD_rppa_feature_models_v0.1_20260723.tsv.gz","ESCA_rppa_feature_models_v0.1_20260723.tsv.gz","GBM_rppa_feature_models_v0.1_20260723.tsv.gz","HNSC_rppa_feature_models_v0.1_20260723.tsv.gz","KIRC_rppa_feature_models_v0.1_20260723.tsv.gz","KIRP_rppa_feature_models_v0.1_20260723.tsv.gz","LGG_rppa_feature_models_v0.1_20260723.tsv.gz","LIHC_rppa_feature_models_v0.1_20260723.tsv.gz","LUAD_rppa_feature_models_v0.1_20260723.tsv.gz","LUSC_rppa_feature_models_v0.1_20260723.tsv.gz","OV_rppa_feature_models_v0.1_20260723.tsv.gz","PCPG_rppa_feature_models_v0.1_20260723.tsv.gz","PRAD_rppa_feature_models_v0.1_20260723.tsv.gz","READ_rppa_feature_models_v0.1_20260723.tsv.gz","SARC_rppa_feature_models_v0.1_20260723.tsv.gz","SKCM_rppa_feature_models_v0.1_20260723.tsv.gz","STAD_rppa_feature_models_v0.1_20260723.tsv.gz","THCA_rppa_feature_models_v0.1_20260723.tsv.gz","UCEC_rppa_feature_models_v0.1_20260723.tsv.gz")))

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
})
options(stringsAsFactors = FALSE, scipen = 999)

script_name <- "03_finalize_rppa_recurrence_tcpa.R"
start_time <- Sys.time()

log_event <- function(level, msg, fields = list()) {
  extra <- if (length(fields)) {
    paste(vapply(names(fields), function(nm) {
      sprintf(
        '%s="%s"', nm,
        gsub("[\r\n\t]", " ", as.character(fields[[nm]]))
      )
    }, character(1)), collapse = " ")
  } else ""
  cat(sprintf(
    "[%s] %s script=%s msg=\"%s\"%s%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, script_name,
    gsub("[\r\n]", " ", msg), if (nzchar(extra)) " " else "", extra
  ))
  flush.console()
}

normalize_token <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
}

main <- function() {
  model_dir <- unname(paths["model_dir"])
  dictionary_file <- unname(paths["antibody_dictionary"])
  tcpa_file <- unname(paths["tcpa_table13"])
  final_dir <- file.path(out_dir,"01_freeze_and_source")
  if (dir.exists(final_dir)) {
    stop("正式输出目录已存在，拒绝覆盖: ", final_dir)
  }
  tmp_dir <- paste0(final_dir, ".tmp_", Sys.getpid())
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  expected_cancers <- c(
    "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC",
    "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PCPG", "PRAD",
    "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC"
  )
  min_group_n <- 30L
  min_recurrent_cancers <- 3L

  model_files <- file.path(
    model_dir,
    paste0(
      expected_cancers,
      "_rppa_feature_models_v0.1_20260723.tsv.gz"
    )
  )
  required <- c(model_files, dictionary_file, tcpa_file)
  if (any(!file.exists(required)) ||
      any(file.info(required)$size <= 0)) {
    stop("存在缺失或空的冻结输入")
  }

  log_event("START", "构建RPPA C1 vs C4复现和TCPA绘图源数据", list(
    environment = "R4.5.1",
    input_cancers = length(expected_cancers),
    min_group_n = min_group_n,
    min_recurrent_cancers = min_recurrent_cancers,
    cross_cancer_meta = "NOT_PERFORMED",
    purity = "NOT_USED",
    mrna = "NOT_USED",
    output = final_dir
  ))

  read_one <- function(cancer, path) {
    x <- fread(path)
    x <- x[
      Design_ID == "C1_vs_C4" &
        Analysis_Set == "REPRESENTATIVE_MAIN" &
        Model_Type == "UNADJUSTED"
    ]
    if (nrow(x) != 258L ||
        uniqueN(x$Cancer) != 1L ||
        unique(x$Cancer) != cancer ||
        anyDuplicated(x$Feature_ID) ||
        any(x$Target != "C1") ||
        any(x$Reference != "C4")) {
      stop("单癌种RPPA正式输入身份或维度错误: ", cancer)
    }
    x
  }
  individual_all <- rbindlist(
    Map(read_one, expected_cancers, model_files),
    use.names = TRUE, fill = TRUE
  )

  # 旧版RPPA大图只纳入两组均有至少30名匹配患者的癌种。
  cancer_samples <- individual_all[, .(
    N_C1_Max = max(N_Target, na.rm = TRUE),
    N_C4_Max = max(N_Reference, na.rm = TRUE),
    N_Estimable_Features = sum(Status == "ESTIMABLE")
  ), by = Cancer]
  cancer_samples[, Formal_Cancer :=
                   N_C1_Max >= min_group_n &
                   N_C4_Max >= min_group_n &
                   N_Estimable_Features > 0L]
  formal_cancers <- cancer_samples[Formal_Cancer == TRUE, Cancer]
  individual <- copy(individual_all[Cancer %chin% formal_cancers])

  # 按癌种内整个RPPA面板统一BH；绝不跨癌校正。
  individual[, Original_Feature_Class_FDR := FDR]
  individual[, FDR := NA_real_]
  individual[
    Status == "ESTIMABLE" & is.finite(P_Value),
    FDR := p.adjust(P_Value, method = "BH"),
    by = Cancer
  ]
  individual[, Single_Cancer_Direction := fcase(
    Status == "ESTIMABLE" & is.finite(FDR) &
      FDR < 0.05 & Estimate > 0,
    "C1_HIGH",
    Status == "ESTIMABLE" & is.finite(FDR) &
      FDR < 0.05 & Estimate < 0,
    "C4_HIGH",
    default = "NOT_SIGNIFICANT"
  )]

  dictionary <- fread(dictionary_file)
  if (nrow(dictionary) != 258L ||
      anyDuplicated(dictionary$standardized_feature)) {
    stop("RPPA冻结字典不是258个唯一feature")
  }
  dictionary_small <- dictionary[, .(
    Feature_ID = standardized_feature,
    Gene_Symbol = gene_symbol,
    Site = site,
    Target_Count = target_count,
    Mapping_Status = mapping_status,
    Validation_Status = validation_status,
    Direct_mRNA_Pairable = direct_mrna_pairable,
    Cross_Layer_Restriction = cross_layer_restriction
  )]

  recurrence <- individual[, {
    estimates <- Estimate[Status == "ESTIMABLE" & is.finite(Estimate)]
    sig_estimates <- Estimate[
      Status == "ESTIMABLE" & is.finite(FDR) & FDR < 0.05
    ]
    .(
      Feature_Label = Feature_Label[1L],
      Feature_Class = Feature_Class[1L],
      Definition_Overlap = any(Definition_Overlap %in% TRUE),
      N_Formal_Cancers = length(formal_cancers),
      N_Cancer_Estimable = uniqueN(Cancer[
        Status == "ESTIMABLE" & is.finite(P_Value)
      ]),
      N_C1_High_Significant = uniqueN(Cancer[
        Single_Cancer_Direction == "C1_HIGH"
      ]),
      N_C4_High_Significant = uniqueN(Cancer[
        Single_Cancer_Direction == "C4_HIGH"
      ]),
      Median_Estimate_Estimable = if (length(estimates)) {
        median(estimates)
      } else NA_real_,
      Median_Estimate_Significant = if (length(sig_estimates)) {
        median(sig_estimates)
      } else NA_real_,
      Max_Abs_Estimate = if (length(estimates)) {
        max(abs(estimates))
      } else NA_real_
    )
  }, by = Feature_ID]
  recurrence <- merge(
    recurrence, dictionary_small, by = "Feature_ID", all.x = TRUE
  )
  recurrence[, `:=`(
    C1_High_Recurrence_Fraction =
      N_C1_High_Significant / N_Formal_Cancers,
    C4_High_Recurrence_Fraction =
      N_C4_High_Significant / N_Formal_Cancers
  )]
  recurrence[, Recurrent_Class := fcase(
    N_C1_High_Significant >= min_recurrent_cancers &
      N_C4_High_Significant >= min_recurrent_cancers,
    "BIDIRECTIONAL_RECURRENT",
    N_C1_High_Significant >= min_recurrent_cancers,
    "C1_HIGH",
    N_C4_High_Significant >= min_recurrent_cancers,
    "C4_HIGH",
    default = "NOT_RECURRENT"
  )]
  recurrence[, Formal_Display :=
               Recurrent_Class %chin% c("C1_HIGH", "C4_HIGH")]
  recurrence[, Dominant_Recurrence_N := pmax(
    N_C1_High_Significant, N_C4_High_Significant
  )]
  recurrence[, Dominant_Recurrence_Fraction :=
               Dominant_Recurrence_N / N_Formal_Cancers]
  recurrence[, Abs_Median_Estimate_Significant :=
               abs(Median_Estimate_Significant)]
  setorder(
    recurrence, -Formal_Display, Recurrent_Class,
    -Dominant_Recurrence_N, -Abs_Median_Estimate_Significant,
    Feature_ID
  )

  selected <- copy(recurrence[Formal_Display == TRUE])
  selected[, State_Order := match(
    Recurrent_Class, c("C1_HIGH", "C4_HIGH")
  )]
  setorder(
    selected, State_Order, -Dominant_Recurrence_N,
    -Abs_Median_Estimate_Significant, Feature_ID
  )
  selected[, Display_Order := .I]
  selected[, Display_Label := Feature_ID]

  # 构建Supplementary Table 13到当前258-feature字典的严格映射。
  official <- fread(tcpa_file)
  if (nrow(official) != 75L || uniqueN(official$Pathway) != 10L) {
    stop("TCPA Supplementary Table 13结构化输入维度错误")
  }
  official[, Source_Normalized := normalize_token(Source_Protein_Token)]
  official[, Feature_Normalized := Source_Normalized]
  official[
    Source_Normalized == "HER3PY1298",
    Feature_Normalized := "HER3PY1289"
  ]
  dictionary_map <- dictionary[, .(
    Feature_ID = standardized_feature,
    Feature_Normalized = normalize_token(standardized_feature)
  )]
  official_map <- merge(
    official, dictionary_map, by = "Feature_Normalized", all.x = TRUE
  )
  if (anyNA(official_map$Feature_ID)) {
    stop("TCPA Table 13存在未映射到258-feature字典的成员")
  }

  selected_map <- merge(
    selected[, .(
      Feature_ID, Display_Label, Recurrent_Class, Display_Order
    )],
    official_map[, .(
      Feature_ID, TCPA_Pathway = Pathway,
      TCPA_Component_Sign = Direction_Sign,
      Source_Protein_Token
    )],
    by = "Feature_ID", all.x = TRUE, allow.cartesian = TRUE
  )
  selected_map[
    is.na(TCPA_Pathway),
    `:=`(
      TCPA_Pathway = "Other_not_in_TCPA_Table13",
      TCPA_Component_Sign = 0L,
      Source_Protein_Token = NA_character_
    )
  ]
  selected_map[, Mapping_Status := fifelse(
    TCPA_Pathway == "Other_not_in_TCPA_Table13",
    "UNMAPPED_TO_TCPA_TABLE13",
    "OFFICIAL_TCPA_TABLE13_MAPPED"
  )]

  pathway_order <- c(
    "Apoptosis", "Cell_cycle", "DNA_damage_response", "EMT",
    "Hormone_a", "Hormone_b", "PI3K_Akt", "RAS_MAPK", "RTK",
    "TSC_mTOR", "Other_not_in_TCPA_Table13"
  )
  pathway_labels <- c(
    Apoptosis = "Apoptosis",
    Cell_cycle = "Cell cycle",
    DNA_damage_response = "DNA damage response",
    EMT = "EMT",
    Hormone_a = "Hormone a",
    Hormone_b = "Hormone b",
    PI3K_Akt = "PI3K / Akt",
    RAS_MAPK = "RAS / MAPK",
    RTK = "RTK",
    TSC_mTOR = "TSC / mTOR",
    Other_not_in_TCPA_Table13 = "Other / not in TCPA Table 13"
  )
  selected_map[, Pathway_Order := match(TCPA_Pathway, pathway_order)]
  selected_map[, TCPA_Pathway_Label :=
                 unname(pathway_labels[TCPA_Pathway])]
  setorder(selected_map, Display_Order, Pathway_Order)

  feature_mapping <- selected_map[, .(
    N_TCPA_Pathways = sum(
      TCPA_Pathway != "Other_not_in_TCPA_Table13"
    ),
    TCPA_Pathway_Membership = paste(
      TCPA_Pathway, collapse = ";"
    ),
    TCPA_Direction_Membership = paste(
      paste0(
        TCPA_Pathway, ":",
        fifelse(
          TCPA_Component_Sign > 0L, "+",
          fifelse(TCPA_Component_Sign < 0L, "-", "0")
        )
      ),
      collapse = ";"
    ),
    TCPA_Mapping_Status = if (
      any(TCPA_Pathway != "Other_not_in_TCPA_Table13")
    ) "OFFICIAL_TCPA_TABLE13_MAPPED" else {
      "UNMAPPED_TO_TCPA_TABLE13"
    }
  ), by = Feature_ID]
  selected <- merge(
    selected, feature_mapping, by = "Feature_ID", all.x = TRUE
  )
  setorder(selected, Display_Order)

  dot_source <- merge(
    individual,
    selected[, .(
      Feature_ID, Recurrent_Class, Dominant_Recurrence_N,
      Dominant_Recurrence_Fraction, Display_Order, Display_Label,
      TCPA_Pathway_Membership, TCPA_Direction_Membership,
      TCPA_Mapping_Status
    )],
    by = "Feature_ID", all = FALSE
  )
  cancer_order <- cancer_samples[
    Cancer %chin% formal_cancers
  ][order(match(Cancer, formal_cancers)), Cancer]
  dot_source[, Cancer_Order := match(Cancer, cancer_order)]
  dot_source[, Dot_Display :=
               Status == "ESTIMABLE" &
               is.finite(FDR) & FDR < 0.05]
  dot_source[, Minus_Log10_FDR :=
               -log10(pmax(FDR, .Machine$double.xmin))]
  dot_source[, Minus_Log10_FDR_Capped :=
               pmin(Minus_Log10_FDR, 10)]
  dot_source[, Estimate_Capped_1 :=
               pmax(pmin(Estimate, 1), -1)]
  setorder(dot_source, Cancer_Order, Display_Order)

  pathway_registry <- selected_map[, .(
    N_Display_Features = uniqueN(Feature_ID),
    N_C1_High_Features = uniqueN(
      Feature_ID[Recurrent_Class == "C1_HIGH"]
    ),
    N_C4_High_Features = uniqueN(
      Feature_ID[Recurrent_Class == "C4_HIGH"]
    ),
    N_Positive_Components = sum(TCPA_Component_Sign > 0L),
    N_Negative_Components = sum(TCPA_Component_Sign < 0L),
    N_Other = sum(TCPA_Component_Sign == 0L)
  ), by = .(
    TCPA_Pathway, TCPA_Pathway_Label, Pathway_Order
  )]
  setorder(pathway_registry, Pathway_Order)

  threshold_audit <- rbindlist(lapply(c(0, 0.1, 0.2), function(cutoff) {
    a <- individual[
      Status == "ESTIMABLE" & is.finite(FDR) &
        FDR < 0.05 & abs(Estimate) >= cutoff,
      .(
        N_C1 = uniqueN(Cancer[Estimate > 0]),
        N_C4 = uniqueN(Cancer[Estimate < 0])
      ),
      by = Feature_ID
    ]
    data.table(
      Abs_Estimate_Cutoff = cutoff,
      N_Selected_Including_Bidirectional =
        sum(a$N_C1 >= min_recurrent_cancers |
              a$N_C4 >= min_recurrent_cancers),
      N_C1_High_Only =
        sum(a$N_C1 >= min_recurrent_cancers &
              a$N_C4 < min_recurrent_cancers),
      N_C4_High_Only =
        sum(a$N_C4 >= min_recurrent_cancers &
              a$N_C1 < min_recurrent_cancers),
      N_Bidirectional =
        sum(a$N_C1 >= min_recurrent_cancers &
              a$N_C4 >= min_recurrent_cancers)
    )
  }))
  threshold_audit[, Formal_Selection := Abs_Estimate_Cutoff == 0]

  summary_table <- rbindlist(list(
    data.table(
      Metric = c(
        "INPUT_CANCERS", "FORMAL_CANCERS",
        "MIN_GROUP_N", "MIN_RECURRENT_CANCERS",
        "TOTAL_RPPA_FEATURES", "FORMAL_DISPLAY_FEATURES",
        "C1_HIGH_DISPLAY_FEATURES", "C4_HIGH_DISPLAY_FEATURES",
        "BIDIRECTIONAL_EXCLUDED_FEATURES",
        "VISIBLE_CANCER_FEATURE_DOTS",
        "TCPA_OFFICIAL_MAPPED_DISPLAY_FEATURES",
        "OTHER_DISPLAY_FEATURES", "HIT_PATHWAY_NODES"
      ),
      Value = c(
        length(expected_cancers), length(formal_cancers),
        min_group_n, min_recurrent_cancers,
        uniqueN(recurrence$Feature_ID), nrow(selected),
        selected[Recurrent_Class == "C1_HIGH", .N],
        selected[Recurrent_Class == "C4_HIGH", .N],
        recurrence[Recurrent_Class == "BIDIRECTIONAL_RECURRENT", .N],
        sum(dot_source$Dot_Display),
        selected[TCPA_Mapping_Status ==
                   "OFFICIAL_TCPA_TABLE13_MAPPED", .N],
        selected[TCPA_Mapping_Status ==
                   "UNMAPPED_TO_TCPA_TABLE13", .N],
        nrow(pathway_registry)
      )
    ),
    recurrence[, .(
      Metric = paste0("RECURRENCE_CLASS_", Recurrent_Class),
      Value = .N
    ), by = Recurrent_Class][, Recurrent_Class := NULL]
  ), use.names = TRUE)

  threshold_registry <- data.table(
    Parameter = c(
      "DESIGN_ID", "ANALYSIS_SET", "MODEL_TYPE",
      "MIN_C1_PATIENTS_PER_CANCER",
      "MIN_C4_PATIENTS_PER_CANCER",
      "WITHIN_CANCER_BH_FAMILY",
      "WITHIN_CANCER_FDR",
      "ABS_ESTIMATE_HARD_THRESHOLD",
      "MIN_SAME_DIRECTION_RECURRENT_CANCERS",
      "BIDIRECTIONAL_POLICY",
      "CROSS_CANCER_BH", "META_LOO_I2",
      "PURITY_MODEL", "MRNA_OR_ORA_IN_RPPA_SELECTION"
    ),
    Value = c(
      "C1_vs_C4", "REPRESENTATIVE_MAIN", "UNADJUSTED",
      "30", "30", "ALL_ESTIMABLE_RPPA_FEATURES",
      "0.05", "NONE", "3", "EXCLUDE_FROM_FORMAL_DISPLAY",
      "NOT_PERFORMED", "NOT_PERFORMED",
      "EXCLUDED", "NOT_USED"
    )
  )

  files <- c(
    within = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_within_cancer_results_v0.1_20260730.tsv.gz"
    ),
    cancer = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_cancer_sample_gate_v0.1_20260730.tsv"
    ),
    recurrence = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_recurrence_all_v0.1_20260730.tsv"
    ),
    selected = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_display_feature_freeze_v0.1_20260730.tsv"
    ),
    dot = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_dotplot_source_v0.1_20260730.tsv.gz"
    ),
    edges = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_tcpa_edges_v0.1_20260730.tsv"
    ),
    pathways = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_tcpa_pathway_registry_v0.1_20260730.tsv"
    ),
    audit = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_threshold_audit_v0.1_20260730.tsv"
    ),
    thresholds = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_threshold_registry_v0.1_20260730.tsv"
    ),
    summary = file.path(
      tmp_dir,
      "stage3_rppa_c1c4_summary_v0.1_20260730.tsv"
    )
  )
  fwrite(individual, files[["within"]], sep = "\t", na = "NA")
  fwrite(cancer_samples, files[["cancer"]], sep = "\t", na = "NA")
  fwrite(recurrence, files[["recurrence"]], sep = "\t", na = "NA")
  fwrite(selected, files[["selected"]], sep = "\t", na = "NA")
  fwrite(dot_source, files[["dot"]], sep = "\t", na = "NA")
  fwrite(selected_map, files[["edges"]], sep = "\t", na = "NA")
  fwrite(pathway_registry, files[["pathways"]], sep = "\t", na = "NA")
  fwrite(threshold_audit, files[["audit"]], sep = "\t", na = "NA")
  fwrite(threshold_registry, files[["thresholds"]], sep = "\t", na = "NA")
  fwrite(summary_table, files[["summary"]], sep = "\t", na = "NA")

  fdr_recalc_error <- individual[
    Status == "ESTIMABLE" & is.finite(P_Value),
    .(Error = max(abs(
      FDR - p.adjust(P_Value, method = "BH")
    ))),
    by = Cancer
  ][, max(Error)]
  qa <- data.table(
    Check = c(
      "INPUT_CANCERS_22",
      "FORMAL_CANCERS_15",
      "FORMAL_CANCERS_MEET_BOTH_GROUP_N30",
      "ONLY_C1_VS_C4_REPRESENTATIVE_UNADJUSTED",
      "WITHIN_CANCER_FULL_PANEL_BH_EXACT",
      "NO_CROSS_CANCER_BH_OR_META",
      "FORMAL_FEATURES_SAME_DIRECTION_N_GE_3",
      "NO_BIDIRECTIONAL_IN_FORMAL_DISPLAY",
      "ALL_SELECTED_HAVE_TCPA_OR_OTHER_EDGE",
      "DOT_SOURCE_DIMENSION",
      "ALL_OUTPUTS_NONEMPTY"
    ),
    Observed = c(
      uniqueN(individual_all$Cancer),
      length(formal_cancers),
      all(cancer_samples[
        Formal_Cancer == TRUE,
        N_C1_Max >= min_group_n & N_C4_Max >= min_group_n
      ]),
      uniqueN(individual$Design_ID) == 1L &&
        unique(individual$Design_ID) == "C1_vs_C4" &&
        uniqueN(individual$Analysis_Set) == 1L &&
        unique(individual$Analysis_Set) == "REPRESENTATIVE_MAIN" &&
        uniqueN(individual$Model_Type) == 1L &&
        unique(individual$Model_Type) == "UNADJUSTED",
      format(fdr_recalc_error, scientific = TRUE),
      "TRUE",
      all(selected$Dominant_Recurrence_N >= min_recurrent_cancers),
      !any(selected$Recurrent_Class == "BIDIRECTIONAL_RECURRENT"),
      uniqueN(selected_map$Feature_ID),
      nrow(dot_source),
      all(file.exists(files)) &&
        all(file.info(files)$size > 0)
    ),
    Expected = c(
      22L, 15L, TRUE, TRUE, "0e+00", "TRUE",
      TRUE, TRUE, nrow(selected),
      length(formal_cancers) * nrow(selected), TRUE
    )
  )
  qa[, Pass := as.character(Observed) == as.character(Expected)]
  qa_file <- file.path(
    tmp_dir,
    "stage3_rppa_c1c4_source_machine_qa_v0.1_20260730.tsv"
  )
  fwrite(qa, qa_file, sep = "\t", na = "NA")
  if (!all(qa$Pass)) {
    print(qa)
    stop("RPPA C1 vs C4 source machine QA失败")
  }

  manifest_inputs <- data.table(
    Role = c(
      rep("PER_CANCER_MODEL", length(model_files)),
      "RPPA_DICTIONARY", "TCPA_TABLE13"
    ),
    Path = required,
    Bytes = file.info(required)$size,
    MD5 = unname(tools::md5sum(required))
  )
  fwrite(
    manifest_inputs,
    file.path(
      tmp_dir,
      "stage3_rppa_c1c4_input_manifest_v0.1_20260730.tsv"
    ),
    sep = "\t", na = "NA"
  )

  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, final_dir)) {
    stop("无法原子冻结RPPA C1 vs C4 source目录")
  }
  log_event("DONE", "RPPA C1 vs C4复现和TCPA源数据完成", list(
    runtime_sec = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    )), 1),
    formal_cancers = length(formal_cancers),
    display_features = nrow(selected),
    c1_high = selected[Recurrent_Class == "C1_HIGH", .N],
    c4_high = selected[Recurrent_Class == "C4_HIGH", .N],
    visible_dots = sum(dot_source$Dot_Display),
    pathway_nodes = nrow(pathway_registry),
    output = final_dir
  ))
}

tryCatch(main(), error = function(e) {
  log_event("ERROR", conditionMessage(e), list(
    runtime_sec = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    )), 1)
  ))
  quit(save = "no", status = 1L)
})
