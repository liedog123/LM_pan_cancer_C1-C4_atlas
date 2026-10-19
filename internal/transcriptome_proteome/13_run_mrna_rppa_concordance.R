#!/usr/bin/env Rscript
# Purpose: Summarise reported RNA-protein concordance and discordance.
# Task: mrna_rppa_concordance; step: s13
# Inputs (DCF fields): mrna_model_dir,mrna_recurrence,mrna_display,rppa_within,rppa_recurrence,rppa_display,cancer_gate,antibody_dictionary
# Outputs: new files under output_dir; source inputs are read-only.
# Usage: Rscript run.R --task mrna_rppa_concordance --config /path/to/task.dcf
# Details: docs/reproduction_routes.md, docs/step_index.tsv and docs/software_requirements.md.

usage <- function() cat("Usage: Rscript 13_run_mrna_rppa_concordance.R --config PATH [--check-config]\n       Rscript 13_run_mrna_rppa_concordance.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
config_lines <- readLines(config_path, warn = FALSE)
config_keys <- trimws(sub(":.*$", "", config_lines[
  grepl("^[^[:space:]#][^:]*:", config_lines)]))
if(anyDuplicated(config_keys)) stop("DCF fields must not repeat.")
cfg <- read.dcf(config_path)
fields <- c("mrna_model_dir", "mrna_recurrence", "mrna_display", "rppa_within", "rppa_recurrence", "rppa_display", "cancer_gate", "antibody_dictionary", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
expected_cancers <- c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LGG","LIHC","LUSC","OV","PRAD","SARC","SKCM","STAD","UCEC")
if(!dir.exists(paths["mrna_model_dir"]) || file.access(paths["mrna_model_dir"],5L)!=0L) stop("Model directory must be readable/searchable.")
required_files <- c(paths[setdiff(fields,c("mrna_model_dir","output_dir"))],file.path(paths["mrna_model_dir"],paste0(expected_cancers,"_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz")))

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

script_name <- "13_run_mrna_rppa_concordance.R"
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

collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
  if (length(x)) paste(x, collapse = ";") else ""
}

main <- function() {
  # 独立输出根仅在普通执行创建，源统计与原子提交连续保留。
  mrna_model_dir <- unname(paths["mrna_model_dir"])
  mrna_recurrence_file <- unname(paths["mrna_recurrence"])
  mrna_display_file <- unname(paths["mrna_display"])
  rppa_within_file <- unname(paths["rppa_within"])
  rppa_recurrence_file <- unname(paths["rppa_recurrence"])
  rppa_display_file <- unname(paths["rppa_display"])
  cancer_gate_file <- unname(paths["cancer_gate"])
  dictionary_file <- unname(paths["antibody_dictionary"])
  output_dir <- file.path(out_dir, "01_mrna_rppa_concordance")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(output_dir)) stop("Output subdirectory already exists.")
  tmp_dir <- paste0(output_dir, ".tmp_", Sys.getpid())
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  expected_cancers <- c(
    "BLCA", "BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LGG",
    "LIHC", "LUSC", "OV", "PRAD", "SARC", "SKCM", "STAD", "UCEC"
  )
  mrna_model_files <- file.path(
    mrna_model_dir,
    paste0(
      expected_cancers,
      "_mrna_ovr_c1c4_models_v0.2_20260729.tsv.gz"
    )
  )
  required <- c(
    mrna_model_files, mrna_recurrence_file, mrna_display_file,
    rppa_within_file, rppa_recurrence_file, rppa_display_file,
    cancer_gate_file, dictionary_file
  )
  if (any(!file.exists(required)) ||
      any(file.info(required)$size <= 0)) {
    stop("存在缺失或空的跨层审计冻结输入")
  }

  log_event("START", "开始mRNA–RPPA C1 vs C4一致性/矛盾审计", list(
    environment = "R4.5.1",
    cancers = length(expected_cancers),
    design = "C1_vs_C4",
    mrna_threshold = "FDR<0.05_AND_ABS_ESTIMATE>1.0",
    rppa_threshold = "FDR<0.05_NO_ABS_EFFECT_CUTOFF",
    cross_cancer_bh = "NOT_PERFORMED",
    meta_loo_i2 = "NOT_PERFORMED",
    output = output_dir
  ))

  cancer_gate <- fread(cancer_gate_file)
  formal_cancers <- cancer_gate[Formal_Cancer == TRUE, Cancer]
  if (!setequal(formal_cancers, expected_cancers) ||
      length(formal_cancers) != 15L) {
    stop("RPPA正式15癌种集合与合同不一致")
  }

  dictionary <- fread(dictionary_file)
  if (nrow(dictionary) != 258L ||
      anyDuplicated(dictionary$standardized_feature)) {
    stop("RPPA冻结字典不是258个唯一feature")
  }
  dictionary[, Gene := toupper(trimws(as.character(gene_symbol)))]
  dictionary[, Base_Direct_Eligible :=
               feature_class == "TOTAL" &
               direct_mrna_pairable == TRUE &
               target_count == 1L &
               !is.na(Gene) & nzchar(Gene)]
  dictionary[
    Base_Direct_Eligible == TRUE,
    N_Direct_Total_Features_Per_Gene := .N,
    by = Gene
  ]
  dictionary[
    is.na(N_Direct_Total_Features_Per_Gene),
    N_Direct_Total_Features_Per_Gene := 0L
  ]
  dictionary[, Cross_Layer_Route := fcase(
    feature_class == "PHOSPHO",
    "EXECUTION_LAYER_CONTEXT_PHOSPHO",
    feature_class == "CLEAVED",
    "EXECUTION_LAYER_CONTEXT_CLEAVED",
    Base_Direct_Eligible == TRUE &
      N_Direct_Total_Features_Per_Gene == 1L,
    "DIRECT_TOTAL_UNIQUE_GENE",
    Base_Direct_Eligible == TRUE &
      N_Direct_Total_Features_Per_Gene > 1L,
    "UNEVALUABLE_DUPLICATE_TOTAL_GENE",
    feature_class == "TOTAL",
    "UNEVALUABLE_TOTAL_MAPPING",
    default = "UNEVALUABLE_OTHER"
  )]

  rppa_recurrence <- fread(rppa_recurrence_file)
  rppa_recurrence[, Gene := toupper(trimws(as.character(Gene_Symbol)))]
  rppa_display <- fread(rppa_display_file)
  rppa_display[, Gene := toupper(trimws(as.character(Gene_Symbol)))]
  rppa_display_small <- rppa_display[, .(
    Feature_ID,
    RPPA_Formal_Display = TRUE,
    RPPA_Display_Order = Display_Order,
    RPPA_TCPA_Pathway_Membership = TCPA_Pathway_Membership,
    RPPA_TCPA_Mapping_Status = TCPA_Mapping_Status
  )]
  feature_routing <- merge(
    dictionary[, .(
      Feature_ID = standardized_feature,
      Feature_Label = raw_feature,
      Feature_Class = feature_class,
      Gene,
      Gene_Symbol = gene_symbol,
      Site = site,
      Target_Count = target_count,
      Direct_mRNA_Pairable = direct_mrna_pairable,
      Dictionary_Mapping_Status = mapping_status,
      Validation_Status = validation_status,
      Dictionary_Cross_Layer_Restriction = cross_layer_restriction,
      N_Direct_Total_Features_Per_Gene,
      Cross_Layer_Route
    )],
    rppa_recurrence[, .(
      Feature_ID, Definition_Overlap,
      RPPA_N_Cancer_Estimable = N_Cancer_Estimable,
      RPPA_N_C1_High_Significant = N_C1_High_Significant,
      RPPA_N_C4_High_Significant = N_C4_High_Significant,
      RPPA_Recurrent_Class = Recurrent_Class,
      RPPA_Formal_Display = Formal_Display
    )],
    by = "Feature_ID", all.x = TRUE
  )
  feature_routing <- merge(
    feature_routing, rppa_display_small,
    by = "Feature_ID", all.x = TRUE,
    suffixes = c("", "_Display_Table")
  )
  feature_routing[
    is.na(RPPA_Formal_Display),
    RPPA_Formal_Display := FALSE
  ]
  feature_routing[
    is.na(RPPA_Formal_Display_Display_Table),
    RPPA_Formal_Display_Display_Table := FALSE
  ]
  setorder(feature_routing, Cross_Layer_Route, Gene, Feature_ID)

  pair_registry <- feature_routing[
    Cross_Layer_Route == "DIRECT_TOTAL_UNIQUE_GENE",
    .(
      Gene, Feature_ID,
      RPPA_Feature_ID = Feature_ID,
      Feature_Label, Feature_Class,
      Gene_Symbol, Definition_Overlap,
      Dictionary_Mapping_Status, Validation_Status,
      RPPA_N_Cancer_Estimable,
      RPPA_N_C1_High_Significant,
      RPPA_N_C4_High_Significant,
      RPPA_Recurrent_Class,
      RPPA_Formal_Display,
      RPPA_Display_Order,
      RPPA_TCPA_Pathway_Membership,
      RPPA_TCPA_Mapping_Status
    )
  ]
  if (anyDuplicated(pair_registry$Gene) ||
      anyDuplicated(pair_registry$Feature_ID)) {
    stop("直接total protein–mRNA配对不是gene与feature双向唯一")
  }

  mrna_keep <- c(
    "Cancer", "Feature_ID", "Feature_Label", "Feature_Class",
    "Definition_Overlap", "Design_ID", "Analysis_Set", "Model_Type",
    "Target", "Reference", "Estimate", "SE", "P_Value", "FDR",
    "Status", "Failure_Reason"
  )
  mrna_rows <- rbindlist(lapply(
    seq_along(mrna_model_files),
    function(i) {
      x <- fread(mrna_model_files[i], select = mrna_keep)
      x <- x[
        Design_ID == "C1_vs_C4" &
          Analysis_Set == "REPRESENTATIVE_MAIN" &
          Model_Type == "UNADJUSTED"
      ]
      x[, Gene := toupper(trimws(as.character(Feature_ID)))]
      x <- x[Gene %chin% pair_registry$Gene]
      if (anyDuplicated(x$Gene)) {
        stop("单癌种mRNA直接配对gene重复: ", expected_cancers[i])
      }
      x
    }
  ), use.names = TRUE, fill = TRUE)

  rppa_within <- fread(rppa_within_file)
  rppa_within <- rppa_within[
    Cancer %chin% expected_cancers &
      Feature_ID %chin% pair_registry$Feature_ID
  ]
  rppa_within <- merge(
    rppa_within,
    pair_registry[, .(Gene, Feature_ID)],
    by = "Feature_ID", all.x = TRUE
  )
  if (anyDuplicated(rppa_within[, .(Cancer, Gene)])) {
    stop("RPPA直接配对Cancer × Gene重复")
  }

  mrna_pair <- mrna_rows[, .(
    Cancer, Gene,
    mRNA_Feature_ID = Feature_ID,
    mRNA_Definition_Overlap = Definition_Overlap,
    mRNA_Estimate = Estimate,
    mRNA_SE = SE,
    mRNA_P_Value = P_Value,
    mRNA_FDR = FDR,
    mRNA_Status = Status,
    mRNA_Failure_Reason = Failure_Reason
  )]
  rppa_pair <- rppa_within[, .(
    Cancer, Gene,
    RPPA_Feature_ID = Feature_ID,
    RPPA_Feature_Label = Feature_Label,
    RPPA_Definition_Overlap = Definition_Overlap,
    RPPA_Estimate = Estimate,
    RPPA_SE = SE,
    RPPA_P_Value = P_Value,
    RPPA_FDR = FDR,
    RPPA_Status = Status,
    RPPA_Failure_Reason = Failure_Reason
  )]
  cancer_gene <- CJ(
    Cancer = expected_cancers,
    Gene = pair_registry$Gene,
    unique = TRUE
  )
  cancer_gene <- merge(
    cancer_gene, mrna_pair,
    by = c("Cancer", "Gene"), all.x = TRUE
  )
  cancer_gene <- merge(
    cancer_gene, rppa_pair,
    by = c("Cancer", "Gene"), all.x = TRUE
  )
  cancer_gene <- merge(
    cancer_gene,
    pair_registry[, .(
      Gene, RPPA_Feature_ID_Registry = Feature_ID,
      RPPA_Gene_Definition_Overlap = Definition_Overlap
    )],
    by = "Gene", all.x = TRUE
  )

  cancer_gene[, `:=`(
    mRNA_Evaluable =
      mRNA_Status == "ESTIMABLE" &
      is.finite(mRNA_Estimate) &
      is.finite(mRNA_FDR),
    RPPA_Evaluable =
      RPPA_Status == "ESTIMABLE" &
      is.finite(RPPA_Estimate) &
      is.finite(RPPA_FDR)
  )]
  cancer_gene[, mRNA_Sig_Direction := fcase(
    mRNA_Evaluable & mRNA_FDR < 0.05 & mRNA_Estimate > 1.0,
    "C1_HIGH",
    mRNA_Evaluable & mRNA_FDR < 0.05 & mRNA_Estimate < -1.0,
    "C4_HIGH",
    default = "NOT_SIGNIFICANT"
  )]
  cancer_gene[, RPPA_Sig_Direction := fcase(
    RPPA_Evaluable & RPPA_FDR < 0.05 & RPPA_Estimate > 0,
    "C1_HIGH",
    RPPA_Evaluable & RPPA_FDR < 0.05 & RPPA_Estimate < 0,
    "C4_HIGH",
    default = "NOT_SIGNIFICANT"
  )]
  cancer_gene[, Cancer_Level_Audit_Class := fcase(
    !mRNA_Evaluable | !RPPA_Evaluable,
    "UNEVALUABLE_COVERAGE",
    mRNA_Sig_Direction != "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction != "NOT_SIGNIFICANT" &
      mRNA_Sig_Direction == RPPA_Sig_Direction,
    "CONCORDANT_SIGNIFICANT",
    mRNA_Sig_Direction != "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction != "NOT_SIGNIFICANT" &
      mRNA_Sig_Direction != RPPA_Sig_Direction,
    "DISCORDANT_SIGNIFICANT",
    mRNA_Sig_Direction != "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction == "NOT_SIGNIFICANT",
    "RNA_ONLY_SIGNIFICANT",
    mRNA_Sig_Direction == "NOT_SIGNIFICANT" &
      RPPA_Sig_Direction != "NOT_SIGNIFICANT",
    "PROTEIN_ONLY_SIGNIFICANT",
    default = "BOTH_NONSIGNIFICANT"
  )]
  cancer_gene[, Effect_Sign_Relation := fcase(
    !mRNA_Evaluable | !RPPA_Evaluable,
    "UNEVALUABLE",
    sign(mRNA_Estimate) == sign(RPPA_Estimate),
    "SAME_DIRECTION",
    default = "OPPOSITE_DIRECTION"
  )]
  cancer_gene[, Cancer_Order := match(Cancer, expected_cancers)]
  setorder(cancer_gene, Gene, Cancer_Order)

  mrna_recurrence <- fread(mrna_recurrence_file)
  mrna_recurrence <- mrna_recurrence[
    Design_ID == "C1_vs_C4"
  ]
  mrna_recurrence[, Gene := toupper(trimws(as.character(Gene)))]
  mrna_recurrence_small <- mrna_recurrence[, .(
    Gene,
    mRNA_Recurrence_Class = Recurrence_Class,
    mRNA_N_Evaluable_Cancers = N_Evaluable_Cancers,
    mRNA_N_C1_High_Significant = N_Sig_Up,
    mRNA_N_C4_High_Significant = N_Sig_Down,
    mRNA_Recurrent_C1_High = Recurrent_Up,
    mRNA_Recurrent_C4_High = Recurrent_Down,
    mRNA_Bidirectional_Recurrent = Bidirectional_Recurrent,
    mRNA_Definition_Overlap_Recurrence = Definition_Overlap
  )]
  rppa_recurrence_small <- rppa_recurrence[, .(
    Gene,
    RPPA_Feature_ID = Feature_ID,
    RPPA_Recurrence_Class = Recurrent_Class,
    RPPA_N_Evaluable_Cancers = N_Cancer_Estimable,
    RPPA_Definition_Overlap_Recurrence = Definition_Overlap
  )]
  gene_summary <- merge(
    pair_registry,
    mrna_recurrence_small,
    by = "Gene", all.x = TRUE
  )
  gene_summary <- merge(
    gene_summary,
    rppa_recurrence_small,
    by = c("Gene", "RPPA_Feature_ID"),
    all.x = TRUE
  )

  per_gene_counts <- cancer_gene[, .(
    N_Formal_Cancers = .N,
    N_Both_Evaluable = sum(mRNA_Evaluable & RPPA_Evaluable),
    N_Concordant_Significant =
      sum(Cancer_Level_Audit_Class == "CONCORDANT_SIGNIFICANT"),
    N_Discordant_Significant =
      sum(Cancer_Level_Audit_Class == "DISCORDANT_SIGNIFICANT"),
    N_RNA_Only_Significant =
      sum(Cancer_Level_Audit_Class == "RNA_ONLY_SIGNIFICANT"),
    N_Protein_Only_Significant =
      sum(Cancer_Level_Audit_Class == "PROTEIN_ONLY_SIGNIFICANT"),
    N_Both_Nonsignificant =
      sum(Cancer_Level_Audit_Class == "BOTH_NONSIGNIFICANT"),
    N_Unevaluable =
      sum(Cancer_Level_Audit_Class == "UNEVALUABLE_COVERAGE"),
    Concordant_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "CONCORDANT_SIGNIFICANT"]
    ),
    Discordant_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "DISCORDANT_SIGNIFICANT"]
    ),
    RNA_Only_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "RNA_ONLY_SIGNIFICANT"]
    ),
    Protein_Only_Cancers = collapse_values(
      Cancer[Cancer_Level_Audit_Class == "PROTEIN_ONLY_SIGNIFICANT"]
    )
  ), by = Gene]
  gene_summary <- merge(
    gene_summary, per_gene_counts, by = "Gene", all.x = TRUE
  )

  gene_summary[, mRNA_Formal_Direction := fcase(
    mRNA_Recurrence_Class == "RECURRENT_UP", "C1_HIGH",
    mRNA_Recurrence_Class == "RECURRENT_DOWN", "C4_HIGH",
    default = "NO_UNIQUE_FORMAL_DIRECTION"
  )]
  gene_summary[, RPPA_Formal_Direction := fcase(
    RPPA_Recurrence_Class == "C1_HIGH", "C1_HIGH",
    RPPA_Recurrence_Class == "C4_HIGH", "C4_HIGH",
    default = "NO_UNIQUE_FORMAL_DIRECTION"
  )]
  gene_summary[, `:=`(
    mRNA_Formal_Supported =
      mRNA_Formal_Direction %chin% c("C1_HIGH", "C4_HIGH"),
    RPPA_Formal_Supported =
      RPPA_Formal_Direction %chin% c("C1_HIGH", "C4_HIGH")
  )]
  gene_summary[, Cross_Layer_Support_Class := fcase(
    mRNA_Recurrence_Class == "BIDIRECTIONAL_RECURRENT" |
      RPPA_Recurrence_Class == "BIDIRECTIONAL_RECURRENT",
    "BIDIRECTIONAL_OR_UNRESOLVED",
    mRNA_Formal_Supported & RPPA_Formal_Supported &
      mRNA_Formal_Direction == RPPA_Formal_Direction,
    "CONCORDANT_SUPPORTED",
    mRNA_Formal_Supported & RPPA_Formal_Supported &
      mRNA_Formal_Direction != RPPA_Formal_Direction,
    "DISCORDANT_SUPPORTED",
    mRNA_Formal_Supported & !RPPA_Formal_Supported,
    "RNA_ONLY_SUPPORTED",
    !mRNA_Formal_Supported & RPPA_Formal_Supported,
    "PROTEIN_ONLY_SUPPORTED",
    default = "NO_CROSS_CANCER_SUPPORT"
  )]

  mrna_display <- fread(mrna_display_file)
  mrna_display[, Gene := toupper(trimws(as.character(Gene)))]
  mrna_display_small <- mrna_display[, .(
    mRNA_Selected_For_Cytoscape = any(Selected_For_Cytoscape),
    mRNA_Display_State = collapse_values(State),
    mRNA_Display_Direction = collapse_values(Direction),
    mRNA_Primary_Cluster_ID = collapse_values(Primary_Cluster_ID),
    mRNA_Display_Assignment = collapse_values(Display_Assignment)
  ), by = Gene]
  gene_summary <- merge(
    gene_summary, mrna_display_small,
    by = "Gene", all.x = TRUE
  )
  gene_summary[
    is.na(mRNA_Selected_For_Cytoscape),
    mRNA_Selected_For_Cytoscape := FALSE
  ]
  gene_summary[, Definition_Overlap_Any :=
                 Definition_Overlap %in% TRUE |
                 mRNA_Definition_Overlap_Recurrence %in% TRUE |
                 RPPA_Definition_Overlap_Recurrence %in% TRUE]
  gene_summary[, Interpretation_Boundary := fifelse(
    Definition_Overlap_Any,
    "DEFINITION_OVERLAP_AXIS_CONSISTENCY_ONLY",
    "ORTHOGONAL_COMPATIBILITY_NOT_CAUSAL"
  )]
  setorder(
    gene_summary, Cross_Layer_Support_Class,
    -N_Concordant_Significant, -N_Discordant_Significant, Gene
  )

  execution_context <- feature_routing[
    Cross_Layer_Route %chin% c(
      "EXECUTION_LAYER_CONTEXT_PHOSPHO",
      "EXECUTION_LAYER_CONTEXT_CLEAVED"
    )
  ]
  execution_context[, Interpretation_Boundary :=
                      "EXECUTION_LAYER_CONTEXT_NO_REQUIRED_MRNA_DIRECTION"]

  cancer_class_summary <- cancer_gene[, .(
    N_Cancer_Gene_Cells = .N,
    N_Genes = uniqueN(Gene),
    N_Cancers = uniqueN(Cancer)
  ), by = Cancer_Level_Audit_Class]
  support_class_summary <- gene_summary[, .(
    N_Genes = .N,
    N_mRNA_Display_Genes = sum(mRNA_Selected_For_Cytoscape),
    N_RPPA_Display_Features = sum(RPPA_Formal_Display %in% TRUE),
    N_Definition_Overlap = sum(Definition_Overlap_Any)
  ), by = Cross_Layer_Support_Class]
  display_overlap <- gene_summary[
    mRNA_Selected_For_Cytoscape == TRUE |
      RPPA_Formal_Display == TRUE,
    .(
      Gene, RPPA_Feature_ID,
      Cross_Layer_Support_Class,
      mRNA_Formal_Direction, RPPA_Formal_Direction,
      N_Concordant_Significant, N_Discordant_Significant,
      mRNA_Selected_For_Cytoscape,
      mRNA_Display_State, mRNA_Primary_Cluster_ID,
      RPPA_Formal_Display,
      RPPA_TCPA_Pathway_Membership,
      Definition_Overlap_Any,
      Interpretation_Boundary
    )
  ]

  threshold_registry <- data.table(
    Parameter = c(
      "DESIGN_ID", "FORMAL_CANCERS",
      "DIRECT_MAPPING_FEATURE_CLASS",
      "DIRECT_MAPPING_TARGET_COUNT",
      "DIRECT_MAPPING_GENE_MULTIPLICITY",
      "MRNA_WITHIN_CANCER_FDR",
      "MRNA_ABS_EFFECT",
      "MRNA_RECURRENCE",
      "RPPA_WITHIN_CANCER_FDR",
      "RPPA_ABS_EFFECT",
      "RPPA_RECURRENCE",
      "PHOSPHO_CLEAVED_POLICY",
      "CROSS_CANCER_BH", "META_LOO_I2", "PURITY"
    ),
    Value = c(
      "C1_vs_C4", "15",
      "TOTAL_ONLY", "1", "EXACTLY_ONE_TOTAL_FEATURE_PER_GENE",
      "0.05", ">1.0", "AT_LEAST_30_PERCENT_EVALUABLE_CANCERS",
      "0.05", "NONE", "AT_LEAST_3_OF_15_CANCERS",
      "EXECUTION_LAYER_CONTEXT_ONLY",
      "NOT_PERFORMED", "NOT_PERFORMED", "NOT_USED"
    )
  )

  files <- c(
    routing = file.path(
      tmp_dir,
      "stage3_mrna_rppa_feature_routing_v0.1_20260730.tsv"
    ),
    pairs = file.path(
      tmp_dir,
      "stage3_mrna_rppa_unique_total_pair_registry_v0.1_20260730.tsv"
    ),
    cells = file.path(
      tmp_dir,
      "stage3_mrna_rppa_cancer_gene_audit_v0.1_20260730.tsv.gz"
    ),
    genes = file.path(
      tmp_dir,
      "stage3_mrna_rppa_gene_summary_v0.1_20260730.tsv"
    ),
    execution = file.path(
      tmp_dir,
      "stage3_mrna_rppa_execution_layer_context_v0.1_20260730.tsv"
    ),
    cancer_summary = file.path(
      tmp_dir,
      "stage3_mrna_rppa_cancer_class_summary_v0.1_20260730.tsv"
    ),
    support_summary = file.path(
      tmp_dir,
      "stage3_mrna_rppa_support_class_summary_v0.1_20260730.tsv"
    ),
    display = file.path(
      tmp_dir,
      "stage3_mrna_rppa_display_overlap_audit_v0.1_20260730.tsv"
    ),
    thresholds = file.path(
      tmp_dir,
      "stage3_mrna_rppa_threshold_registry_v0.1_20260730.tsv"
    )
  )
  fwrite(feature_routing, files[["routing"]], sep = "\t", na = "NA")
  fwrite(pair_registry, files[["pairs"]], sep = "\t", na = "NA")
  fwrite(cancer_gene, files[["cells"]], sep = "\t", na = "NA")
  fwrite(gene_summary, files[["genes"]], sep = "\t", na = "NA")
  fwrite(execution_context, files[["execution"]], sep = "\t", na = "NA")
  fwrite(
    cancer_class_summary,
    files[["cancer_summary"]], sep = "\t", na = "NA"
  )
  fwrite(
    support_class_summary,
    files[["support_summary"]], sep = "\t", na = "NA"
  )
  fwrite(display_overlap, files[["display"]], sep = "\t", na = "NA")
  fwrite(
    threshold_registry,
    files[["thresholds"]], sep = "\t", na = "NA"
  )

  expected_cell_classes <- c(
    "CONCORDANT_SIGNIFICANT", "DISCORDANT_SIGNIFICANT",
    "RNA_ONLY_SIGNIFICANT", "PROTEIN_ONLY_SIGNIFICANT",
    "BOTH_NONSIGNIFICANT", "UNEVALUABLE_COVERAGE"
  )
  expected_support_classes <- c(
    "CONCORDANT_SUPPORTED", "DISCORDANT_SUPPORTED",
    "RNA_ONLY_SUPPORTED", "PROTEIN_ONLY_SUPPORTED",
    "NO_CROSS_CANCER_SUPPORT", "BIDIRECTIONAL_OR_UNRESOLVED"
  )
  qa <- data.table(
    Check = c(
      "FORMAL_CANCERS_15",
      "RPPA_FEATURES_258_ROUTED",
      "DIRECT_PAIR_GENES_187",
      "DIRECT_PAIR_GENE_FEATURE_UNIQUE",
      "NO_PHOSPHO_CLEAVED_IN_DIRECT_PAIRS",
      "CANCER_GENE_CELLS_15_BY_187",
      "CANCER_LEVEL_CLASSES_EXHAUSTIVE",
      "GENE_SUPPORT_CLASSES_EXHAUSTIVE",
      "ALL_DIRECT_GENES_HAVE_SUMMARY",
      "EXECUTION_LAYER_FEATURES_60",
      "NO_CROSS_CANCER_BH_META_FIELDS",
      "ALL_OUTPUTS_NONEMPTY"
    ),
    Observed = c(
      length(formal_cancers),
      nrow(feature_routing),
      nrow(pair_registry),
      !anyDuplicated(pair_registry$Gene) &&
        !anyDuplicated(pair_registry$Feature_ID),
      all(pair_registry$Feature_Class == "TOTAL"),
      nrow(cancer_gene),
      all(cancer_gene$Cancer_Level_Audit_Class %chin%
            expected_cell_classes),
      all(gene_summary$Cross_Layer_Support_Class %chin%
            expected_support_classes),
      nrow(gene_summary),
      nrow(execution_context),
      !any(grepl(
        "META|LOO|I2|CROSS_CANCER_FDR",
        names(gene_summary), ignore.case = TRUE
      )),
      all(file.exists(files)) &&
        all(file.info(files)$size > 0)
    ),
    Expected = c(
      15L, 258L, 187L, TRUE, TRUE,
      15L * 187L, TRUE, TRUE, 187L, 60L, TRUE, TRUE
    )
  )
  qa[, Pass := as.character(Observed) == as.character(Expected)]
  qa_file <- file.path(
    tmp_dir,
    "stage3_mrna_rppa_machine_qa_v0.1_20260730.tsv"
  )
  fwrite(qa, qa_file, sep = "\t", na = "NA")
  if (!all(qa$Pass)) {
    print(qa)
    stop("mRNA–RPPA一致性/矛盾审计machine QA失败")
  }

  input_manifest <- data.table(
    Role = c(
      rep("MRNA_CANCER_MODEL", length(mrna_model_files)),
      "MRNA_RECURRENCE", "MRNA_DISPLAY",
      "RPPA_WITHIN_CANCER", "RPPA_RECURRENCE", "RPPA_DISPLAY",
      "RPPA_CANCER_GATE", "RPPA_DICTIONARY"
    ),
    Path = required,
    Bytes = file.info(required)$size,
    MD5 = unname(tools::md5sum(required))
  )
  fwrite(
    input_manifest,
    file.path(
      tmp_dir,
      "stage3_mrna_rppa_input_manifest_v0.1_20260730.tsv"
    ),
    sep = "\t", na = "NA"
  )

  dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp_dir, output_dir)) {
    stop("无法原子冻结mRNA–RPPA审计目录")
  }
  log_event("DONE", "mRNA–RPPA C1 vs C4一致性/矛盾审计完成", list(
    runtime_sec = round(as.numeric(difftime(
      Sys.time(), start_time, units = "secs"
    )), 1),
    direct_pair_genes = nrow(pair_registry),
    cancer_gene_cells = nrow(cancer_gene),
    concordant_supported = gene_summary[
      Cross_Layer_Support_Class == "CONCORDANT_SUPPORTED", .N
    ],
    discordant_supported = gene_summary[
      Cross_Layer_Support_Class == "DISCORDANT_SUPPORTED", .N
    ],
    execution_layer_features = nrow(execution_context),
    output = output_dir
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
