#!/usr/bin/env Rscript
# 路径元信息检查先于研究输入读取与分析依赖加载。
usage <- function() cat("Usage: Rscript 16_prepare_mrna_representatives.R --config PATH [--check-config]\n       Rscript 16_prepare_mrna_representatives.R --help\nDCF-relative paths; metadata checks create no outputs.\n")
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
fields <- c("ora_significant", "recurrent_directional_genes", "significant_term_gene_membership", "upstream_qa", "output_dir")
if(nrow(cfg)!=1L || !setequal(colnames(cfg),fields) || ncol(cfg)!=length(fields)) stop("DCF must contain exactly the documented fields.")
cfg <- setNames(as.character(cfg[1L,fields]),fields)
if(anyNA(cfg) || any(!nzchar(trimws(cfg))) || any(grepl("SET_ME",cfg,fixed=TRUE))) stop("Replace placeholders.")
path_fields <- fields
paths <- vapply(cfg[path_fields],resolve_path,character(1L),base_dir=dirname(config_path))
out_dir <- unname(paths["output_dir"])
input_paths <- paths[names(paths)!="output_dir"]
required_files <- unname(input_paths)

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


options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(data.table))
# 固定人工代表项来自source37；不重新挑选term或重算ORA/BH。
cluster_registry <- data.table(
  Direction = c(rep("UP", 8), rep("DOWN", 2)),
  State = c(rep("C1", 8), rep("C4", 2)),
  Cluster_ID = c(
    sprintf("C1_FC%02d", 1:8),
    sprintf("C4_FC%02d", 1:2)
  ),
  Cluster_Order = 1:10,
  Cluster_Label = c(
    "Cell cycle and chromosome organization",
    "Hypoxia and glycolytic metabolic adaptation",
    "ECM, integrin and cell adhesion",
    "EMT, motility and tissue morphogenesis",
    "Growth-factor, RTK and mTOR/TGF-beta signaling",
    "Angiogenesis, vascular remodeling and hemostasis",
    "Inflammatory, cytokine and IFN-associated remodeling",
    "Protein secretion, Golgi-ER trafficking and stress response",
    "Aerobic respiration and oxidative energy metabolism",
    "Electron transport, oxidative phosphorylation and ATP synthesis"
  ),
  Parent_Domain = c(
    rep(NA_character_, 8),
    rep("Mitochondrial respiratory metabolism", 2)
  ),
  Cluster_Color = c(
    "#CC79A7", "#E69F00", "#009E73", "#D55E00",
    "#7B61A8", "#56B4E9", "#B8860B", "#7A7A7A",
    "#0072B2", "#4C78A8"
  ),
  Semantic_Boundary = c(
    "Cell-cycle-associated bulk expression; not a proliferative dependency.",
    "Hypoxia/glycolysis-associated enrichment; not proof that lactate or low pH drives C1.",
    "Bulk ECM/adhesion association may include stromal composition; not direct matrix causality.",
    "Bulk EMT/motility annotation; not proof of tumor-cell EMT or invasion.",
    "Growth-factor pathway-associated expression; not evidence of pathway addiction.",
    "Vascular/hemostatic association; not proof of angiogenic dependence.",
    "Inflammatory/IFN-associated expression; not direct immune-cell activation or causal signaling.",
    "Secretory/trafficking/stress-associated expression; not proof of secretory dependency.",
    "C4-high mitochondrial respiration annotation supported by a small recurrent gene set.",
    "C4-high electron-transport annotation; related to, not independent from, C4_FC01."
  )
)

reg_term <- function(
    direction, cluster_id, term_key, display_label, term_order,
    semantic_boundary = NA_character_) {
  data.table(
    Direction = direction,
    Cluster_ID = cluster_id,
    Term_Key = term_key,
    Display_Term_Label = display_label,
    Term_Order = term_order,
    Manual_Review_Status = "FROZEN_SAFE_LABEL",
    Semantic_Boundary = semantic_boundary
  )
}

# 显式人工注册表：只选已核验的非C7正式显著term。
term_registry <- rbindlist(list(
  reg_term("UP", "C1_FC01", "HALLMARK::HALLMARK_G2M_CHECKPOINT",
           "G2/M checkpoint-associated expression", 1),
  reg_term("UP", "C1_FC01", "HALLMARK::HALLMARK_MITOTIC_SPINDLE",
           "Mitotic spindle organization", 2),
  reg_term("UP", "C1_FC01", "GO_BP::GOBP_SISTER_CHROMATID_SEGREGATION",
           "Sister chromatid segregation", 3),
  reg_term("UP", "C1_FC02", "HALLMARK::HALLMARK_HYPOXIA",
           "Hypoxia-associated expression", 1),
  reg_term("UP", "C1_FC02", "HALLMARK::HALLMARK_GLYCOLYSIS",
           "Glycolytic metabolism", 2),
  reg_term("UP", "C1_FC02", "GO_BP::GOBP_RESPONSE_TO_OXYGEN_LEVELS",
           "Response to oxygen levels", 3),
  reg_term("UP", "C1_FC02", "GO_BP::GOBP_PYRUVATE_METABOLIC_PROCESS",
           "Pyruvate metabolism", 4),
  reg_term("UP", "C1_FC03", "REACTOME::REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
           "Extracellular matrix organization", 1),
  reg_term("UP", "C1_FC03", "REACTOME::REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
           "Integrin cell-surface interactions", 2),
  reg_term("UP", "C1_FC03", "REACTOME::REACTOME_COLLAGEN_FORMATION",
           "Collagen formation", 3),
  reg_term("UP", "C1_FC04", "HALLMARK::HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "EMT-associated expression", 1),
  reg_term("UP", "C1_FC04", "GO_BP::GOBP_WOUND_HEALING",
           "Wound-healing-associated program", 2),
  reg_term("UP", "C1_FC04", "GO_BP::GOBP_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
           "Regulation of epithelial cell migration", 3),
  reg_term("UP", "C1_FC05", "HALLMARK::HALLMARK_MTORC1_SIGNALING",
           "mTORC1-signaling-associated expression", 1),
  reg_term("UP", "C1_FC05", "HALLMARK::HALLMARK_TGF_BETA_SIGNALING",
           "TGF-beta-signaling-associated expression", 2),
  reg_term("UP", "C1_FC05", "GO_BP::GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
           "Response to TGF-beta", 3),
  reg_term(
    "UP", "C1_FC05",
    "GO_BP::GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS",
    "Regulation of growth-factor response", 4
  ),
  reg_term("UP", "C1_FC06", "GO_BP::GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT",
           "Regulation of vasculature development", 1),
  reg_term("UP", "C1_FC06", "HALLMARK::HALLMARK_ANGIOGENESIS",
           "Angiogenesis-associated expression", 2),
  reg_term(
    "UP", "C1_FC06",
    "REACTOME::REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
    "Platelet activation and aggregation annotation", 3
  ),
  reg_term("UP", "C1_FC07", "HALLMARK::HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "TNF-NF-kappaB-associated expression", 1),
  reg_term("UP", "C1_FC07", "HALLMARK::HALLMARK_INFLAMMATORY_RESPONSE",
           "Inflammatory-response-associated expression", 2),
  reg_term("UP", "C1_FC07", "REACTOME::REACTOME_INTERFERON_SIGNALING",
           "Interferon-signaling-associated expression", 3),
  reg_term("UP", "C1_FC07", "GO_BP::GOBP_MYELOID_LEUKOCYTE_DIFFERENTIATION",
           "Myeloid leukocyte differentiation", 4),
  reg_term("UP", "C1_FC08", "HALLMARK::HALLMARK_PROTEIN_SECRETION",
           "Protein secretion", 1),
  reg_term(
    "UP", "C1_FC08",
    "REACTOME::REACTOME_INTRA_GOLGI_AND_RETROGRADE_GOLGI_TO_ER_TRAFFIC",
    "Intra-Golgi and Golgi-to-ER trafficking", 2
  ),
  reg_term("UP", "C1_FC08", "HALLMARK::HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
           "Unfolded-protein-response-associated expression", 3),
  reg_term(
    "UP", "C1_FC08",
    "GO_BP::GOBP_ENDOPLASMIC_RETICULUM_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT",
    "ER-to-Golgi vesicle transport", 4
  ),
  reg_term("DOWN", "C4_FC01", "GO_BP::GOBP_AEROBIC_RESPIRATION",
           "Aerobic respiration", 1),
  reg_term("DOWN", "C4_FC01", "GO_BP::GOBP_CELLULAR_RESPIRATION",
           "Cellular respiration", 2),
  reg_term(
    "DOWN", "C4_FC02",
    "GO_BP::GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
    "ATP synthesis-coupled electron transport", 1
  ),
  reg_term("DOWN", "C4_FC02", "GO_BP::GOBP_OXIDATIVE_PHOSPHORYLATION",
           "Oxidative phosphorylation", 2),
  reg_term(
    "DOWN", "C4_FC02",
    "REACTOME::REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "Respiratory electron transport", 3
  )
))

first_token <- function(x) {
  vapply(strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE), function(z) {
    if (!length(z) || !nzchar(z[[1]])) NA_character_ else z[[1]]
  }, character(1))
}


prepare_v1_display <- function(ora, genes, membership, upstream_qa) {
  required_ora <- c(
    "Term_Key", "Collection", "Term_Name", "Gene_Set_Size", "Overlap_N",
    "Design_ID", "Direction", "Fold_Enrichment", "FDR",
    "Full_Significant", "Display_Eligible"
  )
  required_gene <- c(
    "Design_ID", "Gene", "Direction", "N_Evaluable_Cancers", "N_Sig",
    "Recurrence", "Median_Estimate", "Median_Abs_Estimate",
    "Definition_Overlap", "ORA_Eligible", "Evidence_Rank"
  )
  required_membership <- c(
    "Design_ID", "Direction", "Gene", "Term_Key", "Collection",
    "FDR", "Full_Significant", "Display_Eligible"
  )
  if (!all(required_ora %in% names(ora))) stop("ORA输入列不完整")
  if (!all(required_gene %in% names(genes))) stop("复现基因输入列不完整")
  if (!all(required_membership %in% names(membership))) {
    stop("term-gene membership输入列不完整")
  }
  if ("Pass" %in% names(upstream_qa) && any(!upstream_qa$Pass)) {
    stop("上游复现/ORA机器QA存在未通过项")
  }

  main_terms <- copy(ora[
    Design_ID == "C1_vs_C4" &
      Full_Significant &
      Display_Eligible &
      Gene_Set_Size <= 300L &
      Overlap_N >= 5L &
      Collection != "C7"
  ])
  setorder(main_terms, Direction, FDR, -Overlap_N, -Fold_Enrichment, Term_Key)
  if (main_terms[Direction == "UP", .N] != 990L) {
    stop("C1-high非C7展示候选数量偏离冻结输入: ", main_terms[Direction == "UP", .N])
  }
  if (main_terms[Direction == "DOWN", .N] != 8L) {
    stop("C4-high非C7展示候选数量偏离冻结输入: ", main_terms[Direction == "DOWN", .N])
  }
  selected_terms <- merge(
    term_registry,
    cluster_registry[
      ,
      .(
        Direction, State, Cluster_ID, Cluster_Order, Cluster_Label,
        Parent_Domain, Cluster_Color,
        Cluster_Semantic_Boundary = Semantic_Boundary
      )
    ],
    by = c("Direction", "Cluster_ID"),
    all.x = TRUE,
    sort = FALSE
  )
  selected_terms <- merge(
    selected_terms,
    main_terms,
    by = c("Direction", "Term_Key"),
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_Manual", "")
  )
  if (nrow(selected_terms) != nrow(term_registry) ||
      any(is.na(selected_terms$Collection))) {
    stop("人工代表term未全部匹配正式候选")
  }
  if (any(selected_terms$Collection == "C7") ||
      any(!selected_terms$Full_Significant) ||
      any(!selected_terms$Display_Eligible)) {
    stop("人工代表term存在非正式主展示记录")
  }
  setorder(selected_terms, Cluster_Order, Term_Order)
  selected_terms[, Term_Node_ID := paste("TERM", Cluster_ID, Term_Key, sep = "::")]

  selected_membership <- unique(membership[
    Design_ID == "C1_vs_C4" &
      Term_Key %chin% selected_terms$Term_Key,
    .(Direction, Gene, Term_Key, Collection)
  ])
  selected_membership <- merge(
    selected_membership,
    selected_terms[
      ,
      .(
        Direction, Term_Key, Cluster_ID, Cluster_Order, Cluster_Label,
        Term_Order, Term_Node_ID, Display_Term_Label
      )
    ],
    by = c("Direction", "Term_Key"),
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(selected_membership$Cluster_ID))) {
    stop("代表term membership映射失败")
  }

  gene_evidence <- copy(genes[Design_ID == "C1_vs_C4" & ORA_Eligible])
  gene_support <- selected_membership[
    ,
    .(
      N_Selected_Terms = uniqueN(Term_Key),
      N_Selected_Clusters = uniqueN(Cluster_ID),
      Selected_Term_Keys = paste(sort(unique(Term_Key)), collapse = ";"),
      Selected_Cluster_IDs = paste(sort(unique(Cluster_ID)), collapse = ";")
    ),
    by = .(Direction, Gene)
  ]
  gene_candidates <- merge(
    gene_support,
    gene_evidence,
    by = c("Direction", "Gene"),
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(gene_candidates$N_Sig))) {
    stop("代表term成员基因未全部匹配复现证据")
  }
  setorderv(
    gene_candidates,
    c(
      "N_Sig", "Recurrence", "Median_Abs_Estimate",
      "N_Selected_Terms", "Gene"
    ),
    c(-1L, -1L, -1L, -1L, 1L),
    na.last = TRUE
  )
  gene_candidates[, Display_Evidence_Rank := seq_len(.N)]

  # 每个term至少3个、每个cluster至少15个C1基因；C4所有可连接基因保留。
  term_seeds <- selected_membership[
    gene_candidates[
      ,
      .(Direction, Gene, Display_Evidence_Rank)
    ],
    on = .(Direction, Gene),
    nomatch = 0L
  ][
    order(Cluster_Order, Term_Order, Display_Evidence_Rank, Gene),
    head(.SD, min(3L, .N)),
    by = .(Direction, Term_Key)
  ][, .(Direction, Gene)]

  cluster_seeds <- selected_membership[
    gene_candidates[
      ,
      .(Direction, Gene, Display_Evidence_Rank)
    ],
    on = .(Direction, Gene),
    nomatch = 0L
  ][
    Direction == "UP",
    .SD[order(Display_Evidence_Rank, Gene)][
      !duplicated(Gene)
    ][seq_len(min(15L, .N))],
    by = Cluster_ID
  ][, .(Direction, Gene)]

  c4_connected <- gene_candidates[Direction == "DOWN", .(Direction, Gene)]
  seed_pairs <- unique(rbindlist(list(term_seeds, cluster_seeds, c4_connected)))
  if (nrow(seed_pairs) > 240L) {
    stop("最小覆盖种子超过240个gene node")
  }

  target_gene_n <- 220L
  seed_keys <- paste(seed_pairs$Direction, seed_pairs$Gene, sep = "\r")
  remaining <- gene_candidates[
    !(paste(Direction, Gene, sep = "\r") %chin% seed_keys)
  ][order(Display_Evidence_Rank)]
  fill_n <- max(0L, target_gene_n - nrow(seed_pairs))
  selected_pairs <- unique(rbindlist(list(
    seed_pairs,
    remaining[seq_len(min(fill_n, .N)), .(Direction, Gene)]
  )))
  if (nrow(selected_pairs) < 200L || nrow(selected_pairs) > 240L) {
    stop("最终展示gene node数量不在200–240: ", nrow(selected_pairs))
  }

  display_genes <- merge(
    selected_pairs,
    gene_candidates,
    by = c("Direction", "Gene"),
    all.x = TRUE,
    sort = FALSE
  )
  setorderv(
    display_genes,
    c(
      "Direction", "N_Sig", "Recurrence", "Median_Abs_Estimate",
      "N_Selected_Terms", "Gene"
    ),
    c(1L, -1L, -1L, -1L, -1L, 1L)
  )
  display_genes[, State := fifelse(Direction == "UP", "C1", "C4")]
  display_genes[, Gene_Node_ID := paste("GENE", Gene, sep = "::")]
  display_genes[, Selected_For_Cytoscape := TRUE]
  display_genes[, Display_Gene_Rank := seq_len(.N)]

  final_membership <- selected_membership[
    display_genes[
      ,
      .(Direction, Gene, Gene_Node_ID)
    ],
    on = .(Direction, Gene),
    nomatch = 0L
  ]
  final_membership <- unique(final_membership)
  final_membership[, Selected_For_Cytoscape := TRUE]

  directional_gene_freeze <- copy(gene_evidence)
  directional_gene_freeze[, State := fifelse(Direction == "UP", "C1", "C4")]
  directional_gene_freeze[, Connected_To_Representative_Term :=
    paste(Direction, Gene, sep = "\r") %chin%
      paste(gene_candidates$Direction, gene_candidates$Gene, sep = "\r")]
  directional_gene_freeze[, Selected_For_Cytoscape :=
    paste(Direction, Gene, sep = "\r") %chin%
      paste(display_genes$Direction, display_genes$Gene, sep = "\r")]
  setorder(directional_gene_freeze, Direction, Evidence_Rank, Gene)


  list(display = display_genes, directional = directional_gene_freeze,
       terms = selected_terms, membership = selected_membership,
       final_membership = final_membership)
}
# 忠实保留source39的C4合簇、Other和最终字段/排序变换。
refine_v2_display <- function(display_v1, directional) {
  c4_cluster_id <- "C4_FC01"
  other_cluster_id <- "OTHER"
  display <- copy(display_v1)
  display[Direction == "DOWN", Selected_Cluster_IDs := c4_cluster_id]
  display[, N_Selected_Clusters := uniqueN(na.omit(unlist(strsplit(
    paste(Selected_Cluster_IDs, collapse = ";"), ";", fixed = TRUE
  )))), by = Gene]
  display[, `:=`(
    Display_Assignment = "FUNCTIONAL_TERM",
    Functional_Term_Mapped = TRUE,
    Other_Reason = NA_character_
  )]

  other_genes <- directional[
    State == "C4" & Recurrent & ORA_Eligible & !Connected_To_Representative_Term
  ]
  if (nrow(other_genes) != 10L) {
    stop("预期10个C4未连接代表term基因，实际为: ", nrow(other_genes))
  }

  new_display <- other_genes[, .(
    Direction,
    Gene,
    N_Selected_Terms = 0L,
    N_Selected_Clusters = 0L,
    Selected_Term_Keys = NA_character_,
    Selected_Cluster_IDs = NA_character_,
    Design_ID,
    Feature_Label,
    Feature_Class,
    Definition_Overlap,
    Design_Type,
    Target,
    Reference,
    Evidence_Role,
    N_Design_Evaluable_Cancers,
    N_Evaluable_Cancers,
    Min_Recurrence_Cancers,
    N_Sig,
    Recurrence,
    Sig_Cancers,
    Median_Estimate,
    Median_Abs_Estimate,
    Recurrent,
    Common_Measurable,
    Bidirectional_Recurrent,
    ORA_Eligible,
    Evidence_Rank,
    Display_Evidence_Rank = Evidence_Rank,
    State,
    Gene_Node_ID = paste0("GENE::", Gene),
    Selected_For_Cytoscape = TRUE,
    Display_Gene_Rank = NA_integer_,
    Display_Assignment = "OTHER",
    Functional_Term_Mapped = FALSE,
    Other_Reason = "FORMAL_RECURRENT_GENE_NO_DISPLAYED_REPRESENTATIVE_TERM"
  )]

  display <- rbindlist(list(display, new_display), use.names = TRUE, fill = TRUE)
  setorder(display, State, -N_Sig, -Recurrence, -Median_Abs_Estimate, Gene)
  display[, Display_Gene_Rank := seq_len(.N), by = State]
  display[, Primary_Cluster_ID := fifelse(
    Display_Assignment == "OTHER",
    other_cluster_id,
    first_token(Selected_Cluster_IDs)
  )]


  display
}
main <- function() {
  ora <- fread(paths["ora_significant"])
  genes <- fread(paths["recurrent_directional_genes"])
  membership <- fread(paths["significant_term_gene_membership"])
  upstream_qa <- fread(paths["upstream_qa"])
  if (!("Pass" %in% names(upstream_qa)) || !nrow(upstream_qa) ||
      anyNA(upstream_qa$Pass) || !all(upstream_qa$Pass)) stop("Upstream QA must be nonempty and all passed.")
  v1 <- prepare_v1_display(ora, genes, membership, upstream_qa)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  final_dir <- file.path(out_dir,"01_mrna_display_preparation")
  stage_dir <- paste0(final_dir,".tmp_",Sys.getpid())
  if(file.exists(final_dir) || file.exists(stage_dir)) stop("Output stage already exists.")
  if(!dir.create(stage_dir,showWarnings=FALSE)) stop("Cannot create temporary stage.")
  # 保留37→39历史TSV序列化边界，仅写两个必要v0.1中间表。
  intermediate <- file.path(stage_dir,"intermediate_v01")
  dir.create(intermediate,showWarnings=FALSE)
  v1_display_path <- file.path(intermediate,"mrna_c1c4_extreme_display_genes_v0.1_20260730.tsv")
  v1_directional_path <- file.path(intermediate,"mrna_c1c4_extreme_directional_gene_freeze_v0.1_20260730.tsv")
  fwrite(v1$display,v1_display_path,sep="\t",quote=TRUE,na="NA")
  fwrite(v1$directional,v1_directional_path,sep="\t",quote=TRUE,na="NA")
  display_v1 <- fread(v1_display_path)
  directional <- fread(v1_directional_path)
  display <- refine_v2_display(display_v1, directional)
  terms <- v1$terms
  # 仅保留代表项输入相关QA，不声称运行网络/社区/legacy QA。
  qa <- rbindlist(list(
    data.table(Check_ID="UPSTREAM_QA_PASS", Pass=all(upstream_qa$Pass)),
    data.table(Check_ID="REPRESENTATIVE_TERMS_33", Pass=nrow(terms)==33L && uniqueN(terms$Term_Key)==33L),
    data.table(Check_ID="TERMS_FORMAL_NON_C7", Pass=all(terms$Full_Significant) && all(terms$Display_Eligible) && all(terms$Collection!="C7")),
    data.table(Check_ID="V01_DISPLAY_200_240", Pass=nrow(v1$display) %between% c(200L,240L)),
    data.table(Check_ID="DISPLAY_FUNCTIONAL_MEMBERSHIP", Pass=all(v1$display$Gene_Node_ID %chin% v1$final_membership$Gene_Node_ID)),
    data.table(Check_ID="TERMS_CONNECTED", Pass=all(terms$Term_Node_ID %chin% v1$final_membership$Term_Node_ID)),
    data.table(Check_ID="TERM_OVERLAP_MEMBERSHIP_MATCH",
      Pass=v1$membership[,.(Membership_N=uniqueN(Gene)),by=.(Direction,Term_Key)][
        terms[,.(Direction,Term_Key,Overlap_N)],on=.(Direction,Term_Key)][,all(Membership_N==Overlap_N)]),
    data.table(Check_ID="FORMAL_C1_3416", Pass=directional[State=="C1" & Recurrent & ORA_Eligible,uniqueN(Gene)]==3416L),
    data.table(Check_ID="FORMAL_C4_16", Pass=directional[State=="C4" & Recurrent & ORA_Eligible,uniqueN(Gene)]==16L),
    data.table(Check_ID="DISPLAY_230", Pass=uniqueN(display$Gene)==230L),
    data.table(Check_ID="DISPLAY_C1_214", Pass=display[State=="C1",uniqueN(Gene)]==214L),
    data.table(Check_ID="DISPLAY_C4_16", Pass=display[State=="C4",uniqueN(Gene)]==16L),
    data.table(Check_ID="C4_FUNCTIONAL6_OTHER10",
      Pass=display[State=="C4" & Display_Assignment=="FUNCTIONAL_TERM",uniqueN(Gene)]==6L &&
        display[State=="C4" & Display_Assignment=="OTHER",uniqueN(Gene)]==10L),
    data.table(Check_ID="OTHER_ONLY_UNMAPPED_C4", Pass=all(display[Display_Assignment=="OTHER",State=="C4" & !Functional_Term_Mapped])),
    data.table(Check_ID="ALL_C4_FORMAL_DISPLAYED", Pass=setequal(directional[State=="C4" & Recurrent & ORA_Eligible,Gene],display[State=="C4",Gene])),
    data.table(Check_ID="DEFINITION_MARKER_PRESERVED", Pass="Definition_Overlap" %in% names(display) && !anyNA(display$Definition_Overlap))
  ))
  if(anyNA(qa$Pass) || !all(qa$Pass)) stop("Representative-input QA failed: ",paste(qa[!Pass,Check_ID],collapse="; "))
  fwrite(display,file.path(stage_dir,"mrna_c1c4_extreme_display_genes_v0.2_20260730.tsv"),sep="	",quote=TRUE,na="NA")
  fwrite(qa,file.path(stage_dir,"mrna_representative_input_machine_qa.tsv"),sep="	",quote=TRUE,na="NA")
  fwrite(data.table(Input=names(input_paths),Path=unname(input_paths)),
         file.path(stage_dir,"mrna_representative_input_sources.tsv"),sep="	",quote=TRUE,na="NA")
  writeLines(c(
    "Source: 37_freeze_stage3_mrna_c1c4_extreme_display_v0.1_20260730.R; 39_refine_stage3_mrna_c1c4_other_legacy_cytoscape_interface_v0.2_20260730.R",
    "Fixed 33 reviewed terms; source gene selection; source C4 merge/Other. Upstream P/FDR and 22-cancer recurrence are inherited.",
    "Output is the full-schema display input for public 13. Jaccard/Louvain, Cytoscape/legacy tables and figures are not generated.",
    "Two necessary v0.1 intermediate TSVs preserve the source37 fwrite/source39 fread boundary.",
    "These run-specific QA checks concern representative input only; no complete historical network QA is claimed."
  ),file.path(stage_dir,"README.txt"))
  writeLines(capture.output(sessionInfo()),file.path(stage_dir,"sessionInfo.txt"))
  if(!file.rename(stage_dir,final_dir)) stop("Cannot commit output stage.")
  cat("COMPLETED: ",final_dir,"
",sep="")
}
tryCatch(main(),error=function(e) { message(conditionMessage(e)); quit(status=1L) })
