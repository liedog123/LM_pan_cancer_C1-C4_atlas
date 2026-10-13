# Acid scoring, AG grouping and processed-input export

Entries 01–03 continue from the frozen acid-signature gene definition, existing unlogged TPM matrices and clinical information. They do not discover or reselect the signature. The scoring source is `1.0gene_select/LM_PH_correlation.R`, whose relevant section computes scores from `final_selected_genes`; Entry 01 excludes the source's pH section; entry 04 below packages its correlation calculation for the fixed formal panel without plotting.

| Order | Entry | Inputs | Main outputs | Source |
|---|---|---|---|---|
| 01 | `01_score_frozen_acid_signature.R` | Frozen genes + clinical TSV + project TPM tree | `clinical_data_list.rdata`, `per_acid_scores_list.rdata` | `1.0gene_select/LM_PH_correlation.R` |
| 02 | `02_assign_acid_groups.R` | 01 scores/clinical list + frozen genes + same TPM tree | `clinical_data_cluster_list.rdata`, `Heatmap/*.png` | `2.0cluster_group/cluster_group.R` |
| 03 | `03_export_processed_input.R` | 01 scores + 02 clinical groups + TPM tree | `pancancer_full_data_list.rdata` | `2.0cluster_group/pancaner_data_merge.R` |
| 04 | `04_correlate_formal_ph_panel.R` | 01 scores + same 30-project TPM tree | `ph_acid_correlation_matrix.rdata`, `ph_acid_correlation_matrix.tsv` | `1.0gene_select/LM_PH_correlation.R`, formal 12 members fixed explicitly |

New outputs are independent copies. They do not replace study inputs, frozen labels or results. The processed object has the structure expected by [shared IG/C1–C4 assignment](../shared_definitions/README.md), but new output contents and numerical equivalence have not been validated.

## Configuration and execution

Use [acid_scores.paths.example.dcf](../../config/acid_scores.paths.example.dcf), [acid_groups.paths.example.dcf](../../config/acid_groups.paths.example.dcf) and [processed_input.paths.example.dcf](../../config/processed_input.paths.example.dcf). Each DCF contains exactly the fields shown in its example and one record. Replace all placeholders. Paths resolve relative to the DCF directory, including existing ancestors and symlinks; each output must be independent and absent or empty.

```sh
Rscript workflows/acid_ag_preprocessing/01_score_frozen_acid_signature.R --config config/acid_scores.local.dcf --check-config
Rscript workflows/acid_ag_preprocessing/02_assign_acid_groups.R --config config/acid_groups.local.dcf --check-config
Rscript workflows/acid_ag_preprocessing/03_export_processed_input.R --config config/processed_input.local.dcf --check-config
# Ordinary execution computes scores, clusters/plots AG, or exports the processed object:
Rscript workflows/acid_ag_preprocessing/01_score_frozen_acid_signature.R --config config/acid_scores.local.dcf
Rscript workflows/acid_ag_preprocessing/02_assign_acid_groups.R --config config/acid_groups.local.dcf
Rscript workflows/acid_ag_preprocessing/03_export_processed_input.R --config config/processed_input.local.dcf
```

Each entry supports standalone `--help`. Only `--check-config` returns before dependency loading, serialized-input loading, clinical TSV reads and output/plot creation. It checks file/path metadata only; it does not verify schema, score ordering or values. No `--test` or cancer filter is provided.

## Input schemas and source rules

`final_genes` is an Rdata file containing `final_selected_genes` with a `Gene` column. Supply the study's frozen nine-gene definition separately. The source uses those entries directly; no new gene-selection algorithm or different signature is introduced.

`tpm_root` has project directories `<TCGA-PROJECT>/`, each containing `<TCGA-PROJECT>_mrna_expr_tpm.rdata` with object `mrna_expr_tpm`: genes in rows (gene identifiers in row names) and samples in columns, on unlogged TPM scale. 01 and 02 retain the same fixed 30-project source list, without a new cancer filter. Their metadata checks require the corresponding 30 files. 03 checks available matching TPM file metadata without loading eligible-project names; required project coverage is checked during ordinary execution by the source loader.

For 01, `clinical_table` is a tab-delimited file read with `data.table::fread`, with `sample` and any downstream clinical fields. Duplicate clinical sample rows retain the first occurrence. TPM column names are truncated to 15 characters and periods become hyphens. Tumor columns use TCGA sample-type positions 14–15 < 10. Clinical samples are intersected with tumor expression columns in clinical-table order. Scoring applies `log2(TPM + 1)`, `ssgseaParam(normalize = TRUE)` and `gsva`. Output objects are named `clinical_data_list` and `acid_scores_list`, keyed by the original `TCGA-` project names.

For 02, `acid_scores` supplies `acid_scores_list`, `clinical_data` supplies `clinical_data_list`, and `final_genes` supplies the frozen definition. It intersects clinical/expression samples, selects signature genes, computes `log2(TPM + 1)` and row-wise z-scores, then Euclidean distance, `ward.D` and `cutree(k = 3)`. Cluster labels are ranked by mean acid score to assign High/Medium/Low. These map downstream to AG1/AG2/AG3.

The original implementation coerces the score matrix's first row to an unnamed numeric vector and splits it positionally by the clustering labels. It does not reorder scores by sample names. Matching upstream 01 outputs and the same clinical/TPM inputs are therefore required to preserve score and clustering sample order. No name-based reordering or numerical guard was introduced. The original clustering heatmap is retained, with output redirected to the independent `Heatmap/` directory.

For 03, `clinical_groups` contains `clinical_data_cluster_list` with `sample` and `Group`; `acid_scores` contains `acid_scores_list`. Eligibility retains the original requirement that High and Low groups both contain > 30 rows before final three-input intersection. TPM tumor filtering is preserved. Common samples are sorted, clinical rows aligned with `match`, and existing source-listed clinical columns retained. Cancer keys remove `TCGA-`; each output element contains `clinical_data`, unlogged `expression_data` and `acid_scores`.

## Explicit interface adaptations

Each Rdata is loaded in an isolated environment and only its expected named object is retrieved, so extra serialized objects cannot overwrite configured paths. Private clinical/TPM paths, load calls, output saves and the heatmap directory are replaced with the listed interface variables.

The scoring entry retains source line 9 and lines 15–91, excluding unused pH-panel definitions and downstream pH correlations/plots. AG and merge retain their analysis sections. Line-ending spaces are removed as a formatting adaptation; internal whitespace and statistical expressions are retained.

The upstream `fread` path produces data.table clinical entries. The merge source uses data.frame-style variable-column selection, which fails on a data.table in the selected test environment. For public chaining, 03 applies `lapply` in original list order: only items inheriting `data.table` are converted with `as.data.frame`; existing data.frames are returned unchanged. The original merge column-selection statement, filtering, sorting, values and Group labels remain unchanged. A tiny synthetic DT/DF check verified the adapter's names, row order, column values and Group preservation and successful source column selection. This interface issue does not establish that frozen study results are erroneous.

Source: `1.0gene_select/LM_PH_correlation.R` (source line 9 plus lines 15–91). Current source SHA-256: `29daa75a2f5b3d59f787a6d644de073d2a8e2609aa03965a8daf90b23112c641`; public `01_score_frozen_acid_signature.R` SHA-256: `ed59c4e3934e115c26d5778ac0a7edfb62a56b1b186fb04b699521a321529e0d`. Retained blocks were compared after the listed 5 path/load substitutions and line-ending-space removal; input-object retrieval and 03's type adapter are in the execution header.

Source: `2.0cluster_group/cluster_group.R` (source lines 9–133). Current source SHA-256: `821119d7fdeb33d26972a3c0d731c8c909922fd98c6024e5a6e7032cc1b4c194`; public `02_assign_acid_groups.R` SHA-256: `13222f37353701de85d2fcbc0485c85d63afbe0e8dacd63d24676e5d2782b784`. Retained blocks were compared after the listed 4 path/load substitutions and line-ending-space removal; input-object retrieval and 03's type adapter are in the execution header.

Source: `2.0cluster_group/pancaner_data_merge.R` (source lines 10–99). Current source SHA-256: `1af20e3e6dd8761981496b4fea88f97671b2ad419d845de4ee2dc2bd41044ffd`; public `03_export_processed_input.R` SHA-256: `9fa39bc3f04aac8cbb84c42db2bbba0a7fc2f0832a8c44d5b75128342fefadf9`. Retained blocks were compared after the listed 3 path/load substitutions and line-ending-space removal; input-object retrieval and 03's type adapter are in the execution header.

The comparison checks the retained source blocks after precisely listed path/load/type/line-ending adaptations; it is not a claim that an entire original script is byte-identical. Current script fingerprints are not historical-run hashes.

## Access and checks

No TPM, clinical/patient information, frozen nine-gene object, third-party source or frozen assets are distributed. Access/release and licensing remain subject to author decisions; no license is granted here. This workflow starts from existing processed TPM files, not FASTQ or raw TCGA preprocessing.

Checks used R 4.5.1 parsing, synthetic metadata-only configurations and the small synthetic type adapter example, plus source-block comparisons, tracking and public-text scans. No research object/clinical rows were loaded and no scoring, clustering, heatmap or real export was run. Real-input schemas, positional ordering, numerical reproduction, output contents, dependency completeness and end-to-end reproduction remain unverified. Initial historical software versions are unknown. See [environment](../../environment/README.md) and [input provenance](../../docs/inputs_and_provenance.md).

## Fig. 1B formal pH-panel correlation (04)

04 is a continuation from existing nine-gene `acid_scores_list`, not signature selection or rescoring. It uses [fig1b_ph_correlation.paths.example.dcf](../../config/fig1b_ph_correlation.paths.example.dcf), with exactly `acid_scores`, `tpm_root` and `output_dir`. The same standalone CLI and independent empty-output safeguards apply:

```sh
Rscript workflows/acid_ag_preprocessing/04_correlate_formal_ph_panel.R --help
Rscript workflows/acid_ag_preprocessing/04_correlate_formal_ph_panel.R --config config/fig1b_ph_correlation.local.dcf --check-config
Rscript workflows/acid_ag_preprocessing/04_correlate_formal_ph_panel.R --config config/fig1b_ph_correlation.local.dcf
```

Supply the matched 01 score object (or its existing study equivalent), and the same unlogged gene-row/sample-column TPM tree. The fixed 30-project list is unchanged. No clinical table, gene-definition object or scoring dependency is needed by 04; it uses base R and `stats::cor.test`.

The analysis fixes PDK1, ACSS2, DHCR7, ELOVL5, FDFT1, HMGCR, HMGCS1, IDI1, LDLR, MSMO1, MTHFR and SQLE in the source's initial order. The original `LM_PH_correlation.R` then overwrites this definition with seven different genes. That overwrite is deliberately excluded under the author's explicit decision to package the formal 12-member panel. This is not an assertion that the whole current source generates the frozen 12-gene result, or evidence of the historical generating run. The original file remains unchanged.

Current Fig. 1 README and V1.6 contract identify `02_plot_fig1b_ph_associated_correlation.R` as the formal consumer; `03_make_fig1B_lowph_panel1_only.R` is legacy audit. The formal consumer requires a 12×30 matrix assembled from `Fig1B_correlation_display_long.tsv`. Only its `Gene` and `Gene_Order` fields were read to confirm the membership: the same 12 genes, displayed as PDK1, IDI1, SQLE, LDLR, DHCR7, MSMO1, HMGCR, HMGCS1, FDFT1, MTHFR, ACSS2 and ELOVL5. No frozen correlation values were used in checking 04, and the display ordering is not applied to the analysis matrix.

04 retains the source correlation block: truncate TPM names to 15 characters and replace periods with hyphens; intersect expression columns with existing score columns in expression-column order; use unlogged TPM and the original score subsetting; select available pH rows; call `cor.test(..., method = "spearman")` with all other defaults. It retains rho only, removes `TCGA-` from matrix columns, and introduces no P-value output, BH, extra tumour selection, new filtering, score reordering or heatmap. Missing gene rows retain the source's initialized NA; source warnings/errors for ties, insufficient pairs, constant values or incompatible score shapes are not silently repaired. Suitable paired upstream inputs remain necessary.

The source loader is parameterized to `tpm_root` and isolated named-object retrieval, matching 01. Added save/TSV exports contain the same `correlation_matrix`, in initial member order and fixed project order. These matrices are analysis output, not the frozen display long table or a plot. The 04 implementation and its checks do not replace frozen source data.

Checks for 04 use R 4.5.1 parsing/help, synthetic metadata/path configurations and a tiny synthetic correlation example; no real TPM/score object, research correlation or figure was run. Real-input schema/order, numerical agreement with frozen results, historical run provenance and complete reproducibility remain unverified. Upstream inputs and their access/release permissions remain subject to the author. Unlike upstream study objects, the explicitly approved formal pH members are embedded in 04.


## Fig. 1C core calculations (05 and 06)

These are two distinct routes, not a 05 -> 06 chain. Both support standalone
`--help`, strict single-record DCF `--config PATH`, and metadata-only
`--check-config` with the independent absent/empty output-directory safeguards
above. Check mode reads no serialized object or expression table and loads no
analysis dependency. No cancer filter, test subset or plotting is exposed.

05 uses [integrated_acid_scores.paths.example.dcf](../../config/integrated_acid_scores.paths.example.dcf):
`clinical_groups` contains existing `clinical_data_cluster_list`; `final_genes`
contains `final_selected_genes$Gene`; `integrated_expression` is the gzipped
EB++AdjustPANCAN gene-row/sample-column tab-delimited matrix, whose first
`sample` column identifies genes; `output_dir` is independent. Source:
`2.0cluster_group/pan_cancer_ssgsea_acid_score.R`. AG assignments are reused.

05 retains first duplicated gene rows, intersects clinical sample identifiers
with expression columns, keeps clinical order and first duplicated clinical
sample, and scores the matrix directly on its original scale. No log or TPM
transformation is introduced. PKM -> PKM2 is used only when PKM is absent and
PKM2 exists. `ssgseaParam(normalize = TRUE)` and `gsva(..., verbose = FALSE)`
are retained. The source's test-only first-five-cancer branch is removed and
metadata mode is fixed to `full`; no other selection rule changes.

Outputs retain `pancancer_integrated_ssgsea_acid_scores.rdata` (named object
`pancancer_integrated_ssgsea_acid_scores`, with `score_matrix`, `sample_scores`,
`metadata`) and `pancancer_integrated_ssgsea_acid_scores_by_sample.csv`
(`sample`, `Acid_Score`, `Group`, `Cancer`). The original printed cancer
summary is also exported as `pancancer_integrated_ssgsea_acid_summary.csv`
(`Cancer`, `N`, `Median_Acid`). The panel-source's numerical
`group_by(Cancer)/summarise/arrange(desc(Tumor_Median))/Order` block adds
`Fig1C_top_cancer_group_medians.csv` (`Cancer`, `Tumor_Median`,
`High_Median`, `Low_Median`, `Tumor_N`, `High_N`, `Low_N`, `Order`).
This is analysis order; frozen display ordering and validation flags are not
applied. No new threshold is added.

06 uses [acid_tumor_normal.paths.example.dcf](../../config/acid_tumor_normal.paths.example.dcf):
`validation_data` must contain the matching prepared `full_validation_data`
with `Cancer_Type`, `Sample_Type` (Tumor/Normal), and `Acid_Score`;
`output_dir` is independent. It extracts only Acid_Score from Analysis II of
`../LM_pan_cancer_v1/5.pancancer_cluster/6.1clinical_validation_analysis.R`.
The per-module loop is explicitly restricted to Acid_Score; other historical
immune, stage and grade analyses are outside this entry.

The original cancer iteration, two-type eligibility, formula
`wilcox.test(sub_df[[mod]] ~ sub_df$Sample_Type)`, all Wilcoxon defaults,
Tumor-minus-Normal mean difference, and raw-P significance stars (<0.001,
<0.01, <0.05) are retained. No BH family is introduced or reduced.
`Stats_Tumor_vs_Normal.csv` contains Acid_Score rows only:
`Cancer`, `Module`, `Diff_Mean`, `P_Value`, `Significance`; it is not the
legacy all-module table. `Fig1C_bottom_tumor_normal_medians.csv` independently
exports per-Cancer/per-Sample_Type medians via the numerical
`group_by/summarise/pivot_wider/left_join` block from
`final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/04_make_fig1C_acid_axis_tvn.R`,
with `Cancer`, Tumor/Normal median columns and `Diff_Mean`, `P_Value`,
`Significance` (three statistical columns, no `Module`). `Diff_Mean` remains a mean difference and never represents a median
difference. This panel source is a legacy migrated computational source, not a
claim to reproduce the current formal plot.

The upstream boundary is
`../LM_pan_cancer_v1/5.pancancer_cluster/6.0prepare_clinical_validation_data.R`:
its original cancer whitelist, >=10 normal samples, TCGA sample codes <20,
`log2(TPM+1)`, and joint normalized ssGSEA of Acid plus the historical five
immune modules must be preserved. Neither current LM12 nor Acid-only rescoring
is an equivalent replacement, because joint normalization defines this input.
That entire preparation chain is not packaged in this batch. Matching prepared
object acquisition/release remains pending author publication; a placeholder
path does not establish input availability.

```sh
Rscript workflows/acid_ag_preprocessing/05_score_integrated_acid_matrix.R --config config/integrated_acid_scores.local.dcf --check-config
Rscript workflows/acid_ag_preprocessing/06_compare_acid_tumor_normal.R --config config/acid_tumor_normal.local.dcf --check-config
# Ordinary invocation without --check-config computes and writes independent results.
```

05 requires data.table, dplyr and GSVA; 06 requires dplyr, tidyr and base stats.
Input loads are isolated named-object retrievals. Current tests are limited to
R 4.5.1 parse/help, synthetic metadata/path cases, static retained-block
comparison and small synthetic non-scoring calculations. No research table,
object or schema is inspected, no GSVA/ssGSEA is loaded or executed, and no
figure or study analysis is run. Actual input schema/order, numerical agreement,
output contents, dependency completeness and end-to-end reproduction remain
unverified; this test environment is not a historical generation version.
