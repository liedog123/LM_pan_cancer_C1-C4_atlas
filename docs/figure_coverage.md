# Formal-result coverage and remaining core-analysis delivery

Public release: code/configuration only. References described as locally packaged or supplied below belong to internal preparation and are not distributed here. Software license and destination are recorded in [release scope](../docs/release_scope.md).

Updated: 2026-10-03

The workspace contains 90 entries in 15 workflows. It is a core-analysis package, not a mandatory figure-by-figure reproduction package. Study-specific computation, complete statistical families and usable analysis inputs are required. Pure plotting, layout, workbook styling and manual assembly are optional. Existing plotting entries remain available.

This map uses current figure contracts, workflow documentation and targeted source-code inspection. It is not a new statistical audit or a claim of numerical or visual reproduction. The first coverage pass did not inspect study data or workbooks. The Fig. 1B follow-up checked only the formal table's gene identity/order fields; it did not compare correlation values or load patient/expression inputs.

The [output-to-result map](result_output_mapping.md) gives the current detailed mapping and conditional run order for Fig. 1–7 / Supplementary Tables 1–7, including all four complete-result companions. It distinguishes analysis outputs from selected publication tables and does not claim newly reproduced values. Table 6 v1.1's ICB increment is included; prepared descriptive exports remain distinct from direct workflow outputs.

## Coverage by result module

| Formal scope | Existing public entries | Remaining core delivery |
|---|---|---|
| Fig. 1; Supplementary Figs. 1–2; Table 1 | [Signature selection](../workflows/signature_selection/README.md) 01–03, [Acid/AG](../workflows/acid_ag_preprocessing/README.md) 01–06 and [AG survival](../workflows/ag_survival/README.md) 00–03 | Fig. 1B is included in 04; Fig. 1C integrated scoring and separate prepared-input Acid_Score tumour–normal statistics are included in 05/06. The historical validation object's access/scoring background remains an upstream boundary. S1A/B final-nine-gene and 36-pair correlation/BH families are included in signature selection 02/03. Table 1 output/source selection is documented in the linked map; input access remains; no separate heatmap/KM layout entry is required. |
| Fig. 2; Supplementary Figs. 3–6; Table 2 | [Discovery](../workflows/tcga_discovery/README.md), [LM12](../workflows/lm12_construction/README.md), [shared definitions](../workflows/shared_definitions/README.md), [external bulk](../workflows/external_bulk/README.md), [clinical outcomes](../workflows/clinical_outcomes/README.md), [immune ecology](../workflows/immune_ecology/README.md) | No additional core inferential family was identified in the first pass. Table 2 output/source selection and branch order are documented in the linked map; actionable input access remains. GSEA curves, matrices, state/survival displays and sample heatmap renderers are optional; a missing renderer alone is not a missing analysis. |
| Fig. 3; Supplementary Fig. 7 A–F; Table 3 v1.1 | [Transcriptome/proteome](../workflows/transcriptome_proteome/README.md) 01–16 | Entry 16 generates the mRNA representative table consumed by 13 from four 02 recurrence/ORA products; fixed selection and C4/Other assignment are retained without network rendering. Table 3 v1.1 output/source mapping is documented; retain complete statistical backgrounds beyond the selected worksheet views. Network layout, fixed heatmaps and workbook formatting are optional. |
| Fig. 4; Supplementary Fig. 8; Table 4 | [Genomics](../workflows/genomics/README.md) 01–05 | Core mutation, TMB/MSI and both CNV branches are included. Table 4 and its complete CNV companion are mapped; callable/readiness/reference-input access remains open. TP53's descriptive/audit export remains a separate prepared-source boundary. 04/05 already include Fig. 4C/D plotting. TP53 protein-position/lollipop preparation is not automatically required for the core binary-mutation workflow; add it only if that descriptive result is explicitly included in the delivery scope. |
| Fig. 5; Supplementary Figs. 9–10; Table 5 | [Drug response](../workflows/drug_response/README.md) 01–06 | Table 5's complete analytical results and prepared descriptive exports are mapped; required input access and compound/reference definitions remain open. PRISM coverage remains an input, not a packaged output. Four-state patient plot values, class-display summaries and scatter renderers are optional display extensions unless they supply otherwise-unrepresented result calculations. Do not use a C1/C4 manifest as a substitute for all-four-state data if that extension is chosen. |
| Fig. 6; Supplementary Figs. 11 A–F and 12 A–C; Table 6 v1.1 | [Single cell](../workflows/single_cell/README.md) 01–18 and [ICB](../workflows/icb/README.md) 01–08 | Table 6 v1.1 and both complete companions are mapped. Processed-object/metadata/reference access and prepared descriptive exports remain boundaries. Current Fig. 6F/G consume single-cell 15; another historical melanoma programme model is not required. Plotting consumers and corrected display-label renderers are optional; complete upstream statistics remain required. |
| Fig. 7; Supplementary Fig. 13; Table 7 | [CD8](../workflows/cd8_perturbation/README.md) 01–04 analysis and 05–12 data-panel entries | Core analysis and all eight data-panel scripts are already included. Existing optional entries 10/11 retain explicit prepared trajectory/audited-sample input boundaries; resolve access or clearly delimit their usability, not the status of counts/DE itself. Table 7 and its complete gene-DE companion are mapped to the analysis outputs; no new artwork or workbook-layout script is required. |

Code coverage does not establish that required inputs have been distributed or that outputs reproduce the frozen results. Prepared objects, sample/clinical interfaces, reference resources and permissions remain distinct delivery requirements. Table mappings are documentation/source based, not per-sheet numerical comparisons.

## Input boundaries and scoped completion

The [input-access inventory](input_access_and_permissions.md) distinguishes documented external sources, conditional generation routes and unsupplied prepared inputs. Nineteen unchanged reference copies are locally supplied across twenty-three DCFs, including the 142 candidates, complete curated pathway candidates and family dictionary; no further data upload or prepared-input packaging is part of this code-organization closeout. Public visibility does not establish current download availability or redistribution rights. Use [Supplementary Table 2 and the published nine-gene definitions](reference_file_candidates.md#supplementary-material-first) where applicable, without duplicate uploads solely for format conversion.

Prepared patient/expression objects and other matching inputs remain dependencies for running the affected entries, with known source/format details recorded rather than implied delivery. mRNA 02 → 16 provides a representative-table generation route but not its upstream study inputs or the drug cluster registry. [Output mapping](result_output_mapping.md) covers Fig. 1–7 / Tables 1–7 and complete-result companions; [environment evidence](../environment/README.md) distinguishes historical versions from packaging checks. Scoped code organization can close with the reader entry and limited static check; ordinary-run compatibility, numerical reproduction, actual input access, CD8 data submission, provider permissions, license and publication remain separate decisions.

## Supplementary Fig. 1A/B: independent correlation families

Signature selection 02 uses the final nine genes, raw stored TPM and the first tumour-program score row on the original all-sample intersection. It retains Pearson correlation and nine-test BH within each cancer. Entry 03 independently uses all 36 gene pairs with Spearman correlation (`exact = FALSE`) and 36-test BH within each cancer; the score object defines sample intersections, not correlation values. Both retain 21-project coverage, ordering, P/FDR and significance fields. Their TSVs provide the complete 189- and 756-row statistical grids, not heatmaps. The 142-candidate screening family does not replace either family. Matching prepared inputs and numerical equivalence remain separate delivery/verification requirements.

## Fig. 3 cross-layer representative annotation

Transcriptome/proteome 16 extracts source 37's fixed 33-term registry and gene-selection rules, followed by source 39's C4 cluster merge and addition of unmapped formal C4 genes to Other. It consumes the significant ORA, directional recurrent-gene, term-gene membership and upstream-QA tables produced by 02. Supplied P/FDR and recurrence values are retained; there is no new DE, ORA or BH family.

The output mrna_c1c4_extreme_display_genes_v0.2_20260730.tsv supplies 13's mrna_display input. Entry 13 uses this for representative/cluster annotations and overlap summaries, not to determine core concordance classes. The preparation does not recreate Jaccard/Louvain audit fields, Cytoscape nodes/edges, legacy CSVs or plots. Matching upstream inputs and study numerical equivalence remain unverified.

## Fig. 1B: fixed 12-gene analysis

Acid/AG 04 continues from existing nine-gene acid scores and matching per-project TPM. It uses PDK1, ACSS2, DHCR7, ELOVL5, FDFT1, HMGCR, HMGCS1, IDI1, LDLR, MSMO1, MTHFR and SQLE. These members match the formal Fig. 1B identity fields. The public numerical matrix follows source definition order; the formal heatmap has a separate display order.

The source is the correlation section of `1.0gene_select/LM_PH_correlation.R`. Its initial 12-member definition is explicitly retained for the public entry; the subsequent seven-member reassignment is omitted. This is a documented scope adaptation, not a claim that the complete source file or a historical run is identical. The original script and frozen data are unchanged. Correlation values have not been compared with the freeze.

## Fig. 1C: two independent input routes

Acid/AG 05 retains direct EB++AdjustPANCAN integrated expression, the source gene/sample deduplication and ordering, PKM-to-PKM2 fallback and normalized ssGSEA. It reuses existing clinical AG labels; it neither transforms the integrated matrix with an additional log nor recalculates AG. Full-scope numerical objects/tables and descriptive summaries are exported independently.

Acid/AG 06 starts from the matching `full_validation_data` object, not from 05. It retains the source unpaired/default Wilcoxon formula, raw-P significance labels and Tumor-minus-Normal mean difference for Acid_Score only. Group medians are reported separately, not labelled as mean differences. Legacy immune-module, stage/grade and drawing branches are not included.

The validation object's upstream source used a cancer whitelist, sample-type codes below 20 and at least ten normal samples, with per-cancer log2(TPM+1) and joint normalized ssGSEA of Acid plus the historical five immune modules. That preparation is an explicit input boundary, not a newly packaged or executed workflow. Do not substitute current LM12, the integrated 05 scores or an Acid-only rescoring. Input access remains unresolved; packaged computation and numerical reproduction are different claims.

## Focused source pointers

Paths are relative to the original research project and identify provenance, not files promised inside the standalone release.

- Fig. 1C integrated route: `2.0cluster_group/pan_cancer_ssgsea_acid_score.R`. Legacy input route: `final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/04_make_fig1C_acid_axis_tvn.R`; current display consumer is `03_plot_fig1c_pan_cancer_acid_distribution.R` in the same directory.
- Tumour–normal route: legacy `../LM_pan_cancer_v1/5.pancancer_cluster/6.0prepare_clinical_validation_data.R` and `6.1clinical_validation_analysis.R`. A matching prepared-input continuation is permitted; do not transplant a different scoring background.
- S1A/B: `r01_f1s1_01_prepare_and_plot_final_nine_panels.R` and `r01_f1s1_02_prepare_and_plot_pairwise_synergy_and_counts.R` under `final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/`. The statistical sections are packaged in signature selection 02/03; drawing remains optional.
- Fig. 1B: `1.0gene_select/LM_PH_correlation.R` and the current consumer `final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/02_plot_fig1b_ph_associated_correlation.R`.
- Existing CD8 prepared boundaries are documented in its [source mapping](../workflows/cd8_perturbation/source_mapping.md).

## Completion sequence and boundaries

The identified Fig. 1 computation gaps and mRNA representative-input generation route are packaged. The reader entry and one limited documentation/configuration/file-inventory check close the current code-organization scope, reusing existing interface and reference checks. Missing inputs are documented execution dependencies; data acquisition/submission, permissions, licensing, actual publication and ordinary-run/numerical validation remain separate work. Do not turn optional display work or further data uploads into prerequisites for this scoped completion.

Retain full upstream statistical families and supplied P/FDR values even when a figure shows a subset. Do not rerun or re-audit earlier entries merely because this map changes. Real-data schema checks, numerical comparisons and rendering need a separately agreed scope. No release license, data-distribution permission or GitHub publication is granted by this workspace.
