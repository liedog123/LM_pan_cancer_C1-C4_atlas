# Workflow guide

Choose the result you want to examine, obtain matching inputs, then use the corresponding module. Numbered arrows describe dependencies within that module; they are not a repository-wide run-all sequence. Existing matching intermediates can replace earlier producers when the documented interface matches.

## Entry roles

| Role | Reading priority |
|---|---|
| Main analysis | Primary computations and preparations needed for the main routes below. A finalizer that calculates statistics or controls the full family stays here. |
| Supporting | Context, sensitivity, cross-platform or post-hoc branches and their required preparations. Required when reproducing the corresponding reported result. |
| Display/export | Prepared-result display or format conversion. Some entries reconstruct curves or calculate descriptive summaries; the label does not imply zero computation or an interchangeable input. |

These roles organize navigation, not scientific importance, deletion eligibility or permission to reduce test families. Individual scripts can contain primary and sensitivity outputs together; their complete implementations are retained.

## Before running

Read the module README for the exact DCF fields, input schemas and permitted starting intermediates. Obtain [required inputs](input_access_and_permissions.md) and [reference resources](../references/README.md) separately; this public release contains code/configuration only. Use a separate absent/empty output directory per entry. Configuration metadata checks do not validate research values. Check [dependencies and verification evidence](../environment/README.md) for the chosen module.

Current navigation grouping: **66 main-analysis/preparation entries, 17 supporting entries and 7 display/export entries**. All 90 scripts and their 90 example DCFs are retained.

## Signature and states

| Module | Main route | Other entries |
|---|---|---|
| [Signature selection](../workflows/signature_selection/README.md) | 01 selection; 02 and 03 are independent final-nine correlation analyses | — |
| [Acid/AG preprocessing](../workflows/acid_ag_preprocessing/README.md) | 01 → 02 → 03; 01 → 04; 05 and 06 use separate prepared inputs | — |
| [TCGA discovery](../workflows/tcga_discovery/README.md) | 01 → 02 | — |
| [LM12 construction](../workflows/lm12_construction/README.md) | 01 → 02 | — |
| [Shared definitions](../workflows/shared_definitions/README.md) | 01 → 02; 03 converts outputs from both | Display/export: 03 |

Selection can be skipped when examining downstream results with the adopted frozen definitions. Signature 02/03 use the formal nine genes independently of a new 01 selection. Acid 05 uses the supplied integrated matrix/groups; 06 uses the historical validation object and does not consume 05 output. New assignments or canonical-format exports do not replace frozen study labels.

## Bulk validation and clinical context

| Module | Main route | Other entries |
|---|---|---|
| [External bulk](../workflows/external_bulk/README.md) | 01 → 02 | — |
| [AG survival](../workflows/ag_survival/README.md) | 00 → 01 → 02; 03 uses a matching survival list separately | — |
| [Clinical outcomes](../workflows/clinical_outcomes/README.md) | 01 | — |
| [Immune ecology](../workflows/immune_ecology/README.md) | 01 → 02 | — |

External bulk evaluates Acid-linked enrichment directions, not external C1–C4 classification. AG survival retains its own High/Low clinical workflow; C1–C4 clinical models use a distinct patient-representative interface. Immune exports need both 01 core and 02 evidence outputs. AG 03 calculates pooled log-rank statistics and is not classified as a display-only step.

## Multiomics

| Module | Main route | Other entries |
|---|---|---|
| [Transcriptome/proteome](../workflows/transcriptome_proteome/README.md) | mRNA: 01 → 02; RPPA: 04 → 03; lncRNA: 04 → 05 → 06 → 07 → 08; edges: 09 → (10, 11 for E1/E3) → 12 | Supporting: 13, 14, 15, 16 |
| [Genomics](../workflows/genomics/README.md) | 01–05 are independent analysis branches | — |

For mRNA–RPPA concordance, run mRNA 01 → 02 → 16 → 13 and supply separate RPPA 04 → 03 products. Subtype context follows 14 → 15. Edge GSEA uses 10 for E2 and separate 11 runs for E1/E3, with the configured ranking scopes retained before 12 combines them. Complete statistical families and negative/non-estimable results remain in their original consumers. Genomic branches 01–05 use separate prepared inputs.

## Predicted drug response

| Module | Main route | Other entries |
|---|---|---|
| [Drug response](../workflows/drug_response/README.md) | 01; 05 → 06 | Supporting: 02, 03, 04 |

01 starts after prediction/GLDS correction. Supporting 02 supplies internal OOF products to 03, which also needs 01 discovery counts and an independently prepared PRISM coverage table; matching existing support products can be used. 04 is an independent shared-PC1 sensitivity. Main 05 → 06 supplies the complete programme–feature association families. Predictions and cell-line support do not establish clinical efficacy.

## Cellular context

| Module | Main route | Other entries |
|---|---|---|
| [Single cell](../workflows/single_cell/README.md) | States: 04 → 05 → (01, 02) → 03; pseudobulk: 06 → 07 → 08 → 09 → 10 → 11; melanoma: 14 → 15 | Supporting: 12, 13, 16, 17, 18 |

State mapping 03 does not generate the independent G1 mapping/gates needed by pseudobulk 06–11. Keep stability/evidence steps 09/10 in the pathway route. Supporting 12 and 13 are independent Bassez sensitivities; do not run 12 → 13. Communication follows 16 → 17 → 18 and retains all forty state/replicate networks. Independent melanoma 14/15 also supplies the CD8 bridge inputs.

## Immunotherapy context

| Module | Main route | Other entries |
|---|---|---|
| [ICB](../workflows/icb/README.md) | (01 → 02) + (03 → 04) → 05 | Supporting: 06, 07, 08 |

The fixed three-cohort route is (01 → 02) + (03 → 04) → 05. Supporting Liu 06 → 07 and Gide 08 are separate post-hoc routes; do not pool their tests into the fixed families. Gide 08 also uses 01 patient scores and transcriptome/proteome 14 programme memberships. 05 requires an absent output directory with an existing parent.

## CD8 time course

| Module | Main route | Other entries |
|---|---|---|
| [CD8 perturbation](../workflows/cd8_perturbation/README.md) | 01 QC and 02 DE independently use counts; 02 → 03 → 04 | Supporting: 05, 06; Display/export: 07, 08, 09, 10, 11, 12 |

01 QC and 02 DE start independently from the same counts/design; review QC before DE. 03 needs separate melanoma 14/15 products, then 04 continues the phenotype comparison. Supporting 05 uses 02 DE and 03 leading edges; 06 uses 02 DE. Display 07 uses 01 PCA coordinates, 08 uses 04 GSEA tables, and 09 reconstructs running ES from 02 ranks/03 references while reusing supplied NES/P/FDR. Display 10 needs an independently prepared trajectory table; 11 needs audited metadata plus 01 filtered logCPM and 02 DE; 12 uses 04 marker expression. No direct fully connected 02 → 10 generation route is claimed.

## Complete entry index

<details>
<summary>Expand all 90 entries and their roles</summary>

### Signature selection

Input starting point: Adopted candidates, saved DEG/programme scores, TPM; frozen nine genes for 02/03.

| Entry | Role | Purpose |
|---|---|---|
| [01_select_signature_all_sample.R](../workflows/signature_selection/01_select_signature_all_sample.R) | Main analysis | Candidate filtering, correlations and cross-cancer signature selection |
| [02_correlate_final_nine_tumor_program.R](../workflows/signature_selection/02_correlate_final_nine_tumor_program.R) | Main analysis | Formal final-nine tumour-programme correlations |
| [03_correlate_final_nine_gene_pairs.R](../workflows/signature_selection/03_correlate_final_nine_gene_pairs.R) | Main analysis | Formal final-nine gene-pair correlations |

### Acid/AG preprocessing

Input starting point: Frozen acid genes, clinical/TPM inputs; integrated matrix for 05; historical validation object for 06.

| Entry | Role | Purpose |
|---|---|---|
| [01_score_frozen_acid_signature.R](../workflows/acid_ag_preprocessing/01_score_frozen_acid_signature.R) | Main analysis | Frozen acid-signature scoring |
| [02_assign_acid_groups.R](../workflows/acid_ag_preprocessing/02_assign_acid_groups.R) | Main analysis | Within-cancer AG assignment |
| [03_export_processed_input.R](../workflows/acid_ag_preprocessing/03_export_processed_input.R) | Main analysis | Processed-object assembly for state assignment |
| [04_correlate_formal_ph_panel.R](../workflows/acid_ag_preprocessing/04_correlate_formal_ph_panel.R) | Main analysis | Formal twelve-member pH-panel correlations |
| [05_score_integrated_acid_matrix.R](../workflows/acid_ag_preprocessing/05_score_integrated_acid_matrix.R) | Main analysis | Integrated-matrix acid scoring |
| [06_compare_acid_tumor_normal.R](../workflows/acid_ag_preprocessing/06_compare_acid_tumor_normal.R) | Main analysis | Historical tumour–normal comparison |

### TCGA discovery

Input starting point: Matched High/Low processed object; provider ImmPort gene sets.

| Entry | Role | Purpose |
|---|---|---|
| [01_build_high_low_deg.R](../workflows/tcga_discovery/01_build_high_low_deg.R) | Main analysis | High/Low differential-expression rankings |
| [02_run_immport_gsea.R](../workflows/tcga_discovery/02_run_immport_gsea.R) | Main analysis | Per-cancer ImmPort GSEA |

### LM12 construction

Input starting point: Formal GSEA tables plus complete curated candidate terms/family dictionary.

| Entry | Role | Purpose |
|---|---|---|
| [01_build_votes_and_stable_genes.R](../workflows/lm12_construction/01_build_votes_and_stable_genes.R) | Main analysis | Recurrence votes and stable leading-edge membership |
| [02_finalize_lm12_gene_sets.R](../workflows/lm12_construction/02_finalize_lm12_gene_sets.R) | Main analysis | Candidate/family audits and LM12 finalization |

### Shared definitions

Input starting point: Processed pan-cancer object and frozen LM12; assignment outputs do not replace canonical study labels.

| Entry | Role | Purpose |
|---|---|---|
| [01_assign_ig_and_c1c4_from_processed_input.R](../workflows/shared_definitions/01_assign_ig_and_c1c4_from_processed_input.R) | Main analysis | IG and C1–C4 assignment |
| [02_validate_state_axes_from_labels.R](../workflows/shared_definitions/02_validate_state_axes_from_labels.R) | Main analysis | Joint-state axes and label validation |
| [03_export_canonical_labels.R](../workflows/shared_definitions/03_export_canonical_labels.R) | Display/export | Canonical-format label/score export |

### External bulk

Input starting point: Cohort expression/annotation and frozen signature/LM12 definitions.

| Entry | Role | Purpose |
|---|---|---|
| [01_build_external_scores.R](../workflows/external_bulk/01_build_external_scores.R) | Main analysis | Cohort score preparation and coverage |
| [02_run_lm12_gene_level_gsea.R](../workflows/external_bulk/02_run_lm12_gene_level_gsea.R) | Main analysis | External Acid-linked gene-level GSEA and sensitivity |

### AG survival

Input starting point: AG clinical lists, project vector and frozen display order; 03 retains pooled log-rank calculation.

| Entry | Role | Purpose |
|---|---|---|
| [00_prepare_ag_km.R](../workflows/ag_survival/00_prepare_ag_km.R) | Main analysis | Clinical-list preparation and full per-project KM analysis |
| [01_run_ag_cox.R](../workflows/ag_survival/01_run_ag_cox.R) | Main analysis | Within-cancer AG Cox models |
| [02_finalize_ag_cox.R](../workflows/ag_survival/02_finalize_ag_cox.R) | Main analysis | Complete-family Cox BH and final summary |
| [03_export_formal_ag_km.R](../workflows/ag_survival/03_export_formal_ag_km.R) | Main analysis | Formal KM tables and pooled log-rank statistics |

### Clinical outcomes

Input starting point: Representative-patient clinical inputs and canonical C1–C4 labels.

| Entry | Role | Purpose |
|---|---|---|
| [01_run_stage1_clinical_models.R](../workflows/clinical_outcomes/01_run_stage1_clinical_models.R) | Main analysis | C1/C4 clinical Cox, model evaluability and cross-cancer pooling |

### Immune ecology

Input starting point: Patient/state/deconvolution inputs; downstream exports need both output directories.

| Entry | Role | Purpose |
|---|---|---|
| [01_run_stage2_immune_ecology.R](../workflows/immune_ecology/01_run_stage2_immune_ecology.R) | Main analysis | Complete within-cancer immune/stromal comparisons |
| [02_finalize_stage2_evidence.R](../workflows/immune_ecology/02_finalize_stage2_evidence.R) | Main analysis | Evidence/convergence finalization |

### Transcriptome/proteome

Input starting point: Matched RNA/RPPA/readiness/annotation inputs and fixed definitions.

| Entry | Role | Purpose |
|---|---|---|
| [01_run_mrna_models.R](../workflows/transcriptome_proteome/01_run_mrna_models.R) | Main analysis | mRNA within-cancer models |
| [02_run_mrna_recurrence_ora.R](../workflows/transcriptome_proteome/02_run_mrna_recurrence_ora.R) | Main analysis | mRNA recurrence and ORA |
| [03_finalize_rppa_recurrence_tcpa.R](../workflows/transcriptome_proteome/03_finalize_rppa_recurrence_tcpa.R) | Main analysis | RPPA full-panel correction, recurrence and TCPA annotation |
| [04_run_lncrna_rppa_models.R](../workflows/transcriptome_proteome/04_run_lncrna_rppa_models.R) | Main analysis | lncRNA/RPPA model preparation |
| [05_freeze_lncrna_recurrence.R](../workflows/transcriptome_proteome/05_freeze_lncrna_recurrence.R) | Main analysis | lncRNA recurrence finalization |
| [06_run_lncrna_gene_coexpression.R](../workflows/transcriptome_proteome/06_run_lncrna_gene_coexpression.R) | Main analysis | lncRNA–mRNA coexpression |
| [07_build_lncrna_edge_recurrence.R](../workflows/transcriptome_proteome/07_build_lncrna_edge_recurrence.R) | Main analysis | Linked-gene recurrence |
| [08_run_lncrna_linked_mrna_ora.R](../workflows/transcriptome_proteome/08_run_lncrna_linked_mrna_ora.R) | Main analysis | Linked-mRNA ORA |
| [09_run_mrna_boundary_models.R](../workflows/transcriptome_proteome/09_run_mrna_boundary_models.R) | Main analysis | Three-edge mRNA model preparation |
| [10_run_e2_full_rank_gsea.R](../workflows/transcriptome_proteome/10_run_e2_full_rank_gsea.R) | Main analysis | E2 GSEA |
| [11_run_e1_e3_full_rank_gsea.R](../workflows/transcriptome_proteome/11_run_e1_e3_full_rank_gsea.R) | Main analysis | Separate E1/E3 GSEA runs |
| [12_combine_three_edge_gsea.R](../workflows/transcriptome_proteome/12_combine_three_edge_gsea.R) | Main analysis | Combined three-edge evidence |
| [13_run_mrna_rppa_concordance.R](../workflows/transcriptome_proteome/13_run_mrna_rppa_concordance.R) | Supporting | mRNA–RPPA concordance |
| [14_prepare_subtype_context.R](../workflows/transcriptome_proteome/14_prepare_subtype_context.R) | Supporting | Known-subtype patient/programme preparation |
| [15_analyze_subtype_context.R](../workflows/transcriptome_proteome/15_analyze_subtype_context.R) | Supporting | Subtype-adjusted and sensitivity analyses |
| [16_prepare_mrna_representatives.R](../workflows/transcriptome_proteome/16_prepare_mrna_representatives.R) | Supporting | Representative-mRNA input selection for 13 |

### Genomics

Input starting point: Callable mutation, TMB/MSI and CNV/event matrices plus readiness/state inputs.

| Entry | Role | Purpose |
|---|---|---|
| [01_run_mutation_c1c4.R](../workflows/genomics/01_run_mutation_c1c4.R) | Main analysis | Callable mutation comparisons |
| [02_run_tmb_msi_c1c4.R](../workflows/genomics/02_run_tmb_msi_c1c4.R) | Main analysis | TMB/MSI comparisons |
| [03_run_cnv_burden_c1c4.R](../workflows/genomics/03_run_cnv_burden_c1c4.R) | Main analysis | C1/C4 CNV-burden comparisons |
| [04_run_cnv_burden_four_state.R](../workflows/genomics/04_run_cnv_burden_four_state.R) | Main analysis | Four-state CNV omnibus analysis |
| [05_run_cnv_alteration_regions.R](../workflows/genomics/05_run_cnv_alteration_regions.R) | Main analysis | Complete gene-level CNV tests and region summaries |

### Drug response

Input starting point: Corrected drug predictions for 01; matched numeric bundles/reference inputs for 05/06.

| Entry | Role | Purpose |
|---|---|---|
| [01_c1c4_from_corrected_predictions.R](../workflows/drug_response/01_c1c4_from_corrected_predictions.R) | Main analysis | C1/C4 corrected predicted-drug comparisons |
| [02_run_internal_support.R](../workflows/drug_response/02_run_internal_support.R) | Supporting | Repeated OOF internal modelability |
| [03_run_prism_crossplatform.R](../workflows/drug_response/03_run_prism_crossplatform.R) | Supporting | GDSC–PRISM cross-platform support |
| [04_run_shared_axis_sensitivity.R](../workflows/drug_response/04_run_shared_axis_sensitivity.R) | Supporting | Shared-PC1 sensitivity |
| [05_build_shared_numeric_inputs.R](../workflows/drug_response/05_build_shared_numeric_inputs.R) | Main analysis | Shared numeric programme/feature preparation |
| [06_run_programme_feature_associations.R](../workflows/drug_response/06_run_programme_feature_associations.R) | Main analysis | Complete programme–feature association families |

### Single cell

Input starting point: Prepared objects/audits; pseudobulk additionally needs independently prepared G1 mapping/gates.

| Entry | Role | Purpose |
|---|---|---|
| [01_assign_study_compartment_ig_states.R](../workflows/single_cell/01_assign_study_compartment_ig_states.R) | Main analysis | Local IG assignment |
| [02_score_acid_assign_ag_states.R](../workflows/single_cell/02_score_acid_assign_ag_states.R) | Main analysis | Local AG scoring/assignment |
| [03_map_sample_c1c4_states.R](../workflows/single_cell/03_map_sample_c1c4_states.R) | Main analysis | Sample-level C1–C4 mapping |
| [04_score_initial_lm12.R](../workflows/single_cell/04_score_initial_lm12.R) | Main analysis | Initial LM12 scoring |
| [05_refine_and_rescore_lm12.R](../workflows/single_cell/05_refine_and_rescore_lm12.R) | Main analysis | Member refinement/rescoring |
| [06_build_pseudobulk_counts.R](../workflows/single_cell/06_build_pseudobulk_counts.R) | Main analysis | Pseudobulk aggregation |
| [07_build_pseudobulk_model_contract.R](../workflows/single_cell/07_build_pseudobulk_model_contract.R) | Main analysis | Model/gate/sample-set preparation |
| [08_run_primary_pseudobulk_models.R](../workflows/single_cell/08_run_primary_pseudobulk_models.R) | Main analysis | Primary pseudobulk gene models |
| [09_run_pseudobulk_stability.R](../workflows/single_cell/09_run_pseudobulk_stability.R) | Main analysis | Stability and cohort/LOSO evidence |
| [10_summarize_pseudobulk_evidence.R](../workflows/single_cell/10_summarize_pseudobulk_evidence.R) | Main analysis | Integrated gene evidence |
| [11_run_pseudobulk_pathways.R](../workflows/single_cell/11_run_pseudobulk_pathways.R) | Main analysis | Complete pathway inference |
| [12_run_bassez_exclusion_sensitivity.R](../workflows/single_cell/12_run_bassez_exclusion_sensitivity.R) | Supporting | Bassez exclusion sensitivity |
| [13_run_bassez_marker_stability.R](../workflows/single_cell/13_run_bassez_marker_stability.R) | Supporting | Independent Bassez marker stability |
| [14_run_melanoma_pseudobulk_programs.R](../workflows/single_cell/14_run_melanoma_pseudobulk_programs.R) | Main analysis | Melanoma pseudobulk/programme preparation |
| [15_run_melanoma_patient_deg_gsea.R](../workflows/single_cell/15_run_melanoma_patient_deg_gsea.R) | Main analysis | Melanoma patient-level DEG/GSEA |
| [16_build_bassez_balanced_cellchat_input.R](../workflows/single_cell/16_build_bassez_balanced_cellchat_input.R) | Supporting | Balanced CellChat input preparation |
| [17_run_bassez_cellchat_rep1.R](../workflows/single_cell/17_run_bassez_cellchat_rep1.R) | Supporting | CellChat first replicate |
| [18_run_bassez_cellchat_pilot10.R](../workflows/single_cell/18_run_bassez_cellchat_pilot10.R) | Supporting | Complete ten-replicate communication consensus |

### ICB

Input starting point: Harmonized expression/outcomes, prior-treatment mapping, Braun workbook/crosswalk and fixed memberships.

| Entry | Role | Purpose |
|---|---|---|
| [01_prepare_riaz_gide.R](../workflows/icb/01_prepare_riaz_gide.R) | Main analysis | Riaz/Gide assignment and scores |
| [02_analyze_riaz_gide.R](../workflows/icb/02_analyze_riaz_gide.R) | Main analysis | Riaz/Gide response analyses |
| [03_prepare_braun_official.R](../workflows/icb/03_prepare_braun_official.R) | Main analysis | Official Braun assignment/preparation |
| [04_analyze_braun_official.R](../workflows/icb/04_analyze_braun_official.R) | Main analysis | Braun response analyses |
| [05_integrate_three_cohort.R](../workflows/icb/05_integrate_three_cohort.R) | Main analysis | Fixed three-cohort integration and BH |
| [06_prepare_liu_exploratory.R](../workflows/icb/06_prepare_liu_exploratory.R) | Supporting | Separate Liu assignment/preparation |
| [07_analyze_liu_exploratory.R](../workflows/icb/07_analyze_liu_exploratory.R) | Supporting | Separate Liu primary/broad response analyses |
| [08_run_gide_c1_sensitivity.R](../workflows/icb/08_run_gide_c1_sensitivity.R) | Supporting | Post-hoc Gide C1 sensitivity |

### CD8 perturbation

Input starting point: Matched counts/design; 03 also needs melanoma 14/15 products and fixed reference inputs.

| Entry | Role | Purpose |
|---|---|---|
| [01_qc_and_pca.R](../workflows/cd8_perturbation/01_qc_and_pca.R) | Main analysis | Counts QC, filtering/TMM and PCA |
| [02_differential_response.R](../workflows/cd8_perturbation/02_differential_response.R) | Main analysis | Donor-adjusted DE and interaction models |
| [03_fixed_programme_gsea.R](../workflows/cd8_perturbation/03_fixed_programme_gsea.R) | Main analysis | Fixed-programme GSEA and melanoma bridge |
| [04_hcl_control_phenotype.R](../workflows/cd8_perturbation/04_hcl_control_phenotype.R) | Main analysis | HCl-control/phenotype continuation |
| [05_plot_la_temporal_quadrant.R](../workflows/cd8_perturbation/05_plot_la_temporal_quadrant.R) | Supporting | LA descriptive early/late response quadrants |
| [06_plot_hcl_temporal_quadrant.R](../workflows/cd8_perturbation/06_plot_hcl_temporal_quadrant.R) | Supporting | HCl descriptive early/late response quadrants |
| [07_plot_three_condition_pca.R](../workflows/cd8_perturbation/07_plot_three_condition_pca.R) | Display/export | PCA display from existing coordinates |
| [08_plot_three_contrast_program_dotplot.R](../workflows/cd8_perturbation/08_plot_three_contrast_program_dotplot.R) | Display/export | Programme dotplot from supplied NES/P/FDR |
| [09_plot_hypoxia_gsea_curve.R](../workflows/cd8_perturbation/09_plot_hypoxia_gsea_curve.R) | Display/export | Hypoxia running-ES reconstruction; supplied NES/P/FDR retained |
| [10_plot_hypoxia_leading_edge_gene_trajectory.R](../workflows/cd8_perturbation/10_plot_hypoxia_leading_edge_gene_trajectory.R) | Display/export | Leading-edge gene display from a prepared trajectory table |
| [11_plot_hcl_gene_trajectory.R](../workflows/cd8_perturbation/11_plot_hcl_gene_trajectory.R) | Display/export | Paired HCl-minus-NC logCPM display with supplied DE statistics |
| [12_plot_cd8_marker_modules.R](../workflows/cd8_perturbation/12_plot_cd8_marker_modules.R) | Display/export | Marker-expression donor traces and means |

</details>

## Output and implementation details

The [result/output map](result_output_mapping.md) connects complete analytical families to Fig. 1–7 and Supplementary Tables 1–7. Selection for display must not redefine BH families or discard negative/non-estimable results. [Reference implementation notes](../references/README.md) retain the closed historical NA-input check and the documented six/seven-member RTK sensitivity evidence; navigation changes do not alter those implementations.
