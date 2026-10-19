# Reproduction routes

Select the result family below, then copy the corresponding task DCF from the [task index](task_index.tsv). `--plan` inspects structure only. Configure matching inputs before ordinary execution. There is no global run-all command: cross-task prepared inputs are explicit, and absent inputs are not inferred or generated automatically.

## Fig. 1; Supplementary Table 1

| Task | Role |
|---|---|
| `acid_signature_selection` | Adopted multi-stage nine-gene selection; not a replaceable generic DEG call. |
| `signature_correlations` | Signature and tumour-programme correlations reported in Supplementary Fig. 1. |
| `acid_ag_groups` | Acid scores, cancer-wise AG order and processed-input handoff define downstream comparisons. |
| `acid_ph_correlations` | Fixed pH transcript panel and correlation scope. |
| `acid_distributions` | Reported cancer-level and tumour-normal score comparisons. |
| `ag_survival` | Paper-specific AG endpoints, complete-case rules, covariates and cross-cancer FDR; formal KM exports. |

## Fig. 2; Supplementary Table 2

| Task | Role |
|---|---|
| `tcga_immune_discovery` | AG-associated discovery supplies the pathway-definition route. |
| `lm12_definition` | Cross-cancer voting, stable genes and adopted family definitions. |
| `c1c4_state_assignment` | Ordered axes and canonical sample/patient labels are central definitions. |
| `external_bulk_enrichment` | External cohort enrichment is a reported supporting result. |
| `c1c4_clinical` | C1/C4 patient eligibility, clinical models and survival meta-analysis. |
| `immune_ecology` | Within-cancer and cross-cancer immune/stromal tests with fixed state and testing scopes. |

## Fig. 3; Supplementary Table 3

| Task | Role |
|---|---|
| `mrna_programmes` | Reported recurrent mRNAs and their functional representative selections. |
| `lncrna_rppa` | Reported lncRNA recurrence, linked mRNA functions and RPPA features. |
| `boundary_gsea` | Reported fixed adjacent-state contrasts and full testing families. |
| `mrna_rppa_concordance` | RNA-protein comparison remains in Supplementary Table 3. |
| `subtype_context` | Immune-subtype and PAM50 adjustment remains in Supplementary Fig. 5 and Table 3. |

## Fig. 4; Supplementary Table 4

| Task | Role |
|---|---|
| `somatic_mutations` | Reported patient-level mutation selection and Fisher/BH comparisons. |
| `tmb_msi` | Reported burden/MSI transformations, thresholds and eligibility. |
| `cnv_burden` | Reported two-state and four-state CNV burden comparisons. |
| `cnv_regions` | Gene-level copy-number tests and descriptive region aggregation. |

## Fig. 5; Supplementary Table 5

| Task | Role |
|---|---|
| `predicted_drug_response` | Reported C1/C4 inference across all evaluable drugs within each cancer. |
| `drug_crossplatform` | Reported three drugs retain the original five-drug PRISM correction family. |
| `drug_shared_axis` | Reported leave-one-drug-out PC1 correction for three drugs across 18 cancers. |
| `programme_feature_associations` | Reported Fig. 5C associations require complete within-cancer testing families. |

## Fig. 6 and ICB; Supplementary Table 6

| Task | Role |
|---|---|
| `single_cell_states` | Reported AUCell refinement, compartment labels and sample-level voting. |
| `single_cell_pseudobulk` | Count/model eligibility and stability evidence feed the reported programme analyses. |
| `bassez_sensitivity` | Reported study-exclusion and marker sensitivity remains in Supplementary Fig. 10. |
| `melanoma_pseudobulk` | Patient-level CD8 and other cell-type programme comparisons remain reported. |
| `icb_three_cohort` | Outcome-blind state assignment and primary response comparisons are reported. |
| `liu_icb` | Separate post hoc and timing sensitivity remains reported in Table 6. |
| `gide_c1` | Exploratory within-C1 programme-response result remains in Supplementary Fig. 9. |

## Fig. 7; Supplementary Table 7

| Task | Role |
|---|---|
| `cd8_qc` | Sample metadata, filtered expression and PCA support reported figures. |
| `cd8_response` | Donor-adjusted contrasts, programme testing and HCl comparisons support main conclusions. |
| `cd8_figures` | Also produces reported quadrant and trajectory tables; not dispensable preview-only code. |

## Starting boundaries and handoffs

Expression/count matrices, patient metadata, processed Seurat objects, clinical data, provider resources and analysis-prepared registries must match the paper Methods and the fields in [step_index.tsv](step_index.tsv). The supplementary workbooks supply definitions and reported results; they are not drop-in replacements for arbitrary RData/RDS objects or every TSV consumed by the original analyses. See the [supplementary content map](supplementary_inputs.tsv).

Within a task, the [output-link index](output_links.tsv) maps `@sNN/...` producers to consumers. Across tasks, matching prepared outputs must be configured explicitly. Definition construction can be skipped when starting downstream from the adopted definitions. Metadata gate files and stability summaries that are consumed by later steps remain necessary dependencies even when their intermediate records are not separate submitted tables.

Survival uses the supplied `projects` and `clinical_data_cluster_list` objects, preserving the original High/Low group selection. Its preparation writes `tow_group_clinical_data_list.rdata` containing `survival_data_list`. Cox fitting and finalisation retain their covariate/complete-case and FDR rules. Formal KM source data are `km_12panel_plot_data.tsv` and `km_12panel_summary.tsv`; pooled curves use the original all-endpoint complete-case population and are descriptive.

Drug comparisons begin with the matching prepared GLDS-corrected matrices; the original oncoPredict training/prediction and GLDS-correction producers are not included. `predicted_drug_response` writes the complete per-cancer result table plus readiness/audit tables. PRISM is now independent of internal model-performance results: it takes only GDSC fitted dose responses and the two exact PRISM snapshots named in its DCF. It produces `gdsc_prism_crossplatform_summary.tsv` and `gdsc_prism_matched_cellline_pairs.tsv`; all five evaluable drugs remain in the BH family although Table 5 presents three.

`drug_shared_axis` takes the complete per-cancer drug result table, the matched shared C1/C4 patient manifest, the per-cancer prediction-file audit and the corrected prediction directory. The shared patient TSV must contain `Cancer`, `Sample15` and `C_subtype`; the audit must contain `Cancer` and `source_file`. Prediction CSVs have sample IDs in the first column and the 286 drug columns. The complete drug table includes `drug`, `Cancer`, `status`, `significant`, `cliffs_delta_C1_minus_C4`, `more_sensitive_group` and `FDR`. Do not replace the shared-patient population with a different manifest merely because its columns match. Output is `strict_drug_shared_axis_residual_by_cancer.tsv`, retaining the original 18-cancer/three-drug comparison and testing families.

`programme_feature_associations` s05 prepares only the numeric components consumed by s06: RNA programmes, signed TCPA pathways, CNV/genomic features, immune estimates and predicted responses. The unreported MOFA pilot/top-variable features and their pilot/definition input requirements have been removed. The original programme construction, reference populations, feature matrices and s06 testing families remain.

The single-cell route begins from the supplied processed objects and frozen annotations, not FASTQ alignment or cell calling. The ICB routes begin from provider expression and clinical resources specified in each configuration. CD8 begins from the prepared count object and sample metadata; FASTQ preprocessing is outside this package, and CD8 data access remains subject to the accession and ethics status in the manuscript.

Analysis source tables and component plotting functions are included where present in the retained scripts. Final manuscript LaTeX assembly and exact editorial panel geometry are maintained with the manuscript; this repository does not claim to regenerate every final composite PDF.

## Script conventions

R source files use UTF-8 and a common opening comment block: purpose, task/step, DCF input fields, output directory boundary, example invocation and documentation links. Task purpose is given in English for public readers. Existing internal comments, including Chinese methodological notes, are retained; comment translation is not a change in the statistical method. Configurations are path templates. All ordinary analysis outputs must go to independent empty directories.

## Verification boundary

Source comparison preserves the retained scientific bodies except the four documented extractions: survival input preparation, PRISM concordance, shared-axis adjustment and programme-feature numeric preparation. Tests in [checks](../checks/README.md) exercise syntax, task wiring and metadata guards. They do not establish study-input schemas or numerical agreement with frozen results. No original patient data, result package or new participant-level table is bundled here.
