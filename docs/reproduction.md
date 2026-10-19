# Reproduction guide

The repository has one entry point, 18 analysis files, five common files and eight configuration templates. It starts from matching prepared inputs. The result map below identifies the retained task aliases; their implementation and required fields are registered once in [stages.tsv](../config/stages.tsv).

## Configure and select a route

Use R 4.5.1 with the required packages listed per stage in [software_dependencies.tsv](software_dependencies.tsv). The listing is a static dependency map, not a lockfile; package versions should follow the Methods and the environment used for the corresponding results. No package is installed automatically.

```sh
Rscript run.R --list
cp config/single_cell.dcf /path/outside/checkout/single_cell.dcf
Rscript run.R --module single_cell --config /path/outside/checkout/single_cell.dcf --plan
Rscript run.R --module pan_pseudobulk --config /path/outside/checkout/single_cell.dcf --plan
Rscript run.R --task single_cell_pseudobulk --steps s06,s07 --config /path/outside/checkout/single_cell.dcf --plan
```

Replace `SET_ME` with real paths. `--plan` reads configuration structure only and accepts incomplete path templates. `--check-config` checks selected-stage filesystem metadata only, without loading study contents or creating analysis outputs. A configured invocation without either option runs the selected stages. Use `--stage task_id` or `--stage task_id.sNN` to restrict a domain; comma-separated selections run in registry order. `--module` accepts either one configuration domain or one analysis-file name. `--task` preserves the former 35 task names; `--steps` is supported with this alias. Legacy task DCF files must be migrated to the corresponding new module template, not passed unchanged.

Every module DCF contains one `Module`/`Output-Root` header and all registered stages for that domain. Keep the stage identifiers and order. `Input-task.name` header fields provide shared paths; stage fields are exactly those in `stages.tsv`. All fields are paths except `edge_id` (`E1` or `E3`). Paths resolve relative to the DCF location. `@input/task.name` uses a header path; `@task.sNN/file` uses a preceding stage output; `@root/task/sNN/file` uses a prepared handoff beneath `Output-Root`. Paths and filenames of the original outputs remain unchanged. References never imply that an upstream analysis has already been run.

Stages run in separate `Rscript` processes. Stage-specific input gates, sample filters, complete-case rules, covariates, comparison universes, seed rules and testing families remain in their named functions. Functions with different behavior remain separate even when they share a former name. Exact matching helpers are shared and rebound to the calling stage scope. Empty independent output directories protect the inputs and existing results. Run logs use `.task_logs/<selection>/`; a failed stage stops the selected route before later stages.

Definition/discovery routes and optional sensitivity routes are selectable separately. Starting from the adopted frozen definitions does not require rerunning discovery. There is no automatic all-project rebuild or implicit replacement of frozen labels.

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

Expression/count matrices, patient metadata, processed Seurat objects, clinical data, provider resources and analysis-prepared registries must match the paper Methods and the fields in [config/stages.tsv](../config/stages.tsv). The supplementary workbooks supply definitions and reported results; they are not drop-in replacements for arbitrary RData/RDS objects or every TSV consumed by the original analyses. See the [supplementary content map](supplementary_inputs.tsv).

Within each module DCF, `@task.sNN/...` references name their producer; `config/stages.tsv` fixes the execution order. Cross-task `@root/task/sNN/...` paths keep the original prepared-output convention and must refer to matching completed results. Selecting a downstream task does not automatically run its producers. Definition construction can be skipped when starting downstream from the adopted definitions. Metadata gate files and stability summaries that are consumed by later steps remain necessary dependencies even when their intermediate records are not separate submitted tables.

Survival uses the supplied `projects` and `clinical_data_cluster_list` objects, preserving the original High/Low group selection. Its preparation writes `tow_group_clinical_data_list.rdata` containing `survival_data_list`. Cox fitting and finalisation retain their covariate/complete-case and FDR rules. Formal KM source data are `km_12panel_plot_data.tsv` and `km_12panel_summary.tsv`; pooled curves use the original all-endpoint complete-case population and are descriptive.

Drug comparisons begin with the matching prepared GLDS-corrected matrices; the original oncoPredict training/prediction and GLDS-correction producers are not included. `predicted_drug_response` writes the complete per-cancer result table plus readiness/audit tables. PRISM is now independent of internal model-performance results: it takes only GDSC fitted dose responses and the two exact PRISM snapshots specified by the `prism_curve` and `prism_cell_metadata` fields (`secondary-screen-dose-response-top8-filtered.csv` and `secondary-screen-cell-line-info.csv`). It produces `gdsc_prism_crossplatform_summary.tsv` and `gdsc_prism_matched_cellline_pairs.tsv`; all five evaluable drugs remain in the BH family although Table 5 presents three.

`drug_shared_axis` takes the complete per-cancer drug result table, the matched shared C1/C4 patient manifest, the per-cancer prediction-file audit and the corrected prediction directory. The shared patient TSV must contain `Cancer`, `Sample15` and `C_subtype`; the audit must contain `Cancer` and `source_file`. Prediction CSVs have sample IDs in the first column and the 286 drug columns. The complete drug table includes `drug`, `Cancer`, `status`, `significant`, `cliffs_delta_C1_minus_C4`, `more_sensitive_group` and `FDR`. Do not replace the shared-patient population with a different manifest merely because its columns match. Output is `strict_drug_shared_axis_residual_by_cancer.tsv`, retaining the original 18-cancer/three-drug comparison and testing families.

`programme_feature_associations` s05 prepares only the numeric components consumed by s06: RNA programmes, signed TCPA pathways, CNV/genomic features, immune estimates and predicted responses. The unreported MOFA pilot/top-variable features and their pilot/definition input requirements have been removed. The original programme construction, reference populations, feature matrices and s06 testing families remain.

The single-cell route begins from the supplied processed objects and frozen annotations, not FASTQ alignment or cell calling. The ICB routes begin from provider expression and clinical resources specified in each configuration. CD8 begins from the prepared count object and sample metadata; FASTQ preprocessing is outside this package, and CD8 data access remains subject to the accession and ethics status in the manuscript.

Analysis source tables and component plotting functions are included where present in the retained scripts. Final manuscript LaTeX assembly and exact editorial panel geometry are maintained with the manuscript; this repository does not claim to regenerate every final composite PDF.

## Code conventions and checks

Each file begins with its purpose and configuration/usage boundary. Each stage records its source at the pre-refactor commit and retains methodological comments within the calculation body. UTF-8 comments are preserved, including Chinese notes. Shared code uses explicit functions, and configuration contains paths rather than a second copy of statistical parameters.

```sh
Rscript checks/run_checks.R /path/outside/checkout/checks.tsv
```

Checks parse the code, verify all module/task routes and schemas, exercise metadata guards with fake files, and test process isolation and failure-stop behavior using synthetic stages. They do not run scientific analyses or establish numerical agreement on the study data. The consolidation was also compared expression by expression with commit `41670141faa1902e623bebcb725ec79e33b0b563`: all 86 calculation bodies and stage-specific metadata checks are preserved after expanding shared helpers and adjusting the package-root binding.

The source map in `config/stages.tsv` retains `legacy_script` only as a Git provenance key, not a current filesystem link. Original implementations are recoverable from that fixed commit and the internal backup. Current distribution contains no legacy `internal/` or `tasks/` implementation tree.
