# AG survival source mapping

Paths below are relative to the project root.

| Public entry | Source and retained section | Adaptation |
|---|---|---|
| 00_prepare_ag_km.R | 3.0survival/KM.R, initialization through survival_data_list save, before duplicate selected-panel assembly | DCF header; isolated loads; configured KM/results path; source significant_projects serialized and session output added |
| 01_run_ag_cox.R | 3.0survival/COX.R, from the Cox function definition through results save, ending before the legacy plotting block | Strict DCF header; isolated named-object loads; configured output; source gene_list argument retained as unused NULL; warning serialization and session output added |
| 02_finalize_ag_cox.R | final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/06_make_fig1D_cox_summary.R, result_cancers through EOF | Configured input/order/output directories; current E summary basename; source complete 22-result family required before final BH; source summary plots retained |
| 02 diagnostic extraction | final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/00_freeze_fig1_source_manifests.R, extract_model through D model-manifest write | Header source_dir aliases configured source_data_dir for source diagnostic writes; d_frozen renamed to newly generated cox_df; current E manifest basename; no B/C/signature handling included |
| 03_export_formal_ag_km.R | 3.0survival/revisions/yanmin_review_r01/v0.1_20260831/scripts/01_export_km_12panel_figure_ready_v0.2.R, expected_cancers through EOF | Configured survival input/output and subordinate qa directory; isolated named-object load; retained source fixed selection, statistical functions, outputs and checks |

Additional method sources: 3.0survival/KM.R (High/Low preparation and original raw-P selection); KM_total_samples.R (seven-field pooled complete-case filter). Current Fig. 1D/E mapping follows fig1_final_main_and_supplement_contract_v1.6_20260923.md, not historical D=Cox labels.

Source statistical bodies are retained as listed, rather than claiming whole-file byte identity. The finalizer's source order file must cover the complete input family; it must not encode only the eleven KM display cancers. Outputs are independent copies and do not replace the study freeze. Current source-script SHA-256 fingerprints and the exact comparison whitelist are included below.

## Current source-text fingerprints and comparison whitelist

These SHA-256 values identify source script text inspected during packaging on 2026-09-30, not historical execution snapshots or data hashes.

| Source script | SHA-256 |
|---|---|
| 3.0survival/KM.R | 8e096ea3e524e5b16d6b359fdd1772d666fb30a85daf9ff5e7c2e598cfdbd918 |
| 3.0survival/COX.R | 2f2b337ca5b62a565a21a65b5c99bdd0d62f41d02c5c595bdc13828ac03fb57c |
| 3.0survival/KM_total_samples.R | cd9deac8394cd72ecbf006a15cf5bb17c3e4ca197ca9daaa60a827c85dd57701 |
| 3.0survival/revisions/yanmin_review_r01/v0.1_20260831/scripts/01_export_km_12panel_figure_ready_v0.2.R | 0a1280c6262509aa19c4c6aacce85e572649872b8eca79c7fdc822f31b1f665b |
| final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/06_make_fig1D_cox_summary.R | 122aeb83748eccf1f39db502285fc3256fd9fe8b115959b8771dc39917027ee0 |
| final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/00_freeze_fig1_source_manifests.R | 93599a18206b435c474d78dbcd0791b3e035112861e3b12ff7307b32c448cd15 |

Static comparison passed for all five extracted sections: 00 whole retained KM block with only configured output paths and appended significant-project/session serialization; 01 whole retained Cox block with only configured save paths and appended warning/session serialization; 02 final-BH/summary/plot block with only E summary basename and input-identity guard; 02 diagnostic extraction with d_frozen-to-cox_df variable and E manifest basename; 03 expected_cancers through EOF byte-identical. The generic DCF headers and isolated loads are outside these retained sections. Whole scripts are not byte-identical to their sources.

The 02 guard requires 22 non-NULL results, unique normalized cancer names and full membership in a duplicate-free frozen order. It does not redefine any model or shorten the BH family. Original race handling, NA handling, error/warning collection, low-event/separation estimates, fixed KM selection and freeze assertions remain unmodified.
