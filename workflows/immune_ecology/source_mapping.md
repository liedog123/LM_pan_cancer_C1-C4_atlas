# Stage 2 source mapping

Paths below are relative to the originating project; fingerprints identify current source text, not recovered historical execution hashes.

| Public entry | Source script | Source SHA-256 |
|---|---|---|
| 01_run_stage2_immune_ecology.R | 8.0pancancer_c1c4_bulk_downstream_rebuild/02_immune_ecology/scripts/01_run_stage2_immune_ecology_v0.1_20260723.R | 6ee64dc33ef47ad5cb6b20ed70934a053d80d8244085e3501da70ead475a52b5 |
| 02_finalize_stage2_evidence.R | 8.0pancancer_c1c4_bulk_downstream_rebuild/02_immune_ecology/scripts/02_finalize_stage2_evidence_and_gate_v0.1_20260723.R | 1140a57dea22c0170e642723fed621d1ed578533062d7afb8849079b1d03d9bb |

01 text from log_event <- function through EOF is byte-identical. The replacement execution header supplies strict DCF/metadata/output protection, loads original packages after the check early exit, changes script_name, sets smoke_mode=FALSE and selected_cancers=character(), and resolves the seven input files and isolated output directory. No statistical function, threshold, feature, design or family changes.

02 text from log_event <- function through EOF is identical after two explicit path substitutions: read_result reads input_result_dir instead of result_dir; the START log reports input_result_dir. write_result and all final artifacts continue to use result_dir below the new output_dir. The replacement header resolves the eleven expected input basenames, applies the same output protection, loads original packages after the early exit, and assigns separate result/figure/doc directories. No evidence rule, sensitivity merge, Gate logic, report content or plotting body changes. The original manifest collects only existing newly written finalization artifacts; there is no core-table completeness check requiring a shared writable results directory.

Formal authority: final_figures/c1c4_atlas/fig2_c1c4_state_immune_ecology/docs/fig2_final_main_and_supplement_contract_v1.8_20260907.md.

| Formal result | Existing panel/statistical consumer | Stage 2 inputs |
|---|---|---|
| Fig. 2E, CIBERSORT C1/C4 | final_figures/c1c4_atlas/fig2_c1c4_state_immune_ecology/panel_sources/05_prepare_and_plot_fig2E_cibersort_c1_vs_c4_immune_ecology.R | immune_pairwise_effects_by_cancer and immune_evidence_classification |
| Supplementary Fig. 4A, MCPcounter C1/C4 | final_figures/c1c4_atlas/fig2_c1c4_state_immune_ecology/panel_sources/06_prepare_and_plot_fig2S2A_mcpcounter_c1_vs_c4_immune_ecology.R | same two complete tables |
| Supplementary Fig. 4B, ESTIMATE state distributions | final_figures/c1c4_atlas/fig2_c1c4_state_immune_ecology/panel_sources/10_prepare_and_plot_fig2S2B_estimate_c1c4_boxplots.R; current display continuation 30_render_fig2S2B_estimate_meta_removed_candidate_v0.1_20260902.R | stage2_sample_analysis_manifest, reused ESTIMATE scores, pairwise statistics and evidence classification |

These mappings were read from source code/current authority, without opening research inputs or frozen source-data values. The historical panel scripts expect one original results directory; they are not new public entries and require future directed path adaptation to consume separate public core/evidence outputs. Frozen figure and source-data files remain untouched. Supplementary Table 2 export and publication panel redraw remain delivery tasks; the analytical families are retained here in full.
