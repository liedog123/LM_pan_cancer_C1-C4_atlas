# Genomics source mapping

Source paths relative to the originating project:

- 8.0pancancer_c1c4_bulk_downstream_rebuild/04_genomic_instability/scripts/01a_stage4_mutation_c1c4_v0.1.R (324 lines): full body from nonsilent_classes through EOF.
- The same scripts directory's stage4_common_v0.1.R: only log_msg, assert_nonempty, write_tsv and save_session_info are extracted unchanged.
- docs/stage4_mutation_c1c4_analysis_contract_v0.1_20260730.md and docs/stage4_mutation_c1c4_execution_report_v0.1_20260730.md in that stage: accepted scientific and execution boundary.
- scripts/00b_stage4_mutation_callable_manifest_recovery_v0.1.R in that stage: upstream preparation only, not engineered here.

Adaptations are limited to a strict base-R DCF/CLI metadata guard replacing the private source/header; callable_dir for the four accepted input-table paths formerly under out_dir; mutation_events for xena_path; output_dir for all outputs; and the source no-figure diagnostic redirected from stage4_root/figures/01_mutation to output_dir/figures/01_mutation. The full statistical body is preserved after reversing those five path substitutions. No event-class, denominator, sample rule, test, threshold, BH family, summary, negative-result or stop condition is changed. The metadata header uses the existing transcriptome_proteome 14 interface with its Sys.readlink NA guard, adds callable-file metadata and protects the package root.

The upstream official all-record pairing, exact canonical Sample15, Xena/MC3 PASS identity and large-deletion audit must already be accepted. Prepared inputs are not bundled or publicly obtainable through this package yet. The ordinary source's first-row decision semantics and all recovery-QA acceptance semantics are preserved, rather than strengthened by a new scientific rule. No raw-MC3 end-to-end or frozen Fig. 4 reproduction is claimed.

Verification on 2026-10-01: one R4.5.1 parse, one help check and 28 synthetic invalid-schema metadata/path cases passed; successful checks created no outputs, and the nonempty-output sentinel was preserved. Full reversed-body comparison, four unchanged extracted helpers and five contiguous scientific-block comparisons passed. Existing entries were not retested. Logs are maintainer-held task artifacts. Study schema/content/checksum, analysis dependencies, statistical bodies and figures were not loaded or executed.

Current source-script SHA-256 (code only):

01a_stage4_mutation_c1c4_v0.1.R: 1e750a3e4f661ab9a6392961223c9abb7df1103d74fd2485fd19ef62e1e01a54

stage4_common_v0.1.R: 02212c65afb803ced056fff8e911ea08ace45952a36c1ca43a69d359a1f52334

## TMB/MANTIS/MSI (02)

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/04_genomic_instability/scripts/01_stage4_tmb_msi_c1c4_v0.1.R. The complete continuous_endpoint-through-EOF body is retained byte-for-byte. Only eight unchanged stage4_common_v0.1.R helpers are extracted: log_msg, assert_nonempty, assert_cols, write_tsv, save_session_info, wilcox_summary, summarize_cross_cancer and save_figure_pair. The private source/header is replaced by the strict base-R CLI/DCF metadata guard, readiness_dir, and isolated output_dir/results plus output_dir/figures. No statistical-body path reversal is required.

Upstream boundary: scripts/00_stage4_readiness_and_input_freeze_v0.1.R (TMB/MANTIS construction) and the common canonical-roster helpers. Scientific source references: docs/stage4_genomic_instability_c1c4_analysis_contract_v0.1_20260730.md (TMB/MANTIS/MSI clauses) and docs/stage4_c1c4_execution_report_v0.1_20260730.md in the same stage. The old mutation HOLD in the original general contract has been superseded by the accepted mutation-specific contract and is not restored. Fig. 4's current README maps these statistics to frozen Fig. 4B and Supplementary Fig. 8B; all three retained source plots are diagnostics, not frozen-layout exports. The original source scripts, accepted results and frozen assets remain unchanged.


Genomics 02 verification on 2026-10-01: one R4.5.1 parse, one help check and 30 synthetic invalid-schema metadata/path cases passed. Successful checks created no outputs and preserved the nonempty-output sentinel. The complete 233-line continuous_endpoint-through-EOF body, all eight extracted helpers and five contiguous source blocks compared identically. Existing entries were not retested. Interface/static logs are maintainer-held task artifacts. No study content/schema/checksum, analysis package, statistical body or plot was loaded/executed; real schemas, dependency/device completeness and numerical/end-to-end reproduction remain unverified.

Current 02 source-script SHA-256 (code only): 01_stage4_tmb_msi_c1c4_v0.1.R: d5671207d365741c1500f2a2fa0ec7af1c3aac60055f9288ac0b91b3e6917b4d. Common-helper source SHA-256 remains 02212c65afb803ced056fff8e911ea08ace45952a36c1ca43a69d359a1f52334.

## CNV burden (03–04)

03 source: 8.0pancancer_c1c4_bulk_downstream_rebuild/04_genomic_instability/scripts/02_stage4_cnv_burden_c1c4_v0.1.R. Its complete 151-line tryCatch-through-EOF body is byte-identical. Eight unchanged helpers are extracted from stage4_common_v0.1.R: log_msg, assert_nonempty, assert_cols, write_tsv, save_session_info, wilcox_summary, summarize_cross_cancer and save_figure_pair. Adaptations replace the private header with strict base-R CLI/DCF metadata checks; readiness_dir and paired matrices come from configuration, with isolated output_dir/results and output_dir/figures. No body reversal is required.

03 upstream source: scripts/00_stage4_readiness_and_input_freeze_v0.1.R (CNV section around lines 140–195) and the common canonical-roster/matrix-selection helpers. The actual paired-matrix representative selection/alignment is an upstream prerequisite, not MC3 callable recovery. Scientific definitions follow docs/stage4_genomic_instability_c1c4_analysis_contract_v0.1_20260730.md, CNV burden clauses; the old high-level-event analysis is not added.

04 source: final_figures/c1c4_atlas/fig4_genomic_instability/panel_sources/plot_fig4_c_four_state_cnv_burden.R. The complete 314-line options-through-EOF body matches after reversing only the main input/output path block, script_name to the new public basename, and the figure-manifest status from frozen formal panel to public rerun candidate (unfrozen). log_event and sample_type_priority are unchanged. Full all-cancer counts/statistics/group summaries, original significance-only figure selection, zero-significant-cancer stop, complete source drawing code, file basenames and dimensions remain intact. The status adaptation identifies the new output as a candidate; it does not alter the statistical analysis or promote a freeze.

04 authority: the same Fig. 4 directory's README.md and docs/fig4_main_and_supplement_design_contract_v1.0_20260809.md, panel C. This is a separate canonical-plus-matrix computation branch, not a consumer of 03. C1/C4 Supplementary Fig. 8A layout export remains distinct from 03's retained forest. Fig. 4D >0/<0 alteration regions remain outside this increment. Original scripts, Methods, accepted statistics and frozen assets are unchanged.

Genomics 03–04 verification on 2026-10-01: two R4.5.1 parses, two help checks and 42 synthetic invalid-schema metadata/path cases (21 per entry) passed. Success created no outputs and nonempty-output sentinels remained unchanged. Two complete source bodies, ten helper blocks (03 eight; 04 two) and ten contiguous scientific blocks passed static comparisons with the documented whitelist reversals. Only source code was hashed. No research content/schema/checksum, analysis package, statistical body or plot was loaded/executed; existing entries were not retested. Real schemas, dependency/device completeness, numerical reproduction and end-to-end execution remain unverified. Interface/static logs are maintainer-held task artifacts.

Current CNV source-script SHA-256 (code only):

02_stage4_cnv_burden_c1c4_v0.1.R: 26195bb0149246014267475c1a531d945376334e9bb721f1c1b458bd067a59a8

plot_fig4_c_four_state_cnv_burden.R: 67a8b103bcb2aa81de3080a62dacb8f703bd30005b69837f6533eab5e35f0441

## CNV alteration regions and Fig. 4D (05)

Source: final_figures/c1c4_atlas/fig4_genomic_instability/panel_sources/plot_fig4_d_cnv_chromosome_map.R; authority: the same Figure directory's docs/fig4_main_and_supplement_design_contract_v1.0_20260809.md, panel D and its dedicated boundary clauses. Source-script SHA-256 (code only): 45595ea3a424f552b3da17f6031ed89167bfb0e558be502b83ad063a3fef2e31.

The private header is replaced by the inherited 04 strict base-R CLI/DCF/path guard, four configured input paths and isolated output_dir/source_data plus output_dir/panel_exports. fields are manifest, thresholded_cnv, gene_coordinates, cytoband and output_dir. No study data, schema or research checksum was read. The complete 226-line MIN_N-through-EOF scientific/drawing/export body is identical after reversing only these exact message substitutions:

- [1/7] Read frozen C1/C4 manifest and selected CNV matrix → [1/7] Read C1/C4 manifest and selected CNV matrix for an unfrozen public rerun
- [7/7] Freeze source data, QA and exports → [7/7] Write unfrozen public rerun source data, QA and exports

All thresholds, full-manifest frequency filtering, complete cancer-wise Amp/Del BH families, coordinate rules, stops, star/Top-15 display selection, seed, dimensions, basenames and source manifest fields are retained. Clustering uses the immediately preceding gene's End in Start-shift(End)>10e6; it does not cap total region span at 10 Mb. The complete all-gene statistical output precedes significant/display selection. All regions are constructed in memory; only display regions are exported by the source, and no extra export is added. The rerun uses the prepared paired-CNV manifest, not the mutation-callable denominator or 03/04 statistical outputs. Acquisition and release/permissions remain unresolved. New products are unfrozen candidates; original source scripts and frozen assets are unchanged.

Genomics 05 verification on 2026-10-01: one R4.5.1 parse, one help check and 33 synthetic invalid-schema metadata/CLI/path cases passed. Successful checks created no output directories and the nonempty-output sentinel was unchanged. The complete 226-line MIN_N-through-EOF body and seven contiguous blocks compared identically after reversing only two output-identity messages. Existing entries were not retested. No study content/schema/checksum, analysis package, statistical body or figure was loaded or executed; real schemas, dependency/device completeness, numerical equivalence and end-to-end reproduction remain unverified. Logs are maintainer-held task artifacts.
