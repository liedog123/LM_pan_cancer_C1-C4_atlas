# Genomics: mutation, scalar endpoints, CNV burden and alteration regions

This scoped entry tests within-cancer C1 versus C4 nonsilent gene-mutation frequencies from an already accepted callable universe. It represents genomic background association, without establishing a cause of expression states.

Run from the package root:

```sh
Rscript workflows/genomics/01_run_mutation_c1c4.R --help
Rscript workflows/genomics/01_run_mutation_c1c4.R --config config/genomics_mutation.example.dcf --check-config
Rscript workflows/genomics/01_run_mutation_c1c4.R --config config/genomics_mutation.example.dcf
```

Replace all placeholders. The strict single-record DCF has exactly callable_dir, mutation_events and output_dir. Relative paths resolve from the DCF directory. Help and configuration checks use base R, read configuration text and file metadata only, load no analysis dependency or research contents, and create no output. Ordinary execution requires data.table. output_dir must be absent or empty, writable and independent of the package, configuration and input directories; direct, symbolic-link and dangling-link boundaries are checked.

## Prepared inputs and upstream boundary

callable_dir must contain four readable nonempty regular TSVs:

| File | Source-inferred content consumed during ordinary execution |
|---|---|
| mutation_callable_decision_v0.1.tsv | accepted_for_formal_denominator; first value must be TRUE |
| mutation_callable_recovery_qa.tsv | pass; all checks must pass |
| mutation_callable_patient_audit.tsv | mc3_callable_exact_sample15, selected_canonical_sample15, Patient_ID, Cancer, C_subtype |
| mutation_callable_formal_eligibility_by_cancer.tsv | Required presence only; source analysis recomputes coverage from the patient audit |

mutation_events is the frozen Xena mutation event table readable by fread, containing sample, gene and effect. Only exact accepted Sample15 entries enter the denominator. Patients must be unique and have C1/C4 states. These schemas derive from code and were not checked against study files.

The upstream source 00b_stage4_mutation_callable_manifest_recovery_v0.1.R is not packaged. Its accepted callable preparation uses all 11,069 official MC3 final tumour-normal pairing records, not the in_broad_batch2/in_broad_batch2a subsets; collapses records to exact Sample15 keys with Cancer/Patient_ID and matches the canonical representative roster. Official manifest integrity and same-version Xena/MC3 v0.2.8 PUBLIC FILTER=PASS row/sample/locus agreement must already have passed. Xena's large deletion reannotation is audited and excluded because it is outside the frozen eleven-class definition. Only accepted callable samples lacking a contract event can be encoded zero. Public upstream source references include the GDC MC3 pairing resource https://api.gdc.cancer.gov/data/752fa601-8e54-4ce1-9184-9e15ea742e78 and public MAF https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc; they do not supply the study-prepared canonical roster or accepted intermediate tables.

Prepared-input release, exact Xena acquisition/version documentation, an approved accession or controlled handoff mechanism and permissions remain unresolved. Example paths are placeholders. This entry does not provide raw-MC3 end-to-end reproduction.

## Preserved analysis and outputs

All eleven source nonsilent classes are retained. Multiple events for one Sample15/gene collapse to a binary event. Both callable groups require at least 30 patients; each tested gene requires at least five mutated patients across the groups. Two-sided Fisher tests produce odds ratios and 95% confidence intervals. BH covers the complete testable gene family within each cancer. Formal significance requires FDR < 0.05 and absolute C1-minus-C4 frequency difference >= 0.10. Cross-cancer summaries count significant directions, preserve conflicting directions and perform no second BH. Full tests, negative results, coverage, summaries and source QA are retained.

Outputs retain all eleven source basenames: mutation_callable_patient_manifest.tsv; mutation_nonsilent_event_class_counts.tsv; mutation_nonsilent_patient_gene_sparse_binary.tsv; mutation_formal_coverage_by_cancer.tsv; mutation_gene_c1_vs_c4_by_cancer.tsv; mutation_gene_cross_cancer_significant_counts.tsv; mutation_single_cancer_significant_gene_summary.tsv; mutation_significant_gene_cancer_units.tsv; mutation_endpoint_summary.tsv; mutation_analysis_qa.tsv; mutation_analysis_sessionInfo.txt. (Eleven generated files; accepted upstream tables stay in callable_dir.) No figure is generated. The source no_figure_output_generated check now inspects output_dir/figures/01_mutation only; it is an isolated-output diagnostic, not frozen-figure QA.

This entry covers the mutation statistical source for Fig. 4 and Supplementary Table 4, not frozen-layout redraw, TMB/MSI or CNV analyses. No source script, accepted study result or frozen figure is modified. Licensing and release remain pending; no license is granted.

Verification scope: R4.5.1 parse/help, synthetic invalid-schema metadata/path checks and retained-source comparisons only. No study schema/content/checksum, analysis dependency, statistical body or figure is read or executed. Numerical equivalence, dependency completeness, prepared-data acquisition and end-to-end reproducibility remain unverified.

Verification on 2026-10-01: one R4.5.1 parse, one help check and 28 synthetic invalid-schema metadata/path cases passed; successful checks created no outputs, and the nonempty-output sentinel was preserved. Full reversed-body comparison, four unchanged extracted helpers and five contiguous scientific-block comparisons passed. Existing entries were not retested. Logs are maintainer-held task artifacts. Study schema/content/checksum, analysis dependencies, statistical bodies and figures were not loaded or executed.

## TMB/MANTIS/MSI continuation (02)

02 continues from prepared Stage 4 readiness tables and retains all within-cancer scalar-endpoint statistics. Run from the package root after replacing the example placeholders:

```sh
Rscript workflows/genomics/02_run_tmb_msi_c1c4.R --help
Rscript workflows/genomics/02_run_tmb_msi_c1c4.R --config config/genomics_tmb_msi.example.dcf --check-config
Rscript workflows/genomics/02_run_tmb_msi_c1c4.R --config config/genomics_tmb_msi.example.dcf
```

The single-record DCF has exactly readiness_dir and output_dir. Relative paths resolve from the DCF directory. Base-R help/check-config reads configuration and path metadata only, loads no study content or analysis dependency and creates no outputs. All three inputs must be readable nonempty regular files. output_dir must be new or empty and writable, independent of input/config/package paths, with direct, symlink, dangling-link and ancestor checks. Ordinary execution requires data.table and ggplot2, and writes results/ and figures/ under this isolated root.

| Required readiness file | Code-inferred ordinary-run content |
|---|---|
| tmb_mantis_c1c4_patient_data.tsv | Patient_ID, Cancer, C_subtype, TMB_Raw, TMB_Log10, MANTIS_Score, MSI_Class |
| readiness_qa.tsv | pass; all checks must pass |
| stage4_layer_execution_status.tsv | Presence only; the source does not read its contents |

Patient_ID must be unique and labels C1/C4. The retained source stops for negative observed TMB_Raw or non-finite observed MANTIS_Score. TMB_Log10 is an upstream calculation, not recomputed here. Readiness derives from canonical conflict exclusion and deterministic representative selection, matched to Thorsson's original Nonsilent Mutation Rate and Bonneville's original MANTIS Score. TMB_Log10 = log10(Nonsilent Mutation Rate + 1); MSI_Class = MSI-H only for MANTIS > 0.4, MSS for <= 0.4, and missing for missing MANTIS. The legacy 20-cancer joined table and mutation-callable denominator are not substitutes. Upstream preparation is not packaged; source-inferred fields are not real-schema validation. Approved intermediate release/acquisition or controlled handoff remains pending; example paths do not supply inputs.

Each continuous endpoint independently excludes its missing observations and requires >=30 patients in each group. Two-sided Wilcoxon (exact=FALSE) reports medians/IQR, Hodges-Lehmann shift/95% CI and Cliff's delta. MSI independently uses nonmissing MSI_Class, requires >=30 per group and >=5 combined MSI-H events, then uses two-sided Fisher with OR/95% CI and frequency difference. Each predefined endpoint uses raw P<0.05; there is no across-cancer or joint-three-endpoint BH. Full coverage/noneligible rows, HOLD statuses, negative results, opposing directions and complete cross-cancer counts are retained.

results/ contains tmb_log10_c1_vs_c4_by_cancer.tsv, mantis_continuous_c1_vs_c4_by_cancer.tsv, msi_binary_c1_vs_c4_by_cancer.tsv, the corresponding tmb/mantis/msi_cross_cancer_significant_counts.tsv, tmb_msi_qa.tsv and sessionInfo.txt. figures/ retains all three source forest plots as PDF, TIFF and preview PNG. These are source diagnostics; they do not recreate frozen Fig. 4B (TMB dumbbell) or Supplementary Fig. 8B (MSI-H bars). MANTIS full-view and source forests are outside those frozen panels. Complete statistics also support Supplementary Table 4.


Genomics 02 verification on 2026-10-01: one R4.5.1 parse, one help check and 30 synthetic invalid-schema metadata/path cases passed. Successful checks created no outputs and preserved the nonempty-output sentinel. The complete 233-line continuous_endpoint-through-EOF body, all eight extracted helpers and five contiguous source blocks compared identically. Existing entries were not retested. Interface/static logs are maintainer-held task artifacts. No study content/schema/checksum, analysis package, statistical body or plot was loaded/executed; real schemas, dependency/device completeness and numerical/end-to-end reproduction remain unverified.

## CNV burden branches (03 and 04)

These are independent input branches, not a 03 → 04 sequence. 03 is a formal C1/C4 continuation from prepared readiness; 04 computes four-state background statistics and the current Fig. 4C source plot from canonical labels and the continuous matrix. 04 does not consume 03 results. Neither entry implements the >0/<0 alteration-region analysis for Fig. 4D.

Run from the package root after replacing all placeholders:

```sh
Rscript workflows/genomics/03_run_cnv_burden_c1c4.R --help
Rscript workflows/genomics/03_run_cnv_burden_c1c4.R --config config/genomics_cnv_burden_c1c4.example.dcf --check-config
Rscript workflows/genomics/03_run_cnv_burden_c1c4.R --config config/genomics_cnv_burden_c1c4.example.dcf
Rscript workflows/genomics/04_run_cnv_burden_four_state.R --help
Rscript workflows/genomics/04_run_cnv_burden_four_state.R --config config/genomics_cnv_burden_four_state.example.dcf --check-config
Rscript workflows/genomics/04_run_cnv_burden_four_state.R --config config/genomics_cnv_burden_four_state.example.dcf
```

03's strict single-record DCF has exactly readiness_dir, continuous_cnv, thresholded_cnv and output_dir. 04's has exactly canonical_labels, continuous_cnv, cancer_order, state_palette and output_dir. Relative paths resolve from the DCF directory. Base-R help/check-config reads only configuration text and path metadata, returns before analysis dependencies/study contents and creates no outputs. Inputs must be readable nonempty regular files; readiness_dir must be readable/searchable. Each output root must be absent or empty and writable, independent of input/config/package paths. Direct, symlink, dangling-link and non-directory-ancestor checks retain the Sys.readlink NA guard.

### C1/C4 readiness continuation (03)

readiness_dir contains cnv_c1c4_sample_manifest.tsv (Patient_ID, Sample15, Cancer, C_subtype) and readiness_qa.tsv (pass; all checks must pass). Patient_ID must be unique and states C1/C4. continuous_cnv and thresholded_cnv are the paired gene-level GISTIC matrices: a Sample gene-ID column and exact Sample15 sample columns. The source checks identical row gene vectors, finite observed continuous values, thresholded domain -2:2 and complete joined burden metrics. These are source-inferred schemas, not real-input validation.

Upstream 00_stage4_readiness_and_input_freeze_v0.1.R excludes patients with discordant canonical Cancer/state, then selects an available exact Sample15 in each matrix using sample-type priority 01,03,06,02,05,07,08,09 and barcode order. It intersects the selected manifests by Patient_ID, stops if selected Sample15 differs between matrices, and checks matrix structure/gene alignment. This CNV readiness manifest is not the mutation-callable universe. Upstream preparation is not packaged; approved matching readiness and matrices must be acquired separately.

The gene-weighted main endpoint is mean_abs_continuous = mean(abs(continuous score)); it is a proxy, not segment/arm/focal burden. Five thresholded sensitivity fractions retain abs>=1, gain>=1, loss<=-1, +2 and -2 definitions. Each metric retains every cancer row, eligibility, negative/opposing and non-estimable results; formal testing requires >=30 C1 and >=30 C4 patients per cancer. Two-sided Wilcoxon exact=FALSE retains Hodges-Lehmann shift/95% CI, Cliff's delta and raw P<0.05, without new BH. Complete cross-cancer counts use the eligible primary endpoint only.

output_dir/results contains cnv_burden_c1c4_patient_data.tsv, cnv_burden_c1_vs_c4_by_cancer.tsv, cnv_burden_cross_cancer_significant_counts.tsv, cnv_burden_qa.tsv and sessionInfo.txt. output_dir/figures retains cnv_mean_abs_continuous_c1_vs_c4_forest as PDF, TIFF and preview PNG, with source size 7.2 x 6.2 inches. This source forest does not reproduce the frozen Supplementary Fig. 8A violin/box layout. Full statistics support Supplementary Table 4.

### Four-state computation and Fig. 4C plot (04)

canonical_labels supplies Sample, Cancer and C_subtype; cancer_order supplies Cancer; state_palette supplies state and color. continuous_cnv has the same continuous gene-level matrix convention. 04 independently excludes patients with discordant Cancer/state, intersects exact Sample15 with available matrix columns, applies the same source sample-type priority/barcode tie-break, and selects one sample per patient. All four states are retained.

mean_abs_continuous is unchanged. Every state must have >=20 patients for a cancer's within-cancer Kruskal-Wallis omnibus test. Counts, eligibility, full statistics and group summaries are retained for all cancers; raw P<0.05 selects cancers for the figure only. No cross-cancer BH or pairwise tests are added. If no eligible cancer is significant at the source selection step, the script retains its stop: already-written complete statistics/group summaries remain, and plotting stops without lowering the threshold.

output_dir/source_data contains fig4_c_four_state_cnv_sample_manifest.tsv, fig4_c_four_state_cnv_patient_data.tsv, fig4_c_four_state_cnv_omnibus_statistics.tsv, fig4_c_four_state_cnv_group_summary.tsv and fig4_c_four_state_cnv_figure_manifest.tsv. output_dir/panel_exports contains fig4_c_four_state_cnv_burden_violin.pdf/.tiff/.png. The source computation, figure selection, basenames, 170 x 42 mm dimensions and drawing code are retained. The generated figure-manifest status is public rerun candidate (unfrozen); this changes output identity only, not analysis. New products do not replace or authenticate frozen assets. C2/C3 provide four-state background association, not a new formal C1/C4 comparison or causal mechanism.

Prepared readiness, canonical labels, matrices, frozen order and palette are not bundled. Author-approved release/acquisition or controlled handoff and permissions remain unresolved; examples do not supply inputs. Code/plot-source preservation has been checked, while frozen-figure numerical identity and rendered appearance have not been tested.

Genomics 03–04 verification on 2026-10-01: two R4.5.1 parses, two help checks and 42 synthetic invalid-schema metadata/path cases (21 per entry) passed. Success created no outputs and nonempty-output sentinels remained unchanged. Two complete source bodies, ten helper blocks (03 eight; 04 two) and ten contiguous scientific blocks passed static comparisons with the documented whitelist reversals. Only source code was hashed. No research content/schema/checksum, analysis package, statistical body or plot was loaded/executed; existing entries were not retested. Real schemas, dependency/device completeness, numerical reproduction and end-to-end execution remain unverified. Interface/static logs are maintainer-held task artifacts.

## CNV alteration-region analysis and Fig. 4D source plot (05)

05 continues from the accepted C1/C4 CNV manifest and thresholded matrix, independently of the statistics from 03 and 04. It retains the complete gene-level analysis, significant-gene mapping, region construction and source drawing code rather than starting from the displayed regions.

```sh
Rscript workflows/genomics/05_run_cnv_alteration_regions.R --help
Rscript workflows/genomics/05_run_cnv_alteration_regions.R --config config/genomics_cnv_regions.example.dcf --check-config
Rscript workflows/genomics/05_run_cnv_alteration_regions.R --config config/genomics_cnv_regions.example.dcf
```

The strict single-record DCF has exactly manifest, thresholded_cnv, gene_coordinates, cytoband and output_dir. Replace all placeholders; relative paths resolve from the DCF directory. All four inputs must be readable nonempty regular files. Base-R help/check-config reads configuration text and path metadata, then returns before analysis dependencies, study contents and output creation. The independent output root must be absent or writable and empty; the inherited 03/04 guard protects input directories, the configuration file and package root, handles symlinks/dangling links and non-directory ancestors, and rejects root/home output targets. A distinct sibling output under a configuration directory is permitted when it does not overlap these protected targets.

| Input | Source-inferred ordinary-run requirements |
|---|---|
| manifest | Patient_ID, Sample15, Cancer, C_subtype; unique patients, C1/C4 only, exact matrix Sample15 matching. This is the prepared paired-CNV readiness manifest, not mutation callable. |
| thresholded_cnv | Sample gene-ID column plus exact Sample15 columns; unique genes, nonmissing integer CNV values in -2:2. |
| gene_coordinates | GENCODE v23 gene probemap with gene, chrom, chromStart and chromEnd; symbols strip the suffix beginning with a vertical bar; only unique-chromosome mappings on chr1–22/X/Y are retained. |
| cytoband | hg38 cytoband table with five columns in order: chrom, chromStart, chromEnd, band, gieStain; source chromosome lengths/layout are computed from it. |

Schemas are inferred from source code, not checked against study data. Accepted manifest preparation is the upstream paired-matrix readiness process described under 03; this entry does not reconstruct or independently accept that preparation.

Amplification is >0 and deletion <0, including low- and high-level calls. The alt_freq >0.01 filter is calculated over the entire manifest-selected matrix before excluding cancers with small groups. Each eligible cancer requires >=30 C1 and >=30 C4 patients. Fisher tests and BH cover all retained genes separately for Amp and Del within each cancer. Formal events require FDR<0.05 and absolute frequency difference >0.15; no cross-cancer BH is added. The original no-eligible-cancer and no-significant-event stops are retained; complete all-gene statistics are written before the latter stop.

Significant events must map at >=90%. Within Cancer × Type × Bias × Chr, Start minus the immediately preceding End >10e6 starts a new cluster. This is an adjacent-gap rule; a cluster's total span can exceed 10 Mb. These are alteration regions/clusters, not GISTIC focal peaks, arm-level events or independent drivers. Within each cancer the source displays the top 15 non-star regions by mean absolute difference plus all genuinely hit star regions; labels retain their source deduplication. The star list, shapes, colour/layout parameters, ggrepel seed 20260809 and 180 × 72 mm export sizes are unchanged.

output_dir/source_data retains fig4_d_cnv_legacy_rule_all_gene_statistics.tsv.gz, fig4_d_cnv_chromosome_map.tsv, fig4_d_cnv_chromosome_labels.tsv, fig4_d_cnv_chromosome_layout.tsv, fig4_d_cnv_legacy_rule_coverage.tsv, fig4_d_cnv_legacy_rule_mapping_summary.tsv, fig4_d_cnv_legacy_rule_qa.tsv, fig4_d_cnv_chromosome_map_manifest.tsv and fig4_d_cnv_legacy_rule_sessionInfo.txt. The complete all-gene statistical family is exported; significant/mapped events and all constructed regions remain source in-memory objects, and the map TSV contains display regions only. No extra region export is introduced. output_dir/panel_exports retains fig4_d_cnv_chromosome_map.pdf/.tiff/.png, with source 600-dpi TIFF and 300-dpi PNG. Two log messages identify the products as an unfrozen public rerun; the source manifest fields and basenames are unchanged. Generated products neither replace nor authenticate frozen Fig. 4D. Complete gene statistics also support Supplementary Table 4.

Matching readiness, thresholded matrix, GENCODE v23 probemap and hg38 cytobands must be acquired separately. Source resource names are recorded here, but approved acquisition/version instructions, intermediate release/controlled handoff, and permissions remain delivery items. Example paths do not provide inputs. Historical scripts and frozen assets remain unchanged.

Genomics 05 verification on 2026-10-01: one R4.5.1 parse, one help check and 33 synthetic invalid-schema metadata/CLI/path cases passed. Successful checks created no output directories and the nonempty-output sentinel was unchanged. The complete 226-line MIN_N-through-EOF body and seven contiguous blocks compared identically after reversing only two output-identity messages. Existing entries were not retested. No study content/schema/checksum, analysis package, statistical body or figure was loaded or executed; real schemas, dependency/device completeness, numerical equivalence and end-to-end reproduction remain unverified. Logs are maintainer-held task artifacts.
