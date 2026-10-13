# Source mapping

Code hashes only; no study inputs or frozen outputs were opened or hashed.

## 01_run_mrna_models.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/35_run_stage3_mrna_ovr_c1c4_models_v0.2_20260729.R. SHA-256: 22bf351457029399549b23dc6dd62a202143d290f499e36bc098bdfdb71c0b13; lines 1–537.

log_event/helpers and complete run_design/run_cancer, with main statistical and QA/output bodies retained.

Finite adaptations: prepended strict DCF metadata early exit; removed source CLI mode parsing/private roots; forced complete formal scope and separate output root; changed script identity. No numerical transformation, filter, statistical family or threshold change.

Five-design output fields retain historical source labels STATE_PORTRAIT_MAIN (OVR) and JOINT_ENDPOINT_SUPPORT (C1/C4). Current authority is C1/C4 formal; OVR historical support only. These labels do not define current Figure eligibility.

## 02_run_mrna_recurrence_ora.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/36_build_stage3_mrna_ovr_c1c4_recurrence_ora_v0.2_20260729.R. SHA-256: 990c61f773eae397266fb82ed4c43c8066e1bccc9ba7338998224c83dfcfd6b3; lines 1–702.

write_tsv through compute_ora and main recurrence/ORA/QA bodies retained.

Finite adaptations: prepended strict DCF metadata early exit; removed source CLI mode parsing/private roots; forced complete formal scope and separate output root; changed script identity. No numerical transformation, filter, statistical family or threshold change.

Manifest Output paths rebound by basename to model_dir after exact 22-file family validation; GMT paths come from DCF. Final generated report sentence reflects current C1/C4 authority.

Five-design output fields retain historical source labels STATE_PORTRAIT_MAIN (OVR) and JOINT_ENDPOINT_SUPPORT (C1/C4). Current authority is C1/C4 formal; OVR historical support only. These labels do not define current Figure eligibility.

## 03_finalize_rppa_recurrence_tcpa.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/46_build_stage3_rppa_c1c4_tcpa_source_v0.1_20260730.R. SHA-256: 217636e07b48f29b013ecfd59d349602bfe90dcda8efb6206ec23bc47daee36e; lines 1–613.

expected_cancers through EOF retained verbatim.

Finite adaptations: prepended strict DCF metadata early exit; removed source CLI mode parsing/private roots; forced complete formal scope and separate output root; changed script identity. No numerical transformation, filter, statistical family or threshold change.

RPPA model fitting is upstream of 03 and is now provided by 04. Full-panel per-cancer BH replaces old class-specific FDR; TCPA membership is annotation, not ORA.

## Verification performed

R4.5.1 (2025-06-13): three parse checks, three help checks and twelve synthetic metadata/early-exit/output-isolation checks passed. File-field directory, input/output overlap and nonempty-output cases were rejected; positive checks created no output directories. Seven retained source-block comparisons passed (35 helper segment excludes removed parse_args). The final 02 basename correction was separately reparsed/help/metadata tested. Source report final narrative is explicitly adapted; no study input or frozen results were read.

01_run_mrna_models.R source lines 13–308 retained-block comparison PASS.
01_run_mrna_models.R source lines 358–459 retained-block comparison PASS.
01_run_mrna_models.R source lines 460–538 retained-block comparison PASS.
02_run_mrna_recurrence_ora.R source lines 38–164 retained-block comparison PASS.
02_run_mrna_recurrence_ora.R source lines 218–397 retained-block comparison PASS.
02_run_mrna_recurrence_ora.R source lines 442–619 retained-block comparison PASS.
03_finalize_rppa_recurrence_tcpa.R source lines 59–614 retained-block comparison PASS.

The report authority adaptation replaces source 36 line 658 (the final next-step narrative sentence). Historical Evidence_Role fields remain unchanged. Checks establish source-block preservation and metadata interfaces only, not numerical equivalence.


## 04_run_lncrna_rppa_models.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/01_run_stage3_feature_models_v0.1_20260723.R. SHA-256: 706a25f2aebe177741d7547143176768e443d664d3da5157b460d953e6a897ae; 490 lines.

Source log_event/norm helpers (lines 13–28), design_specs (40–59), and align_matrix_to_patient through run_one_model (61–243) are retained verbatim. The complete lncRNA/RPPA preparation and model calls (390–445) are retained verbatim. Historical branches inside these helpers are retained but not executed.

Finite adaptations: prepend the strict DCF regular-file/directory metadata interface, including all 22 lncRNA Rdata basenames; remove old parse_args/private roots, mRNA loading/execution and smoke/subset selectors; select C1_vs_C4 only, REPRESENTATIVE_MAIN only and UNADJUSTED only in run_layer_models; preserve its post-fit masking, class-specific BH, schema, basename and coverage body. Require complete 22-cancer patient scope. Adapt source scope logging/QA from three layers/seven designs to two layers/one formal branch. Use an independent empty output root with temporary sibling and atomic commit; rebind layer-manifest Output to final paths before writing. No feature filter, sample normalization/alignment, empirical-Bayes background, model calculation, threshold, transform or source BH family changes. RPPA's source Feature_Class FDR remains intermediate and is finalized by public 03.

Only source text/metadata and synthetic invalid-schema files are used for packaging checks. Real input schema, model numerical equivalence, dependency completeness and end-to-end execution remain unverified. lncRNA display/export continuations remain outside 04; 05–08 are documented below.


04 verification performed on 2026-09-30: R 4.5.1 (2025-06-13) parse and help passed. Six synthetic metadata cases passed: complete 22-file invalid-schema early exit; rejection of a directory in a file field, input/output overlap, symlink overlap, nonempty output and incomplete 22-file scope. Assertions confirmed that no new output directory was created and an existing sentinel remained unchanged. Four retained source blocks and the run_layer scope-only reverse comparison passed. Source-to-public basename/schema/final-manifest-path compatibility is a static inference, not a tested research-data handoff. No study object, statistical manifest, gene-set member or research checksum was read; no analysis body or existing 01–03 entry was executed.

## 05_freeze_lncrna_recurrence.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/40_freeze_stage3_lncrna_c1c4_recurrence_v0.1_20260730.R. SHA-256: 49b808060e45751fecd21966e4a56298ac69cb560b55d3279f3bd6629a82ef90; 236 lines.

Finite common adaptations: strict DCF help/metadata early exit, regular-file/directory validation, symlink resolution, independent empty output root, ordinary-only parent creation, source temporary-sibling atomic commit, and public script identity (06–08). Original sources and old public 01–04 are unchanged.

Private input/output roots are replaced with model_dir and an independent output_dir/01_lncrna_recurrence. The complete read_one-through-EOF body is retained; output_dir parent creation occurs only during ordinary execution.

## 06_run_lncrna_gene_coexpression.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/41_run_stage3_lncrna_c1c4_gene_coexpression_v0.1_20260730.R. SHA-256: 08b41d8e6b8d7f39911267b445d3f6d2675e8193ccf192c527105c5543c534c9; 493 lines.

Finite common adaptations: strict DCF help/metadata early exit, regular-file/directory validation, symlink resolution, independent empty output root, ordinary-only parent creation, source temporary-sibling atomic commit, and public script identity (06–08). Original sources and old public 01–04 are unchanged.

Private roots and the source paths list are replaced by the five DCF inputs and output_dir/01_gene_level_coexpression. All original helper calculations are retained; complete formal 22-cancer identity is checked. The output manifest uses literal-prefix substring relocation rather than regular-expression substitution, preserving paths containing regex metacharacters. No serialization, candidate/filter, patient alignment, statistical or output-basename change.

## 07_build_lncrna_edge_recurrence.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/42_build_stage3_lncrna_c1c4_edge_recurrence_v0.1_20260730.R. SHA-256: e4869a0cc97231df30526fd720a22e9a35d075f1982ede714939a31ae16bc94a; 449 lines.

Finite common adaptations: strict DCF help/metadata early exit, regular-file/directory validation, symlink resolution, independent empty output root, ordinary-only parent creation, source temporary-sibling atomic commit, and public script identity (06–08). Original sources and old public 01–04 are unchanged.

Private roots are replaced by coexpression_dir, recurrent_candidates and output_dir/02_edge_recurrence. After loading the manifest, complete 22-cancer identity and three Cancer-specific basenames are validated before paths are rebound inside coexpression_dir. The lnc_availability-through-EOF statistical/output body is retained verbatim. No recurrence/background/threshold or serialization change.

## 08_run_lncrna_linked_mrna_ora.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/43_run_stage3_lncrna_c1c4_linked_mrna_ora_v0.1_20260730.R. SHA-256: d7507ac121f2cbcf3a42d4c50a4af41ed9de36ed9a016507c5c9c592eeb54246; 481 lines.

Finite common adaptations: strict DCF help/metadata early exit, regular-file/directory validation, symlink resolution, independent empty output root, ordinary-only parent creation, source temporary-sibling atomic commit, and public script identity (06–08). Original sources and old public 01–04 are unchanged.

Private roots are replaced by edge_recurrence_dir and output_dir/03_linked_mrna_ora; five resource paths come from DCF. read_gmt helper and Version-through-EOF statistical/output body are retained verbatim. Collection versions/evidence roles, all terms/families, resource serialization and ordinary-run MD5 recording remain unchanged.

05–08 verification on 2026-09-30:
R4.5.1 (2025-06-13): four parse and four help checks; 28 synthetic invalid-schema metadata/safety cases; seven retained-source-block comparisons; a literal-prefix path-adapter example. No analysis dependency or study content was loaded, no study file was hashed and no statistical body or figure ran. Actual schemas, dependencies, numerical equivalence and end-to-end reproduction remain unverified.
Source-block comparisons cover 40 read_one–EOF; 41 norm_sample–before main, annotation–before path_columns, output_paths–EOF; 42 lnc_availability–EOF; 43 read_gmt–before main and Version–EOF. Check scope is source preservation and metadata safety, not statistical equivalence. Source-required checksums and frozen QA assertions remain executable source logic, not newly verified facts.

Check stdout: lncrna_pack_checks_20260930_3a8132065a101.log (maintainer-held internal verification artifact, not distributed; not a required runtime input).
Independent root-agent source-diff review additionally confirmed: 05's second expected_cancers-through-EOF body is unchanged; 06's original log_event/helpers (8,340 characters) and 07's log_event helper (513 characters) are unchanged; remaining differences in 06–08 are limited to the listed CLI, paths/output, 22-cancer identity and manifest IO adaptations. This is static code evidence, not research-data execution.

## 09–12 local boundary statistical continuations

Common finite adaptations: strict DCF CLI and metadata early exit; symlink-resolved input/output isolation; complete fixed model family; independent empty output root; retained source temporary-sibling/atomic commit and serialization. No statistical transformation, filter, empirical-Bayes background, grouping, threshold, seed/fallback, BH family or vote rule is changed.
09 removes source smoke/cancer selectors and fixes formal ALL with a complete-22-cancer patient guard. Paths are supplied by DCF and final_dir is output_dir/01_within_cancer_models. Source edge models and full-family statistical outputs remain retained.
10/11 replace private roots/resource paths; after ordinary manifest loading, complete unique Cancer identity and its fixed basename are validated, then Output is rebound into model_dir/individual. No study-manifest content is read during metadata checks. 11 moves the source E1/E3 CLI choice to the strict edge_id enum; each edge has its own output root.
12 supplies the three edge directories and definition registry explicitly, including the source all-edge machine_qa requirement. It consumes provided E2 results and retains the original Reuse_Status/report labels. Those labels do not prove the supplied file is the historical frozen object. Historical pending eligibility/report dates and promotion text are source-stage labels; current Fig. 3B is frozen and its 16 display rows belong to the separate panel consumer.

SOURCE 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/01_run_fig3_mrna_boundary_models_v0.1_20260806.R
SHA256 a7145cad2d3c890fb525ee5d2ef833a87f869e79067c59da869f152fe7dbd675  8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/01_run_fig3_mrna_boundary_models_v0.1_20260806.R
BLOCK_PASS 09_run_mrna_boundary_models.R source=35-197 public=86-248
BLOCK_PASS 09_run_mrna_boundary_models.R source=224-245 public=262-283
BLOCK_PASS 09_run_mrna_boundary_models.R source=246-325 public=285-364
SOURCE 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/03_e2_full_rank_gsea/v0.1_20260806/01_run_e2_full_rank_gsea_v0.1_20260806.R
SHA256 f7f3a95dac41944b1fe52831e37110e3b4faccc25efffb5c622f907c38fe8dc4  8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/03_e2_full_rank_gsea/v0.1_20260806/01_run_e2_full_rank_gsea_v0.1_20260806.R
BLOCK_PASS 10_run_e2_full_rank_gsea.R source=10-126 public=72-188
BLOCK_PASS 10_run_e2_full_rank_gsea.R source=156-360 public=208-412
SOURCE 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/04_three_boundary_full_rank_gsea/v0.1_20260806/01_run_single_edge_full_rank_gsea_v0.1_20260806.R
SHA256 342ef15faaeb59d0b9627c258592d38e4eca0ad613dcb7023cfdd99f58e6224b  8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/04_three_boundary_full_rank_gsea/v0.1_20260806/01_run_single_edge_full_rank_gsea_v0.1_20260806.R
BLOCK_PASS 11_run_e1_e3_full_rank_gsea.R source=10-131 public=73-194
BLOCK_PASS 11_run_e1_e3_full_rank_gsea.R source=177-384 public=226-433
SOURCE 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/04_three_boundary_full_rank_gsea/v0.1_20260806/02_combine_three_edge_gsea_v0.1_20260806.R
SHA256 151bb82ae912b8c1f7da276493b09b11e00e233e5920c1058d31037cb97b3966  8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/exploratory/02_fig3_mrna_boundary_pairwise_ora/v0.1_20260806/04_three_boundary_full_rank_gsea/v0.1_20260806/02_combine_three_edge_gsea_v0.1_20260806.R
BLOCK_PASS 12_combine_three_edge_gsea.R source=7-26 public=69-88
BLOCK_PASS 12_combine_three_edge_gsea.R source=59-160 public=108-209
BLOCK_PASS 12_combine_three_edge_gsea.R source=163-185 public=211-233
BLOCK_PASS 12_combine_three_edge_gsea.R source=190-242 public=237-289
TOTAL_RETAINED_BLOCKS 11

Verification: four parse, four help and 41 synthetic metadata/path-safety cases passed; eleven contiguous retained-source blocks passed. Logs: interface_checks.log and static_checks.log. Source-only hashes; no research data/manifest/resource was read or hashed and no analysis or existing entry ran. Actual schemas, dependency completeness and numerical/end-to-end equivalence remain unverified.

## 13_run_mrna_rppa_concordance.R

Source: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/48_run_stage3_mrna_rppa_c1c4_audit_v0.1_20260730.R (737 lines). SHA-256: bf13f77e216a8318607c3f4c4db5d19a10ae917f85a2cd97cd34c48ba4cf3c6b. Contract: docs/stage3_mrna_rppa_c1c4_concordance_audit_contract_v0.1_20260730.md in the same source module. The complete expected_cancers-through-EOF body and log_event/collapse_values helper block are retained verbatim.

Finite adaptations: strict nine-field DCF help/metadata early exit, regular-file/15-model-family validation, symlink-resolved input/output isolation, independent empty output_dir, public script identity, and eight explicit input paths. Ordinary execution creates output_dir and uses 01_mrna_rppa_concordance as the source output directory, preserving the original temporary sibling and atomic rename. The input manifest records resolved supplied input paths, so it contains no temporary-output relocation. No model, transform, mapping, significance, recurrence, classification, QA assertion, output basename or serialization logic changes.

The 15-cancer within-layer pairing scope does not reduce or recompute supplied mRNA recurrence's upstream 22-cancer/evaluable denominator. Full-258 routes, unique TOTAL mapping, complete 15 x gene cells, all six per-cancer and gene-support classes, definition/display annotations, phospho/cleaved context, nine tables and original QA/input manifest remain. Historical R4.5.1 logging and identity assertions are retained without claiming fresh historical-environment or study-QA validation.

Verification on 2026-09-30: one parse, one help and 23 synthetic invalid-schema metadata/path-safety cases passed in R4.5.1. Two contiguous source-block comparisons passed, including the entire expected_cancers-through-EOF statistical/output/error-handling body. Source code alone was hashed. Logs: interface_checks.log and static_checks.log (maintainer-held artifacts). No research contents/schema/checksum, dependency loading, statistical body or plotting was executed; existing entries were not retested. Real schemas, dependencies and numerical/end-to-end reproduction remain unverified.

13_run_mrna_rppa_concordance.R source lines 91-737 = public lines 113-759: PASS.
13_run_mrna_rppa_concordance.R source lines 11-32 = public lines 75-96: PASS.

## 14_prepare_subtype_context.R and 15_analyze_subtype_context.R

Source: 13.0c1c4_biological_context_and_icb/scripts/01_prepare_A.R. SHA-256: ebf5b319730860cc622e056b7de74d206de7998c5085d2cc15013cefba33b959; 102 lines.

Source: 13.0c1c4_biological_context_and_icb/scripts/02_analyze_A.R. SHA-256: 2c62ce9be3c332e2995bc1e9d3233b71bfb22e53671054b420fab17eabf6d0b7; 94 lines.

Finite adaptations: prepend strict base-R CLI/DCF metadata checks; replace nine preparation inputs and two independent analysis inputs; change the log script identifier; move unchanged source basenames into separate temporary stages with atomic commit. All mapping, member selection, scoring, model, eligibility, BH, descriptive and sensitivity calculations are retained. Source scripts and study outputs are unchanged. Source module README and the S7 freeze contract supersede the design plan's historical unexecuted status. The result mapping is Supplementary Fig. 7 D–F and the complete Supplementary Table 3 increment; the mixed A/B candidate plot source is excluded. Static source-block comparisons and current metadata-interface checks do not establish numerical reproduction.

14–15 verification on 2026-09-30: two R4.5.1 parses, two help checks, 30 synthetic invalid-schema metadata/path-safety cases, two complete reversed-body comparisons and five contiguous unchanged source-block comparisons passed. Metadata success created no output directories; the nonempty-output sentinel was preserved. Checks covered strict DCF fields/duplicates, relative paths and reordered fields, regular-file types, missing/empty inputs, direct and symlink overlap, dangling symbolic links, root/config overlap, nonempty outputs and non-directory ancestors. Interface and static logs are maintainer-held task artifacts. No study contents/schema or checksum, analysis package, programme scoring/model/BH/descriptive statistical body or figure was loaded/executed. Existing entries were not retested. Actual schemas, dependency completeness, numerical equivalence and end-to-end reproduction remain unverified.


## 16_prepare_mrna_representatives.R

Sources: 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/scripts/37_freeze_stage3_mrna_c1c4_extreme_display_v0.1_20260730.R (SHA-256 1c22f3c87fa1aa0113ad9fd2faff89041c55d0b9b8a61946de07c6864cc7d55c); scripts/39_refine_stage3_mrna_c1c4_other_legacy_cytoscape_interface_v0.2_20260730.R in the same module (SHA-256 8e334026d53498d979eb94e06ce3871bd55d14322c09996a2bbe2a4668436bbc). The current source module README and docs/stage3_mrna_c1c4_extreme_ora_display_contract_v0.2_20260730.md establish representative input authority.

Contiguous verbatim blocks: source 37 lines 118–262 fixed cluster/reg_term/33-term registry; 373–410 required schemas, upstream-QA guard and candidate qualification/counts; 500–677 representative matches, candidate evidence order, term/cluster/C4 seeds, fill, display schema and directional freeze; source 39 lines 72–77 first_token and 203–267 complete display copy/C4 merge/Other/final ordering/primary-cluster body. Five source-block comparisons passed on the final script.

Finite adaptations: strict five-field DCF CLI/help/metadata safety; four explicit public-02 product paths; prepare_v1_display/refine_v2_display pure wrappers; nonempty complete upstream-QA requirement; independent temporary output and atomic commit; public script identity. Only necessary two v0.1 display/directional tables are written using source fwrite quote=TRUE/na=NA and read through fread, preserving the original 37-to-39 TSV boundary; final v0.2 basename/full fields/serialization are retained. Representative-input-only QA, supplied-path manifest, brief source explanation and sessionInfo are generated.

Source 37's 998-term Jaccard/Louvain, keyword audit and plotting/network outputs and source 39's legacy-input/CSV/network blocks are omitted. Jaccard/Louvain populates selected-term audit fields only and does not enter display-gene candidate ranking, seed, fill or 39 display refinement; no omitted audit columns are invented. The fixed source registry and eligibility remain, not an automatic reselection. Relevant source-derived count, membership and definition-marker assertions remain; complete historical network-QA counts are not claimed. No transformation, statistical test, BH family, direction recurrence, evidence rank rule, target cutoff, Other rule, column selection or scientific logic changes.

Verification: R4.5.1 (2025-06-13) parse/help; 18 synthetic invalid-schema metadata cases; 33-term/10-cluster registry check; source selection seed/fill on synthetic 214-UP/6-DOWN candidates; original TSV serialization/readback with source refinement into 230 synthetic records; C4 merge, ten Other records, first_token and required six consumer fields; five contiguous source-block comparisons; public-script private-root absence. Metadata successes created no output directory and preserved the nonempty-output sentinel. Final logs: tmp/mrna_representatives_final_20261003_5r8j6ce2/interface_checks.log and static_checks.log (internal maintainer artifacts, not runtime inputs). Only source code was hashed. No study input/schema/checksum, model/ORA/BH/statistical main, existing entry or figure was loaded or executed. Actual schemas, numerical equivalence, full dependency completeness and end-to-end reproduction remain unverified.

## Locally supplied bulk and programme references

Eight adopted references are locally packaged unchanged and configured through ../references/. Source R code and the differing historical TCPA mappings are preserved. The [reference interfaces](../../references/README.md) distinguish Fig. 3D post-statistical annotation from six-member drug RTK scoring, and identify a possible BH-mediated effect on Fig. 5C's first five RPPA rows; A 2026-10-06 R 4.5.1 local sensitivity comparison added HER3_pY1289 with positive sign to the historical six-member RTK score. The historical baseline passed numerical checks. Seven RTK tests and one non-RTK test crossed FDR < 0.05; RTK significant tests changed from 13 to 10 and non-RTK from 136 to 135. Fig. 5C's first five rows retained all 100 evaluable rho values, significance markers and significant-cancer counts (11/9/8/5/5), although FDR values changed. RPPA orthogonal-candidate membership was unchanged. This is a conditional mapping scenario, not an adopted correction or full package reproduction; TrackB/MOFA effects were not evaluated and formal assets remain unchanged. A 2026-10-06 read-only check verified exclusion of the sole NA from pathway inputs for the 120 nonempty formal E1/E2/E3 ranks across both ranking scopes. This closes that historical input concern, not future-input validation or numerical reproduction. The gene/cluster producer boundary is unchanged. Other prepared inputs and external permissions remain independent prerequisites.
