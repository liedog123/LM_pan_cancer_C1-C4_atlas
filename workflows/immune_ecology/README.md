# Bulk immune and stromal ecology

These two entries continue Stage 2 from prepared patient-level inputs. 01 retains the complete cancer-level statistical family and cross-cancer summaries; 02 retains the sensitivity-aware evidence classification and cross-method convergence actually consumed by the Fig. 2 immune panels. Neither entry reruns deconvolution or establishes numerical reproduction of frozen results.

## Invocation and execution order

Run from the package root after replacing all example placeholders:

    Rscript workflows/immune_ecology/01_run_stage2_immune_ecology.R --config config/immune_ecology.paths.example.dcf --check-config
    Rscript workflows/immune_ecology/01_run_stage2_immune_ecology.R --config /path/to/core.dcf
    Rscript workflows/immune_ecology/02_finalize_stage2_evidence.R --config config/immune_evidence.paths.example.dcf --check-config
    Rscript workflows/immune_ecology/02_finalize_stage2_evidence.R --config /path/to/evidence.dcf

01 writes results/ below its independent output_dir. Set 02 input_results to that results/ directory (or separately obtained compatible formal Stage 2 core tables), and choose a different independent output_dir for 02. Downstream consumers must take core tables from 01 and evidence/convergence tables from 02; no merged copy is made automatically.

Both entries accept only --config PATH with optional --check-config, or --help. The strict DCF contains one record and exactly the listed fields. Relative paths resolve from the DCF directory through existing ancestors and symlinks. Inputs must be readable nonempty regular files; output must be absent or empty, writable and independent of input parents and configuration. Root/home targets, nonempty output, overlapping output and symlink overlaps are rejected. Metadata checks read the DCF and file metadata only, load no analysis dependencies or research contents, and create no output directories. Normal invocations run the statistical or evidence/plotting body. The public entries fix smoke_mode=FALSE and selected_cancers=character() for the complete source range, and do not expose source smoke/cancer-subset CLI switches.

## Core inputs

| 01 field | Required content |
|---|---|
| canonical_labels | CSV with unique Sample, Cancer, C_subtype and AGIG. Sample keys are normalized to the first 15 TCGA barcode characters. |
| patient_manifest | TSV with unique Patient_ID and Sample, Cancer, C_subtype, AGIG, Acid_z, IG_z, discordant_multi_sample_state, age10, sex and purity_z. |
| estimate_scores | CSV with Sample, Cancer, StromalScore, ImmuneScore, ESTIMATEScore and TumorPurity. |
| cibersort_scores | CSV with Sample, Cancer, all 22 LM22 features, Absolute_score_(sig_score)_CIBERSORT, P-value, Correlation and RMSE. |
| mcpcounter_scores | CSV with Sample, Cancer and the ten source MCPcounter features. |
| acid_genes | Rdata containing final_selected_genes with Gene; used only for definition overlap audit. |
| lm12_genes | Long TSV containing one of gene, Gene, symbol or Symbol; used only for definition overlap audit. |
| output_dir | Independent absent/empty directory for new results/. |

The patient manifest already embodies deterministic representative-sample selection. The existing clinical_outcomes entry exports stage1_patient_analysis_manifest_v0.1_20260722.tsv with these transformed covariates; source code establishes this schema connection, but no real-input schema or numerical agreement was tested here. Frozen canonical labels remain the formal reference. ESTIMATE/MCPcounter old Subtype columns are dropped and scores rejoined to canonical patient keys; source Cancer is checked. CIBERSORT QC fields are not biological endpoints. Fixed feature names, normalization, join, transformation and failure rules are preserved in the script.

## Retained statistics and evidence

01 covers 4 ESTIMATE, 23 CIBERSORT absolute (22 cells plus total abundance) and 10 MCPcounter features. Designs are C1-C4 omnibus, three AGIG atomic edges, C1 versus C4, continuous Acid and continuous IG. Representative-main and EXCLUDE_DISCORDANT sets are both retained; CIBERSORT also retains ALL_MATCHED and CONFIDENT_P_LT_0.05 sets and state-dependent confidence-retention audits.

Source Kruskal-Wallis/epsilon-squared, Wilcoxon/median difference/rank-biserial/Hedges' g, Spearman/Fisher z, rank inverse-normal transformation and age/sex or age/sex/purity HC3 models are retained. ESTIMATE does not adjust its derived scores for purity. All sample/complete-rate/group-size gates, numerical HC3 downgrades, failure and non-estimable rows, feature-family BH, custom REML, heterogeneity, prediction intervals, direction consistency, leave-one-cancer-out and cancer-equal-weight state summaries are retained. Adjusted coefficient BH additionally separates Effect_Level as in the source; global adjusted BH and meta BH remain distinct. No displayed-cell, C1/C4-only or significant-cancer subset redefines a family.

02 combines primary meta results with age/sex, purity and CIBERSORT-confidence sensitivity, applies the source FDR/direction/I2/LOO evidence rules, and records cross-method directional convergence and cancer-specific evidence. All original classifications and negative/context-dependent outcomes are retained. Its original Gate 2 table, report, source heatmaps, convergence plot and confidence-retention plot remain source diagnostics and exports. Their fixed historical date, version suffix and completion wording are copied source conventions, not evidence of a newly passed Gate or reproduction of current frozen figures. They are generated only in a normal run, which was not performed during packaging.

01 writes 14 source tables plus stage2_core_output_manifest. 02 requires exactly eleven version-suffixed TSV basenames: immune_random_effects_meta, immune_leave_one_cancer_out, immune_pairwise_effects_by_cancer, immune_continuous_axis_associations_by_cancer, immune_omnibus_by_cancer, immune_adjusted_models_by_cancer, immune_input_and_feature_registry, cibersort_confidence_retention_audit, immune_state_profile_equal_weight, stage2_sample_analysis_manifest and immune_coverage_and_evaluability_by_cancer, each ending _v0.1_20260723.tsv. 02 writes evidence, cancer-specific evidence, convergence, Gate, figure manifest and output manifest in its new results/, source PDF/TIFF/TSV exports in figures/, and its source diagnostic report in docs/. Its output manifest lists new finalization artifacts only, not copies of its input tables.

## Sources, figure mapping and checks

See [source mapping](source_mapping.md) for exact source paths, fingerprints, body adaptations and formal panel consumers. Fig. 2E and Supplementary Fig. 4A use core pairwise statistics and final evidence classifications. Supplementary Fig. 4B uses the patient manifest, reused ESTIMATE scores and Stage 2 statistics; its current frozen panel removes the cross-cancer meta display while retaining its statistical background. The two entries provide the full analytical family, not the current frozen panel layout or final supplementary-table exporter. Source plots do not reproduce the canonical palette/assembly.

R4.5.1 checks cover parsing, synthetic nonempty-file metadata CLI/safety cases and source text comparison only. Synthetic files intentionally contain invalid research schemas. No study table, Rdata, score, gene membership or research-input checksum was read; no model, REML/LOO, evidence classification, deconvolution or plot was run. Real schemas, dependency/device completeness, numerical reproduction and end-to-end equivalence remain unverified.

## Input access and release boundary

No patient data, scores, signatures or third-party resources are distributed. Canonical/definition copies can be generated by the existing shared_definitions/signature/LM12 workflows only after obtaining their upstream inputs; equivalence to the study freeze is unverified. The patient manifest can continue through clinical_outcomes from its separately obtained readiness and clinical inputs. The three deconvolution score tables currently require an author-approved intermediate-data release or controlled handoff with sample keys and algorithm/version documentation; a public download/accession is not yet supplied. Generating new tables through IOBR requires compatible expression and the source method interfaces (cibersort_abs with abs.method="sig.score", estimate, mcpcounter), third-party terms and a separately authorized upstream implementation; these entries do not implement or test that route. IOBR signature resources are loaded from the installed package for overlap audit, not redistributed. This unresolved input access/permission boundary prevents claiming a complete independently runnable release. Code/data licensing and publication remain pending; no license is granted by this workspace.
